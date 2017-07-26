# -*- coding: utf-8 -*-
"""
@author: J. Gerardo Moreno N.
"""
import json
import numpy as np
import geopandas as gpd
from folium import Map, GeoJson, LayerControl
from time import time 
from os import listdir
from pyproj import Proj, transform
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point
from geojson import FeatureCollection, Feature, Polygon as PolyJson

#-----------------------------FILTER_GEOJSON------------------------------#
# DESCRIPCIÓN: Filtra un GeoJSON por los puntos de interés, calcula el    #
#              centroide de los polígonos y cambia el formato de coords.  #
# PARÁMETROS:                                                             #
#   ENTRADA: data: GeoJSON que se desea filtrar.                          #
#   SALIDA:  result: GeoJSON filtrado.                                    #
#            coords: Coordenadas de result.                               #
#-------------------------------------------------------------------------#
def filter_geojson(cali):
    inpProj = Proj(init='epsg:25830')
    outProj = Proj(init='epsg:4326')
    
    coords = []
    
    result = {}
    result['type'] = cali['type']
    result['crs'] = cali['crs']
    result['features'] = []

    for feature in cali['features']:
        if ((feature['properties']['califi'] == "EDA" and  feature['properties']['uso'] == "SJL")
            or (feature['properties']['califi'] == "PID")
            or (feature['properties']['califi'] == "GTR"  and  feature['properties']['tipoca'] == "4")
            or (feature['properties']['califi'] == "GSP"  and (feature['properties']['tipoca'] == "4" or feature['properties']['tipoca'] == "3"))
            or (feature['properties']['califi'] == "TER"  and  feature['properties']['tipoca'] == "3")
            or (feature['properties']['califi'] == "E/SP" and  feature['properties']['tipoca'] == "1P")):
            
            feature['properties']['nombre'] = ""
            feature['properties']['tiempo'] = 0
            feature['properties']['poblacion'] = 0
            feature['properties']['tweets'] = 0
            feature['properties']['trafico'] = 0
            
            del feature['properties']['clase']
            del feature['properties']['tipouso']
            del feature['properties']['origen']
            del feature['properties']['ficha_es']
            del feature['properties']['ficha_va']
            
            polygon = Polygon(feature['geometry']['coordinates'][0])
            point = polygon.envelope.centroid 
            x1, y1 = point.x, point.y
            x2, y2 = transform(inpProj, outProj, x1, y1)
            coords.append([x2, y2])
            feature['geometry']['coordinates'] = [x2, y2]
            feature['geometry']['type'] = "Point"
            result['features'].append(feature)

    result['crs']['properties']['name'] = "urn:ogc:def:crs:EPSG::4326" 
    return result, coords
#-------------------------------------------------------------------------#

#------------------------VORONOI_FINITE_POLYGONS--------------------------#
# DESCRIPCIÓN: Reconstruye las regiones infinitas del diagrama de Voronoi #
#              a regiones finitas.                                        #
# PARÁMETROS:                                                             #
#   ENTRADA: vor: Diagrama de Voronoi en 2D.                              #
#            radius: Distancia para los puntos infinitos (opcional).      #
#   SALIDA:  regions: Indice de vértices en cada región de Voronoi.       #
#            vertices: Coordenadas para cada vértice de Voronoi.          #
#-------------------------------------------------------------------------#
def voronoi_finite_polygons(vor, radius = None):
    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()*2

    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            new_regions.append(vertices)
            continue

        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                continue

            t = vor.points[p2] - vor.points[p1] 
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)
#-------------------------------------------------------------------------#

#------------------------CLIP_VORONOI-------------------------------------#
# DESCRIPCIÓN: Recorta el Voronoi con respecto a una caja.                # 
# PARÁMETROS:                                                             #
#   ENTRADA: regions: Indice de vértices en cada región de Voronoi.       #
#            vertices: Coordenadas para cada vértice de Voronoi.          #
#            box: Caja para recortar el Voronoi.                          #
#   SALIDA:  voro: GeoJSON Voronoi.                                       #
#-------------------------------------------------------------------------#
def clip_voronoi(regions, vertices, box):
    crs = {
    "type": "name",
    "properties": {
        "name": "urn:ogc:def:crs:EPSG::4326"
        }
    }
    
    feature_list = []
    for region in regions:
        vertice_list = []
        polygon = vertices[region]
        poly = Polygon(polygon)
        poly = poly.intersection(box)
        polygon = [[p[0], p[1]] for p in poly.exterior.coords]
        vertice_list.append(polygon)
        temp = PolyJson(vertice_list)
        feature = Feature(geometry=temp, properties={})
        feature_list.append(feature)

    feature_collection = FeatureCollection(feature_list, crs=crs)
    return feature_collection
#-------------------------------------------------------------------------#

#---------------------------ADD_POBL--------------------------------------#
# DESCRIPCIÓN: Agrega la población total a las calificaciones.            # 
# PARÁMETROS:                                                             #
#   ENTRADA: voro: GeoDataFrame Voronoi.                                  #
#            pobl: GeoDataFrame Población.                                #
#            cali: GeoJson Calificaciones.                                #
#-------------------------------------------------------------------------#
def add_pobl(pobl_df, voro_df, cali):
    pobl_df = pobl_df.to_crs({'init': 'epsg:4326'})
    pobl_df.crs = {'init' :'epsg:4326'}
    voro_df.crs = pobl_df.crs # Para quitar el warning CRS does not match !
    pobl_with_voro = gpd.sjoin(pobl_df, voro_df, how="inner", op='intersects')
        
    for data in pobl_with_voro.iterrows():
        cali['features'][data[1][6]]['properties']['poblacion'] = cali['features'][data[1][6]]['properties']['poblacion'] + int(data[1][4])
  
    return cali    
#--------------------------------------------------------------------------#

#---------------------------ADD_TWEETS-------------------------------------#
# DESCRIPCIÓN: Agrega la actividad de redes sociales a las calificaciones. # 
# PARÁMETROS:                                                              #
#   ENTRADA: tweets: GeoDataFrame Tweets.                                  #
#            voro:   GeoDataFrame Voronoi.                                 #
#            cali:   GeoJSON Calificaciones.                               #
#--------------------------------------------------------------------------#               
def add_tweets(tweets_df, voro_df, cali):    
    voro_df.crs = tweets_df.crs # Para quitar el warning CRS does not match !
    tweets_with_voro = gpd.sjoin(tweets_df, voro_df, how="inner", op='intersects')
    
    for data in tweets_with_voro.iterrows():
        cali['features'][data[1][1]]['properties']['tweets'] = cali['features'][data[1][1]]['properties']['tweets'] + 1
    
    return cali
#--------------------------------------------------------------------------#

#---------------------------ADD_TRAFFIC------------------------------------#
# DESCRIPCIÓN: Agrega el tráfico a las calificaciones.                     # 
# PARÁMETROS:                                                              #
#   ENTRADA: traficos: Lista de GeoDataFrames de tráfico.                  # 
#            voro: GeoDataFrame Voronoi.                                   #
#            cali: GeoJSON Calificaciones.                                 #
#--------------------------------------------------------------------------#
def add_traffic(traficos_df, voro_df, cali):
    for trafico_df in traficos_df:
        trafico_df = trafico_df.to_crs({'init': 'epsg:4326'})
        trafico_df.crs = {'init' :'epsg:4326'}
        voro_df.crs = trafico_df.crs # Para quitar el warning CRS does not match !
        trafico_with_voro = gpd.sjoin(trafico_df, voro_df, how="inner", op='intersects')
        
        for data in trafico_with_voro.iterrows():
            cali['features'][data[1][4]]['properties']['trafico'] = cali['features'][data[1][4]]['properties']['trafico'] + int(data[1][3])    
    
    return cali                                                                                          
#--------------------------------------------------------------------------#

#-----------------------------ADD_TIME-------------------------------------#
# DESCRIPCIÓN: Agrega el tiempo de espera a las calificaciones.            # 
# PARÁMETROS:                                                              #
#   ENTRADA: tiempo: GeoJSON Tiempo.                                       #
#            cali:   GeoJSON Calificaciones.                               #
#--------------------------------------------------------------------------#
def add_time(tiempo, cali):
    for index, tiempo_feature in enumerate(tiempo['features']):
        nombre = tiempo_feature['properties']['nombre']
        tiempo = tiempo_feature['properties']['tiempo']
        cali['features'][index]['properties']['nombre'] = nombre
        cali['features'][index]['properties']['tiempo'] = tiempo
        
    return cali
#--------------------------------------------------------------------------#

def main():
    # Puntos de interés
    print('Leyendo archivo calificaciones...')
    with open('opendata/calificaciones.JSON', 'r') as input_file:
        cali = json.load(input_file)
    
    print('Filtrando calificaciones...')
    cali, coords = filter_geojson(cali)
    
    print('Creando Voronoi...')
    vor = Voronoi(coords)
    
    print('Reconstruyendo puntos infinitos...')
    regions, vertices = voronoi_finite_polygons(vor)
    
    print ('Recortando Voronoi...')
    min_x = vor.min_bound[0] - 0.1
    max_x = vor.max_bound[0] + 0.1
    min_y = vor.min_bound[1] - 0.1
    max_y = vor.max_bound[1] + 0.1
    
    box = Polygon([[min_x, min_y], [min_x, max_y], [max_x, max_y], [max_x, min_y]])
    
    voro = clip_voronoi(regions, vertices, box)
    
    voro_df = gpd.GeoDataFrame(gpd.GeoSeries(Polygon(v['geometry']['coordinates'][0]) for v in voro['features']), columns=['geometry'])
    
    # Población
    print('Leyendo archivo manzanas_pob...')
    pobl_df = gpd.read_file('opendata/manzanas_pob.json')
      
    print('Agregando datos de población...')
    cali = add_pobl(pobl_df, voro_df, cali)
    
    # Tiempo medio
    print('Leyendo archivo tiempo_medio...')
    with open('opendata/tiempo.json','r') as input_file:
        tiempo = json.load(input_file)
        
    print('Agregando datos de tiempo...')
    cali = add_time(tiempo, cali)
    
    # Tráfico
    print('Leyendo archivos de tráfico...')
    traficos_df = []
    for filename in listdir('opendata/trafico'):
        traficos_df.append(gpd.read_file('opendata/trafico/'+filename))
    
    print('Agregando datos de tráfico...')
    cali = add_traffic(traficos_df, voro_df, cali)
    
    # Tweets 
    print('Leyendo archivo valencia_tweets...')
    with open('opendata/valencia_tweets.json','r') as input_file:
        tweets = json.load(input_file)
    
    tweets_df = gpd.GeoDataFrame(gpd.GeoSeries(Point(t['coordinates']['coordinates']) for t in tweets), columns=['geometry'])
          
    print('Agregando datos de redes sociales...')
    cali = add_tweets(tweets_df, voro_df, cali)
    
    # Guardar los archivos de salida
    print('Guardando puntos de interés...')
    with open('puntos_de_interes.json', 'w') as output_file:
        json.dump(cali, output_file, indent=3)
    
    print('Guardando voronoi...')
    with open('voronoi.json', 'w') as output_file:
        json.dump(voro, output_file, indent=3)
      
    # Visualizar la información en mapa web
    print('Generando mapa web...')
    valencia = [39.4561165311493, -0.3545661635]
    mapa = Map(location = valencia, tiles = 'OpenStreetMap', zoom_start = 10)
    GeoJson(open('voronoi.json'), name = 'Diagrama de Voronoi').add_to(mapa)
    GeoJson(open('puntos_de_interes.json'), name = 'Puntos de Interés').add_to(mapa)
    LayerControl().add_to(mapa)
    mapa.save('valencia.html')
    
    print('Listo')

if __name__ == "__main__":
    # Calcular tiempo de ejecución
    tiempo_inicial = time()
    main()
    tiempo_final = time()
    
    tiempo_ejecucion = tiempo_final - tiempo_inicial
    print('Tiempo de ejecución: ', '%.2f'% (tiempo_ejecucion/60), 'minutos')