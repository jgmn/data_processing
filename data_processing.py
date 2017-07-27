# -*- coding: utf-8 -*-
"""
@author: J. Gerardo Moreno N.
"""
import json
import numpy as np
import geopandas as gpd
import pandas as pd
import warnings
from folium import Map, GeoJson, LayerControl
from time import time 
from os import listdir
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point

#------------------------VORONOI_FINITE_POLYGONS------------------------------#
# DESCRIPCIÓN: Reconstruye las regiones infinitas del diagrama de Voronoi     #
#              a regiones finitas.                                            #
# PARÁMETROS:                                                                 #
#   ENTRADA: vor: Diagrama de Voronoi en 2D.                                  #
#            radius: Distancia para los puntos infinitos (opcional).          #
#   SALIDA:  regions: Indice de vértices en cada región de Voronoi.           #
#            vertices: Coordenadas para cada vértice de Voronoi.              #
#-----------------------------------------------------------------------------#
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
#-----------------------------------------------------------------------------#

#------------------------CLIP_VORONOI-----------------------------------------#
# DESCRIPCIÓN: Recorta el Voronoi con respecto a una caja.                    # 
# PARÁMETROS:                                                                 #
#   ENTRADA: regions: Indice de vértices en cada región de Voronoi.           #
#            vertices: Coordenadas para cada vértice de Voronoi.              #
#            box: Caja para recortar el Voronoi.                              #
#   SALIDA:  voro_df: GeoDataFrame Voronoi.                                   #
#-----------------------------------------------------------------------------#
def clip_voronoi(regions, vertices, box):
     poly_list = []
     for region in regions:
        polygon = vertices[region]
        poly = Polygon(polygon)
        poly = poly.intersection(box)
        polygon = [[p[0], p[1]] for p in poly.exterior.coords]
        poly_list.append(Polygon(polygon))

     voro_df = gpd.GeoDataFrame(gpd.GeoSeries(poly_list), columns = ['geometry'])
     voro_df.crs = {'init' :'epsg:4326'}
     
     return voro_df
#-----------------------------------------------------------------------------#

#---------------------------ADD_POBL------------------------------------------#
# DESCRIPCIÓN: Agrega la población total a los puntos de interés.             # 
# PARÁMETROS:                                                                 #
#   ENTRADA: pobl_df: GeoDataFrame Población.                                 #
#            voro_df: GeoDataFrame Voronoi.                                   #
#            pdi_df: GeoDataFrame Puntos de Interés.                          #
#   SALIDA:  cali: GeoJson Calificaciones.                                    #
#-----------------------------------------------------------------------------#
def add_pobl(pobl_df, voro_df, pdi_df):
    pobl_df = pobl_df.to_crs({'init': 'epsg:4326'})
    pobl_df.crs = {'init' :'epsg:4326'}
    voro_df.crs = pobl_df.crs # Para quitar el warning CRS does not match !
    pobl_with_voro = gpd.sjoin(pobl_df, voro_df, how = "inner", op = 'intersects')
    pdi_df['poblacion'] = 0
    
    for data in pobl_with_voro.iterrows():
        pdi_df['poblacion'].iloc[[data[1][6]]] = pdi_df['poblacion'].iloc[[data[1][6]]] + int(data[1][4])
        
    return pdi_df    
#-----------------------------------------------------------------------------#

#---------------------------ADD_TWEETS----------------------------------------#
# DESCRIPCIÓN: Agrega la actividad de redes sociales a los puntos de interés  # 
# PARÁMETROS:                                                                 #
#   ENTRADA: tweets_df: GeoDataFrame Tweets.                                  #
#            voro_df: GeoDataFrame Voronoi.                                   #
#            pdi_df: GeoDataFrame Puntos de Interés.                          #
#   SALIDA:  pdi_df: GeoDataFrame Puntos de Interés.                          #
#-----------------------------------------------------------------------------#               
def add_tweets(tweets_df, voro_df, pdi_df):    
    voro_df.crs = tweets_df.crs # Para quitar el warning CRS does not match !
    tweets_with_voro = gpd.sjoin(tweets_df, voro_df, how = "inner", op = 'intersects')
    pdi_df['tweets'] = 0
    
    for data in tweets_with_voro.iterrows():
        pdi_df['tweets'].iloc[[data[1][1]]] = pdi_df['tweets'].iloc[[data[1][1]]] + 1
    
    return pdi_df
#-----------------------------------------------------------------------------#

#---------------------------ADD_TRAFFIC---------------------------------------#
# DESCRIPCIÓN: Agrega el tráfico a los puntos de interés.                     # 
# PARÁMETROS:                                                                 #
#   ENTRADA: traficos_df: Lista GeoDataFrames Tráfico.                        # 
#            voro_df: GeoDataFrame Voronoi.                                   #
#            pdi_df: GeoDataFrame Puntos de Interés.                          #
#   SALIDA: pdi_df: GeoDataFrame Puntos de Interés.                           #
#-----------------------------------------------------------------------------#
def add_traffic(traficos_df, voro_df, pdi_df):
    for trafico_df in traficos_df:
        trafico_df = trafico_df.to_crs({'init': 'epsg:4326'})
        trafico_df.crs = {'init' :'epsg:4326'}
        voro_df.crs = trafico_df.crs # Para quitar el warning CRS does not match !
        trafico_with_voro = gpd.sjoin(trafico_df, voro_df, how = "inner", op = 'intersects')
        pdi_df['trafico'] = 0
        
        for data in trafico_with_voro.iterrows():
            pdi_df['trafico'].iloc[[data[1][4]]] = pdi_df['trafico'].iloc[[data[1][4]]] + int(data[1][3])  
    
    return pdi_df                                                                                          
#-----------------------------------------------------------------------------#

#-------------------------DEL_CLUSTER_POINTS----------------------------------#
# DESCRIPCIÓN: Elimina los clusteres de puntos en los puntos de interés.      # 
# PARÁMETROS:                                                                 #
#   ENTRADA: pdi_df: GeoDataFrame Puntos de Interés.                          #
#   SALIDA:  pdi_df: GeoDataFrame Puntos de Interés sin clusteres.            #
#            coords: Coordenadas de los puntos de interés.                    #
#-----------------------------------------------------------------------------#
def del_cluster_points(pdi_df):
    coords = []
    pdi_df = pdi_df.to_crs({'init': 'epsg:25830'})
    for x in pdi_df.iterrows():
        dist = []
        indices = []
        for y in pdi_df.iterrows():
            p1 = x[1][1]
            p2 = y[1][1]
            d = p1.distance(p2)
            dist.append(d)
        pdi_df_copy = pdi_df.copy()
        pdi_df_copy['distancias'] = dist
        pdi_df_subset = pdi_df_copy[pdi_df_copy['distancias'] <= 150]
        if not pdi_df_subset.empty:
            pdi_df_subset
            pdi_df_subset['total'] = pdi_df_subset.sum(axis = 1)
            max_index = pdi_df_subset['total'].idxmax()
            pdi_df_subset.drop(max_index, inplace = True)
            indices = list(pdi_df_subset.index.values)
            if indices:
                pdi_df.drop(indices, inplace = True)
        del pdi_df_copy
        del pdi_df_subset 
    
    pdi_df = pdi_df.to_crs({'init': 'epsg:4326'})
    coords = ([[r[1][1].x, r[1][1].y] for r in pdi_df.iterrows()])
    
    return pdi_df, coords
#-----------------------------------------------------------------------------#

def main():
    # Paso 1: Eliminar clúster de puntos en base al tiempo mediovque pasan las 
    # personas en los puntos de interés (PDI). El archivo tiempo.json está 
    # filtrado por los PDI y contiene el nombre y tiempo de los puntos.
    print("Leyendo archivo tiempo...")
    pdi_df = gpd.read_file('opendata/tiempo.json')
    pdi_df = pdi_df.sort_index()    
    print("Eliminando clúster de puntos...")
    pdi_df, coords = del_cluster_points(pdi_df)
    
    # Paso 2: Generar el diagrama de Voronoi para el procesamiento de datos.
    vor = Voronoi(coords)
    regions, vertices = voronoi_finite_polygons(vor)
    min_x = vor.min_bound[0] - 0.1
    max_x = vor.max_bound[0] + 0.1
    min_y = vor.min_bound[1] - 0.1
    max_y = vor.max_bound[1] + 0.1
    box = Polygon([[min_x, min_y], [min_x, max_y], [max_x, max_y], [max_x, min_y]])
    voro_df = clip_voronoi(regions, vertices, box)
    
    # Paso 3: Agregar los datos de población, tráfico y tweets.
    # Población
    print('Leyendo archivo manzanas_pob...')
    pobl_df = gpd.read_file('opendata/manzanas_pob.json')
      
    print('Agregando datos de población...')
    pdi_df = add_pobl(pobl_df, voro_df, pdi_df)
    
    # Tráfico
    print('Leyendo archivos de tráfico...')
    traficos_df = []
    for filename in listdir('opendata/trafico'):
        traficos_df.append(gpd.read_file('opendata/trafico/'+filename))
    
    print('Agregando datos de tráfico...')
    pdi_df = add_traffic(traficos_df, voro_df, pdi_df)
    
    # Tweets 
    print('Leyendo archivo valencia_tweets...')
    with open('opendata/valencia_tweets.json','r') as input_file:
        tweets = json.load(input_file)
    
    tweets_df = gpd.GeoDataFrame(gpd.GeoSeries(Point(t['coordinates']['coordinates']) for t in tweets), columns=['geometry'])
          
    print('Agregando datos de redes sociales...')
    pdi_df = add_tweets(tweets_df, voro_df, pdi_df)
    
    # Paso 4: Guardar los archivos para el funcionamiento del algoritmo genético.
    print('Guardando puntos de interés...')
    pdi_df.to_file('puntos_de_interes.json', driver = "GeoJSON")
    
    print('Guardando voronoi...')
    voro_df.to_file('voronoi.json', driver = "GeoJSON")
     
    # Paso 5: Visualizar la información en mapa web.
    print('Generando mapa web...')
    valencia = [39.4561165311493, -0.3545661635]
    mapa = Map(location = valencia, tiles = 'OpenStreetMap', zoom_start = 10)
    GeoJson(open('voronoi.json'), name = 'Diagrama de Voronoi').add_to(mapa)
    GeoJson(open('puntos_de_interes.json'), name = 'Puntos de Interés').add_to(mapa)
    LayerControl().add_to(mapa)
    mapa.save('valencia.html')
    
    print('Listo')

if __name__ == "__main__":
    pd.options.mode.chained_assignment = None # Para quitar el warning SettingWithCopyWarning
    warnings.filterwarnings("ignore", category = RuntimeWarning) # Para quitar el warning RuntimeWarning
    # Calcular tiempo de ejecución
    tiempo_inicial = time()
    main()
    tiempo_final = time()
    tiempo_ejecucion = tiempo_final - tiempo_inicial
    print('Tiempo de ejecución: ', '%.2f'% (tiempo_ejecucion/60), 'minutos')