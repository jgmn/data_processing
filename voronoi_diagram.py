import json
import folium
import numpy as np
# import geopandas as gpd
from os import listdir
from pyproj import Proj, transform
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point, LineString
from geojson import FeatureCollection, Feature, Polygon as PolyJson

def voronoi_finite_polygons(vor, radius = None):
#------------------------VORONOI_FINITE_POLYGONS--------------------------#
# DESCRIPCIÓN: Reconstruye las regiones infinitas del diagrama de Voronoi #
#              a regiones finitas.                                        #
# PARÁMETROS:                                                             #
#   ENTRADA: vor: Diagrama de Voronoi en 2D.                              #
#            radius: Distancia para los puntos infinitos (opcional).      #
#   SALIDA:  regions: Indice de vértices en cada región de Voronoi.       #
#            vertices: Coordenadas para cada vértice de Voronoi.          #
#-------------------------------------------------------------------------#
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

def get_coords(geojson):
#------------------------GET_COORDS---------------------------------------#
# DESCRIPCIÓN: Extrae las coordenadas de un archivo GeoJSON Points.       #
# PARÁMETROS:                                                             #
#   ENTRADA: geojson: Archivo GeoJSON.                                    #
#   SALIDA:  coords: Coordenadas del archivo GeoJSON.                     #
#-------------------------------------------------------------------------#
    coords = []
    for feature in geojson['features']:
        coordinate = feature['geometry']['coordinates']
        x, y = coordinate[1], coordinate[0]
        coords.append([x,y])
    return coords
#-------------------------------------------------------------------------#

def clip_voronoi(regions, vertices, box):
#------------------------CLIP_VORONOI-------------------------------------#
# DESCRIPCIÓN: Recorta el Voronoi con respecto a una caja y genera un     #
#              archivo Voronoi con los polígonos recortados.              # 
# PARÁMETROS:                                                             #
#   ENTRADA: regions: Indice de vértices en cada región de Voronoi.       #
#            vertices: Coordenadas para cada vértice de Voronoi.          #
#            box: Caja para recortar el Voronoi.                          #
#-------------------------------------------------------------------------#
    vorJSON = open('voronoi.JSON', 'w')
    feature_list = []
    for region in regions:
        vertice_list = []
        polygon = vertices[region]
        poly = Polygon(polygon)
        poly = poly.intersection(box)
        polygon = [[p[1], p[0]] for p in poly.exterior.coords]
        vertice_list.append(polygon)
        temp = PolyJson(vertice_list)
        feature = Feature(geometry=temp, properties={})
        feature_list.append(feature)

    feature_collection = FeatureCollection(feature_list)
    print (feature_collection, file=vorJSON)
    vorJSON.close()
#-------------------------------------------------------------------------#

def add_pobl(voro, pobl, cali):
#---------------------------ADD_POBL--------------------------------------#
# DESCRIPCIÓN: Agrega la población total a las calificaciones.            # 
# PARÁMETROS:                                                             #
#   ENTRADA: voro: GeoJSON Voronoi.                                       #
#            pobl: GeoJSON Población.                                     #
#            cali: GeoJSON Calificaciones.                                #
#-------------------------------------------------------------------------#
    for pobl_feature in pobl['features']:
        p1 = Polygon(pobl_feature['geometry']['coordinates'][0])
        for index, voro_feature in enumerate(voro['features']):
            p2 = Polygon(voro_feature['geometry']['coordinates'][0])
            if(p1.intersects(p2)):
                cali['features'][index]['properties']['poblacion'] = cali['features'][index]['properties']['poblacion'] + int(pobl_feature['properties']['pob_total'])
    return cali
#--------------------------------------------------------------------------#
               
def add_tweets(voro, tweets, cali):
#---------------------------ADD_TWEETS-------------------------------------#
# DESCRIPCIÓN: Agrega la actividad de redes sociales a las calificaciones. # 
# PARÁMETROS:                                                              #
#   ENTRADA: voro:   GeoJSON Voronoi.                                      #
#            tweets: GeoJSON Tweets.                                       #
#            cali:   GeoJSON Calificaciones.                               #
#--------------------------------------------------------------------------#
    for tweet in tweets:
        p1 = Point(tweet['coordinates']['coordinates'])
        for index, voro_feature in enumerate(voro['features']):
            p2 = Polygon(voro_feature['geometry']['coordinates'][0])
            if(p1.intersects(p2)):
                cali['features'][index]['properties']['tweets'] = cali['features'][index]['properties']['tweets'] + 1
    return cali
#--------------------------------------------------------------------------#

def add_traffic(voro, cali):
#---------------------------ADD_TRAFFIC------------------------------------#
# DESCRIPCIÓN: Agrega el tráfico a las calificaciones.                     # 
# PARÁMETROS:                                                              #
#   ENTRADA: voro: GeoJSON Voronoi.                                        #
#            cali: GeoJSON Calificaciones.                                 #
#--------------------------------------------------------------------------#
    path_folder = 'traffic_data'
    for file in listdir(path_folder):
        with open(path_folder+'/'+file,'r', encoding="utf8") as input_file:
            data = json.load(input_file)

        traffic = convert_coords(data)

        for traffic_feature in traffic['features']:
            l1 = LineString(traffic_feature['geometry']['coordinates'])
            for index, voro_feature in enumerate(voro['features']):
                p2 = Polygon(voro_feature['geometry']['coordinates'][0])
                if(l1.intersects(p2)):
                    cali['features'][index]['properties']['trafico'] = cali['features'][index]['properties']['trafico'] + int(traffic_feature['properties']['lectura'])
    return cali                                                                                               
#--------------------------------------------------------------------------#

#---------------------------CONVERT_COORDS---------------------------------#
# DESCRIPCIÓN: Convierte las coordenadas de un archivo GeoJSON.            # 
# PARÁMETROS:                                                              #
#   ENTRADA: data:   GeoJSON que se desea convertir.                       #
#   SALIDA:  result: GeoJSON con las nuevas coordenadas.                   #
#--------------------------------------------------------------------------#
def convert_coords(data):
    inpProj = Proj(init='epsg:25830')
    outProj = Proj(init='epsg:4326')

    result = {}
    result['type'] = data['type']
    result['crs'] = data['crs']
    result['features'] = []

    for feature in data['features']:
        coordinates = []
        for coordinate in feature['geometry']['coordinates']:
            x1, y1 = coordinate[0], coordinate[1]
            x2, y2 = transform(inpProj, outProj, x1, y1)
            coordinates.append([x2, y2])
        feature['geometry']['coordinates'] = coordinates
        result['features'].append(feature)

    result['crs']['properties']['name'] = "urn:ogc:def:crs:EPSG::4326"
    return result
#--------------------------------------------------------------------------#

# MAIN
print('Creando el mapa...')
valencia = [39.4561165311493, -0.3545661635]
mapVor = folium.Map(location = valencia, zoom_start = 13)

# Puntos de interés
print('Leyendo archivo calificaciones...')
path_file = 'calificaciones.JSON'
with open(path_file, 'r') as input_file:
    cali = json.load(input_file)

print('Extrayendo puntos de interés...')
coords = get_coords(cali)

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

clip_voronoi(regions, vertices, box)

# Población
print('Leyendo archivo manzanas_pob...')
path_file = 'manzanas_pob.JSON'
with open(path_file,'r') as input_file:
    pobl = json.load(input_file)

# Tweets 
print('Leyendo archivo valenciatweets...')
path_file = 'valenciatweets.JSON'
with open(path_file,'r') as input_file:
    tweets = json.load(input_file)
    
# Voronoi
print('Leyendo archivo voronoi...')
path_file = 'voronoi.JSON'
with open(path_file,'r') as input_file:
    voro = json.load(input_file)

# Procesamiento de datos
print('Agregando datos de población...')
cali = add_pobl(voro, pobl, cali)

print('Agregando datos de redes sociales...')
cali = add_tweets(voro, tweets, cali)

print('Agregando datos de tráfico...')
cali = add_traffic(voro, cali)

print('Guardando caliicaciones...')
with open('calificaciones.JSON', 'w') as output_file:
    json.dump(cali, output_file, indent=3)

print('Agregando Voronoi al mapa...')
folium.GeoJson(open('voronoi.JSON'), name='Diagrama de Voronoi').add_to(mapVor)
folium.LayerControl().add_to(mapVor)

print('Guardando mapa...')
mapVor.save('valencia.html')

print('Listo')
