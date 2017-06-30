import json
import folium
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from shapely.geometry import Polygon
from geojson import FeatureCollection, Feature, Polygon as polyjson

def voronoi_finite_polygons_2d(vor, radius=None):
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

def get_coords(cali):
    coords = []
    for feature in cali['features']:
        #description = "califi: "+ feature['properties']['califi'] + " " + "tipoca: "+ feature['properties']['tipoca']
        coordinate = feature['geometry']['coordinates']
        x, y = coordinate[1], coordinate[0]
        #folium.Marker([x,y], popup = description).add_to(mapVor) 
        coords.append([x,y])

    return coords

def clip_voronoi(regions, vertices, box):
    feature_list = []
    for region in regions:
        temp = []
        polygon = vertices[region]
        poly = Polygon(polygon)
        poly = poly.intersection(box)
        polygon = [(p[1], p[0]) for p in poly.exterior.coords]
        temp = polyjson([polygon])
        feature = Feature(geometry=temp, properties={})
        feature_list.append(feature)

    return feature_list

# Main
print('Creando el mapa...')
valencia = [39.4561165311493, -0.3545661635]
mapVor = folium.Map(location = valencia, zoom_start = 13)

print('Leyendo archivo calificaciones...')
path_input_file = 'calificaciones.JSON'
with open(path_input_file,"r") as input_file:
    cali = json.load(input_file)

print('Leyendo archivo manzanas_pob...')
path_input_file = 'manzanas_pob.JSON'
with open(path_input_file,"r") as input_file:
    pobl = json.load(input_file)

print('Extrayendo coordenadas...')
coords = get_coords(cali)

print('Creando Voronoi...')
vor = Voronoi(coords)

print('Procesando puntos infinitos...')
regions, vertices = voronoi_finite_polygons_2d(vor)

print ('Recortando Voronoi...')
min_x = vor.min_bound[0] - 0.1
max_x = vor.max_bound[0] + 0.1
min_y = vor.min_bound[1] - 0.1
max_y = vor.max_bound[1] + 0.1
box = Polygon([[min_x, min_y], [min_x, max_y], [max_x, max_y], [max_x, min_y]])
feature_list = clip_voronoi(regions, vertices, box)

print('Escribiendo archivo de salida...')
vorJSON = open('voronoi.JSON', 'w')
feature_collection = FeatureCollection(feature_list)
print (feature_collection, file=vorJSON)
vorJSON.close()

print('Agregando Voronoi al mapa...')
folium.GeoJson(open('voronoi.JSON'), name='Diagrama de Voronoi').add_to(mapVor)
folium.LayerControl().add_to(mapVor)

print('Guardando mapa...')
mapVor.save('valencia.html')

print('Listo')
