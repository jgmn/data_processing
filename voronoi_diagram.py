import json
import folium
from scipy.spatial import Voronoi, voronoi_plot_2d
from geojson import FeatureCollection, Feature, Polygon

def add_markers(data):
    coords = []
    for feature in data['features']:
        description = "califi: "+ feature['properties']['califi'] + " " + "tipoca: "+ feature['properties']['tipoca']
        coordinate = feature['geometry']['coordinates']
        x, y = coordinate[1], coordinate[0]
        folium.Marker([x,y], popup = description).add_to(mapVor)
        coords.append([x,y])
    return coords

def create_voronoi(coords):
    vorJSON = open('voronoi.JSON', 'w')
    feature_list = []
    vor = Voronoi(coords)
    for region in vor.regions
        vertice_list = []
        for indice in region:
            vertice = vor.vertices[indice]
            vertice = (vertice[1], vertice[0])
            vertice_list.append(vertice)
        polygon = Polygon([vertice_list])
        feature = Feature(geometry=polygon, properties={})
        feature_list.append(feature)

    feature_collection = FeatureCollection(feature_list)
    print (feature_collection, file = vorJSON)
    vorJSON.close()
    
# Main
print('Creating empty map...')
valencia = [39.4561165311493, -0.3545661635]
mapVor = folium.Map(location = valencia, zoom_start = 13)

print('Reading JSON file...')
path_input_file = 'calificaciones.JSON'
with open(path_input_file,"r") as input_file:
    data = json.load(input_file)

print('Adding markers to the map...')
coords = add_markers(data)

print('Creating Voronoi diagram...')
create_voronoi(coords)

print('Adding Voronoi diagram to the map...')
folium.GeoJson(open('voronoi.JSON'), name='Diagrama de Voronoi').add_to(mapVor)
folium.LayerControl().add_to(mapVor)

print('Saving map...')
mapVor.save('valencia.html')

print('Ready')
