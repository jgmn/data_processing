import json
import folium
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import box
from scipy.spatial import Voronoi, voronoi_plot_2d
from geojson import FeatureCollection, Feature, Polygon

def add_markers(data):
    coords = []
    for feature in data['features']:
        description = "califi: "+ feature['properties']['califi'] + " " + "tipoca: "+ feature['properties']['tipoca']
        coordinate = feature['geometry']['coordinates']
        x, y = coordinate[1], coordinate[0]
        #folium.Marker([x,y], popup = description).add_to(mapVor) 
        coords.append([x,y])
    return coords

def clip_voronoi(regions, vertices, bounding_box):
    feature_list = []
    for region in regions:
        vertice_list = []
        for indice in region:
            vertice = vertices[indice]
            vertice = (vertice[1], vertice[0])
            vertice_list.append(vertice)
        polygon = Polygon([vertice_list])
        
        feature = Feature(geometry=polygon, properties={})
        feature_list.append(feature)

    feature_collection = FeatureCollection(feature_list)
    print (feature_collection, file = vorJSON)
    vorJSON.close()

def voronoi_finite_polygons_2d(vor, radius=None):

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all([v >= 0 for v in vertices]):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices).tolist()
    
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
vor = Voronoi(coords)
regions, vertices = voronoi_finite_polygons_2d(vor)
bounding_box = box()
clip_voronoi(regions, vertices, bounding_box)

"""
# colorize
for region in regions:
    polygon = vertices[region]
    plt.fill(*zip(*polygon), alpha=0.4)

plt.plot(points[:,0], points[:,1], 'ko')
plt.axis('equal')
plt.xlim(vor.min_bound[0] - 0.1, vor.max_bound[0] + 0.1)
plt.ylim(vor.min_bound[1] - 0.1, vor.max_bound[1] + 0.1)

plt.show()
"""
"""print('Adding Voronoi diagram to the map...')
folium.GeoJson(open('voronoi.JSON'), name='Diagrama de Voronoi').add_to(mapVor)
folium.LayerControl().add_to(mapVor)

print('Saving map...')
mapVor.save('valencia.html')

print('Ready')"""
