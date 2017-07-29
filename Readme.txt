UNIVERSIDAD POLITÉCNICA DE VALENCIA                           
Program name: data_processing. 						     
Description: The program process all the information required to generate the
             voronoi and points of interest json files. For data visualization     
             create a Valencia's map.           
Autor: Jesús Gerardo Moreno Nieblas.				              
Date: 23/06/2017.

NOTES
Before to execute data_processing.py decompress opendata.rar on your folder. 

Also, you need to install some dependencies to run it successfully: 
- numpy
- geopandas
- folium
- scipy
- shapely

I strongly recommend you use Anaconda Distribution to avoid errors in the 
installation of libraries. I know what I'm telling you friend. 

Finally, data_processing program generates 3 files: 
* puntos_de_interes.json: Contains each Valencia's points of interest with population, traffic, tweets and time data. 
* voronoi.json: Contains Voronoi's polygons. 
- valencia.html: To visualize points of interest and voronoi json files on a web map.

The files with an asterisk symbol (*) are the input files to genetic algorithm. 