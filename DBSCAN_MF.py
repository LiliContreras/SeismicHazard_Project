#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: lili

Creation of seismic zones by applying HDBSCAN algorithm

"""
#---------------------------------------------------------------------------
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import hdbscan
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon
from shapely.geometry import Polygon as ShapelyPolygon
from scipy.spatial.distance import cdist
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
#To close all remaining figures
plt.close('all')
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
#Functions:

#This function selects points inside the spatial window    
def seleccionar_puntos_en_ventana(x_min, x_max, y_min, y_max, data):
    # Filter points within window
    seleccionados = data[(data['longs'] >= x_min) & (data['longs'] <= x_max) &
                         (data['lat'] >= y_min) & (data['lat'] <= y_max)]
    return seleccionados

#This function is to move the window along the spatial range
def actualizar_ventana(x, y,df):
    # Update rectangle position
    rect.set_xy((x, y))
    
    # Calls the other function and selects points inside window
    seleccionados = seleccionar_puntos_en_ventana(x, x + window_width, y, y + window_height, df)
    #This is used for visualization only:
    colores = []
    for lon, lat in zip(longitudes, latitudes):
        if ((lon >= x) and (lon <= x + window_width) and
            (lat >= y) and (lat <= y + window_height)):
            colores.append('blue')  # Los puntos seleccionados se ponen en azul
        else:
            colores.append('gray')  # Los demás se mantienen en gris
    
    # Changes colors of selected points
    sc.set_color(colores)
    plt.draw()  # Visualization purposes only
    return seleccionados
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# Load data
directorio = '/home/lili/Documents/Maestria/Tesis/Mecanismos_Focales_Catalogo/'
archivo = "2000-2023.csv" #Catalog
file = directorio + archivo
df = pd.read_csv(file)  #Reag catalog
latitudes = df['lat']
longitudes = df['longs']
df['grupo_hdbscan']=0  #Creates new column to save clusters
df = df.drop_duplicates()  #Eliminates duplicated events

# Creates figure and axes
fig, ax = plt.subplots()
sc = ax.scatter(longitudes, latitudes, c='gray', s=20) 

# Configures axis
ax.set_xlabel('Longitud')
ax.set_ylabel('Latitud')

#Window dimensions
window_width = 4.0
window_height = 5.0

# BUilds the rectangle of the movile window
rect = patches.Rectangle((min(longitudes), min(latitudes)), window_width, window_height,
                         linewidth=2, edgecolor='red', facecolor='none')
ax.add_patch(rect)

# Creates list to save results
resultado_global = []

# Starts visualization of selected points within window
plt.show(block=False)  # SHows igure

x_positions = np.arange(min(longitudes)-1, max(longitudes)+2 - window_width, 2)  #Steps by 1 and 2 units
y_positions = np.arange(min(latitudes)-1, max(latitudes)+2 - window_height, 2)

#Creates dictionary to save the created clusters
clusters_dict = {}
nclusters = []
cluster_contador = 0

# This cicle is to roam the domain with the window and apply HDBSCAN
for x in x_positions:
    for y in y_positions:
        #Extracts the points inside the window
        puntos_prueba = actualizar_ventana(x, y, df)
        if len(puntos_prueba) > 5:
            #This are the characteristics to create similiar clusters
            X = puntos_prueba[['strike1', 'dip1', 'rake1', 'prof']]
            #Apply algorithm to data
            hdbscan_cluster = hdbscan.HDBSCAN(algorithm='best', approx_min_span_tree=True, metric='euclidean', min_cluster_size=5)
            cluster_labels = hdbscan_cluster.fit_predict(X)
            #This conditional is to numerate the identified clusters
            if len(set(cluster_labels)) > 1:  
                unique_clusters = set(cluster_labels[cluster_labels >= 0])
                cluster_mapping = {}

                #New mapping of the clusters to follow global secuency
                for cluster in unique_clusters:
                    if cluster not in cluster_mapping:
                        cluster_contador += 1
                        cluster_mapping[cluster] = cluster_contador

                #Adjust labels to folow global secuency
                puntos_prueba['grupo_hdbscan'] = [
                    cluster_mapping[label] if label >= 0 else -1 for label in cluster_labels
                ]
                #Adds founded clusters to original dataframe
                df.loc[puntos_prueba.index, 'grupo_hdbscan'] = puntos_prueba['grupo_hdbscan']

                # Saves points on dictionary in each cluster
                for cluster_id in unique_clusters:
                    global_cluster_id = cluster_mapping[cluster_id]
                    cluster_points = puntos_prueba[puntos_prueba['grupo_hdbscan'] == global_cluster_id]
                    if global_cluster_id in clusters_dict:
                        clusters_dict[global_cluster_id] = pd.concat([clusters_dict[global_cluster_id], cluster_points])
                    else:
                        clusters_dict[global_cluster_id] = cluster_points

# To simplify clusters
centros_clusters = {} #Gets the center of each main cluster

for cluster_id in clusters_dict.keys():
    cluster_points = clusters_dict[cluster_id]
    centro = cluster_points[['longs', 'lat']].mean().values  #COmputes the centroid
    centros_clusters[cluster_id] = centro

centros_array = np.array(list(centros_clusters.values()))

# Distance umbral
umbral_distancia = 0.5

#DIstance matrix
distancias = cdist(centros_array, centros_array)

clusters_combinados = []

#COmbines clusters if the distance between them is smaller that the distance umbral
for i in range(len(centros_array)):
    cluster_actual = {list(centros_clusters.keys())[i]}
    for j in range(len(centros_array)):
        if i != j and distancias[i, j] < umbral_distancia:
            cluster_actual.add(list(centros_clusters.keys())[j])
    clusters_combinados.append(cluster_actual)

#Delete duplicated and combine unique clusters
clusters_combinados = [list(c) for c in set(frozenset(c) for c in clusters_combinados)]


#Create combinated clusters dictionary
nuevo_cluster_id = 1  # Empezamos con 1
clusters_combinados_dict = {}

for group in clusters_combinados:
    if len(group) > 1:  # Joins groups with more than 1 cluster
        cluster_points = pd.concat([df[df['grupo_hdbscan'] == cluster_id] for cluster_id in group])
        
        # Asign a label
        clusters_combinados_dict[nuevo_cluster_id] = cluster_points
        
        #Auxiliar variable to asign labels
        nuevo_cluster_id += 1
        
df_combinado = df.copy() #Save copy of original dataframe



nuevo_cluster_id = max(df_combinado['grupo_hdbscan']) + 1
for group in clusters_combinados:
    if len(group) > 1:  # Same procedure
        for cluster_id in group:
            df_combinado.loc[df_combinado['grupo_hdbscan'] == cluster_id, 'grupo_hdbscan'] = nuevo_cluster_id
        nuevo_cluster_id += 1
        
clusters_dict_combinados = {}


cluster_id_consecutivo = 1

# Recorrer cada cluster en el DataFrame combinado
for cluster_combined_id in df_combinado['grupo_hdbscan'].unique():
    # Filtrar las filas correspondientes a este cluster combinado
    cluster_points = df_combinado[df_combinado['grupo_hdbscan'] == cluster_combined_id]
    
    # Agregar los puntos del cluster al diccionario 
    clusters_dict_combinados[cluster_id_consecutivo] = cluster_points
    
    cluster_id_consecutivo += 1

# Visualización de los clusters en un scatter plot con convex hull
plt.figure(figsize=(10, 8))

for cluster_id in set(df['grupo_hdbscan']):
    cluster_points = df[df['grupo_hdbscan'] == cluster_id]
    
    # Asignar color diferente para cada cluster y negro para ruido (-1)
    if cluster_id == -1:
        color = 'black'
        label = 'Noise'
    else:
        color = plt.cm.tab10(cluster_id / 10)
        label = f'Cluster {cluster_id}'
    
    # Dibujar puntos del cluster
    plt.scatter(cluster_points['longs'], cluster_points['lat'], 
                color=color, label=label, alpha=0.6, edgecolors='k')

    # Dibujar convex hull si el cluster tiene suficientes puntos
    if cluster_id != -1 and len(cluster_points) > 3:  # Se necesitan al menos 3 puntos para un hull
        puntos_cluster = cluster_points[['longs', 'lat']].values
        try:
            hull = ConvexHull(puntos_cluster)
            hull_points = puntos_cluster[hull.vertices]
            polygon = Polygon(hull_points)
            plt.plot(*polygon.exterior.xy, color=color, linestyle='--', linewidth=1.5)
        except Exception as e:
            print(f"Can't create convexhull for cluster {cluster_id}: {e}")

plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Clusters HDBSCAN')
#plt.legend()
plt.grid(True)
plt.show()   
    
# Same visualization but for the simplified clusters
plt.figure(figsize=(10, 8))

for cluster_id in set(df_combinado['grupo_hdbscan']):
    cluster_points = df_combinado[df_combinado['grupo_hdbscan'] == cluster_id]
    
    # Asignar color diferente para cada cluster y negro para ruido (-1)
    if cluster_id == -1:
        color = 'black'
        label = 'Ruido'
    else:
        color = plt.cm.tab10(cluster_id / 10)
        label = f'Cluster {cluster_id}'
    
    # Dibujar puntos del cluster
    plt.scatter(cluster_points['longs'], cluster_points['lat'], 
                color=color, label=label, alpha=0.6, edgecolors='k')

    # Dibujar convex hull si el cluster tiene suficientes puntos
    if cluster_id != -1 and len(cluster_points) > 3:  # Se necesitan al menos 3 puntos para un hull
        puntos_cluster = cluster_points[['longs', 'lat']].values
        try:
            hull = ConvexHull(puntos_cluster)
            hull_points = puntos_cluster[hull.vertices]
            polygon = Polygon(hull_points)
            plt.plot(*polygon.exterior.xy, color=color, linestyle='--', linewidth=1.5)
        except Exception as e:
            print(f"Can't create convexhull for cluster {cluster_id}: {e}")

plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Clusters HDBSCAN Simplified')
#plt.legend()
plt.grid(True)
plt.show()     

#This section is no analyze clusters distribution
cluster_to_analyze = 24 #Select a cluster

#COnditional to check if selected cluster exists
if cluster_to_analyze in clusters_dict_combinados:
    cluster_data = clusters_dict_combinados[cluster_to_analyze]
    
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 7), gridspec_kw={'height_ratios': [1, 1]})
    
    # HISTOGRAMA de 'strike'
    axes[0, 0].hist(cluster_data['strike1'], bins=20, color='blue', alpha=0.7)
    axes[0, 0].set_title(f'Cluster {cluster_to_analyze} - Strike histogram')
    axes[0, 0].set_xlabel('Strike (°)')
    axes[0, 0].set_ylabel('Frequency')

    # HISTOGRAMA de 'dip'
    axes[0, 1].hist(cluster_data['dip1'], bins=20, color='green', alpha=0.7)
    axes[0, 1].set_title(f'Cluster {cluster_to_analyze} - Dip histogram')
    axes[0, 1].set_xlabel('Dip (°)')
    #axes[0, 1].set_ylabel('Frequency')

    # HISTOGRAMA de 'rake'
    axes[0, 2].hist(cluster_data['rake1'], bins=20, color='red', alpha=0.7)
    axes[0, 2].set_title(f'Cluster {cluster_to_analyze} - Rake histogram')
    axes[0, 2].set_xlabel('Rake (°)')
    #axes[0, 2].set_ylabel('Frecuencia')

    # Subfigura de la distribución de puntos en vista a profundidad
    axes[1, 0].scatter(cluster_data['longs'], cluster_data['prof'], c=cluster_data['prof'], cmap='viridis', alpha=0.6)
    axes[1, 0].invert_yaxis()  
    axes[1, 0].set_title(f'Cluster {cluster_to_analyze} - Depth distribution')
    axes[1, 0].set_xlabel('Longitude')
    axes[1, 0].set_ylabel('Depth')

    points = cluster_data[['lat', 'longs']].values
    
    # Calcular el convex hull
    hull = ConvexHull(points)
    
    axes[1, 1].scatter(cluster_data['lat'], cluster_data['longs'], c='blue', alpha=0.6)
    
    # Dibujar el polígono que encierra los puntos (Convex Hull)
    for simplex in hull.simplices:
        axes[1, 1].plot(points[simplex, 0], points[simplex, 1], 'k-', lw=2)
    
    
    # Títulos y etiquetas para la vista en planta
    axes[1, 1].set_title(f'Cluster {cluster_to_analyze} - Top view')
    axes[1, 1].set_xlabel('Latitude')
    axes[1, 1].set_ylabel('Longitude')
    axes[1, 1].grid(True)

    axes[1, 2].axis('off')
    
    # Ajustar el espacio entre subgráficas
    plt.tight_layout()
    plt.show()

else:
    print(f"Cluster {cluster_to_analyze} doedn't exists.")

#This code is extra for exporting clusters and make a map on GMT
'''
#Save geometries for GMT
directorio = '/home/lili/Documents/Maestria/Tesis/Mecanismos_Focales_Catalogo/'
archivo = "NewClustersP.txt"
fileN = directorio + archivo
umbral_area = 10.0  # Ajusta este valor según sea necesario

with open(fileN, "w") as file:
    for cluster_id, cluster_data in clusters_dict_combinados.items():
        if len(cluster_data) >= 3:  # Convex hulls necesitan al menos 3 puntos
            hull = ConvexHull(cluster_data[['longs', 'lat']])
            hull_points = cluster_data.iloc[hull.vertices][['lat', 'longs']]
            
            # Crear un polígono de Shapely para calcular el área
            polygon = Polygon(hull_points[['longs', 'lat']].values)
            
            # Verificar si el área del polígono es menor que el umbral
            if polygon.area <= umbral_area:
                # Escribir un marcador para indicar el inicio de un nuevo polígono
                file.write(f"> Cluster {cluster_id}\n")
                
                # Escribir las coordenadas del convex hull
                for _, row in hull_points.iterrows():
                    file.write(f"{row['longs']} {row['lat']}\n")
                
                # Repetir el primer punto para cerrar el polígono visualmente en GMT
                file.write(f"{hull_points.iloc[0]['longs']} {hull_points.iloc[0]['lat']}\n")
'''