import xml.etree.ElementTree as ET
import os
from collections import defaultdict
import numpy as np
import xml.etree.ElementTree as ET
import os 
from collections import defaultdict
import re
import scipy.io
import pandas as pd
import re

def explore_tree(element):
    """
    Objectif : Convertir un élément XML en dictionnaire récursivement.
    Args:
        element (xml.etree.ElementTree.Element): L'élément XML à convertir.
    Returns:
        dict: Dictionnaire représentant l'élément XML avec ses attributs, enfants et texte.
    """
    result = {}
    # Ajouter les attributs s'ils existent
    if element.attrib:
        result["attrib"] = element.attrib

    children = list(element)
    if children:
        for child in children:
            child_result = explore_tree(child)
            if child.tag in result:
                # Si plusieurs balises identiques, les regrouper en liste
                if isinstance(result[child.tag], list):
                    result[child.tag].append(child_result)
                else:
                    result[child.tag] = [result[child.tag], child_result]
            else:
                result[child.tag] = child_result
    else:
        # Ajouter le texte s'il existe
        text = element.text.strip() if element.text else ""
        if text:
            result["text_explanation"] = text
    return result

def dict_to_xml(dict_xml):
    """
    Objectif : Convertir un dictionnaire en élément XML récursivement.
    """
    def build_element(tag, content):
        element = ET.Element(tag)
        if isinstance(content, dict):
            for key, value in content.items():
                if key == "attrib":
                    for attr_key, attr_value in value.items():
                        element.set(attr_key, str(attr_value))
                elif key == "text_explanation":
                    element.text = str(value)
                else:
                    if isinstance(value, list):
                        for item in value:
                            child = build_element(key, item)
                            element.append(child)
                    else:
                        child = build_element(key, value)
                        element.append(child)
        else:
            # Cas où content est une valeur simple
            element.text = str(content)
        return element

    root_tag = list(dict_xml.keys())[0]
    root_content = dict_xml[root_tag]
    return build_element(root_tag, root_content)

def load_files(cell_mat_path, neighboring_cells_path):
    """
    Input :
    cell_mat_path : path to the .mat file with cell information
    neighboring_cells_path : path to the .txt file with neighboring cell information
    
    Output :
    neighboring_cells : dict {cancer_cell_id: [neighbor_ids]}
    filtered_rows : DataFrame with filtered cancer cell information
    total_cells : int, nombre total d'agents
    total_volume : float, volume total des cellules cancéreuses
    """
    neighboring_cells = None 
    filtered_rows = None
    total_cells = None
    total_volume = None
    ids_volume_dict = None
    
    try: 
        # Charger le fichier .mat avec le descriptif (e.g. les labels) de chaque agent 
        mat = scipy.io.loadmat(cell_mat_path)
        cells = mat.get("cells")  
        cells_T = cells.T #pour avoir les labels en colonne
        df = pd.DataFrame(cells_T)
        df_filtered = df.iloc[:, 0:6] #récupérer les 6 premières colonnes (id, (x,y,z), volume, type cellulaire)
        filtered_rows = df_filtered[df_filtered.iloc[:, 5] == 3.0]  # cancer cells 
        
        # IDs et volumes des cellules cancéreuses
        ids_cancer_cells = filtered_rows.iloc[:, 0].astype(int).tolist() #ids
        ids_cancer_cells_volume = filtered_rows.iloc[:, 4].astype(float).tolist() #volumes
        
        # Dictionnaire id -> volume
        ids_volume_dict = {ids_cancer_cells[i]: ids_cancer_cells_volume[i] for i in range(len(ids_cancer_cells))}
        
        # Nombre total et volume total
        total_cells = df_filtered.shape[0] #nombre total d'agents
        total_volume = df_filtered.iloc[:, 4].sum() #volume total des agents

        # Charger les voisins
        with open(neighboring_cells_path, "r") as file:
            lines = file.readlines()
            data = [line.strip().split() for line in lines]

        # Construire le dictionnaire des voisins cancéreux
        neighboring_cells = {} 
        for pos in data: 
            cell_id = int(pos[0].split(":")[0])
            if cell_id in ids_cancer_cells:
                if cell_id not in neighboring_cells:
                    neighboring_cells[cell_id] = []
                if len(pos) > 1: 
                    neighbor_ids = [int(n) for n in pos[1].split(",")]
                    for neighbor in neighbor_ids: 
                        if neighbor in ids_cancer_cells:
                            neighboring_cells[cell_id].append(neighbor)
    except Exception as e: 
        print(f"Error loading file : {cell_mat_path} or {neighboring_cells_path} : {e}")
    return neighboring_cells, total_cells, total_volume, ids_volume_dict

def clustering_cells(neighboring_cells, ids_volume_dict):
    """ 
    Input : 
    neighboring_cells : dictionnaire des voisins proches
    cell_mat : matrice avec IDs des cellules et labels

    Output : 
    clusters : vecteur avec chaque id par cluster dans le tissu conjonctif et nombre total de cellules cancéreuses
    """
    clusters = []

    for cell_id, neighbors in neighboring_cells.items():
        new_cluster = set()
        new_cluster.add(cell_id)
        for n in neighbors:
            new_cluster.add(n)

        merged = []
        for cluster in clusters:
            if not new_cluster.isdisjoint(cluster):
                new_cluster = new_cluster.union(cluster)
            else:
                merged.append(cluster)

        merged.append(new_cluster)
        clusters = merged
    
    clusters_volume = [] #Ajout de l'info du volume de chaque cancer cell présente dans chaque cluster 
    for cluster in clusters : 
        cluster_volume = {}
        for cell in cluster : 
            cluster_volume[cell] = ids_volume_dict[cell]
        clusters_volume.append(cluster_volume)
    return clusters_volume

def cancer_cell_cluster_per_time(cell_mat_path, neighboring_cells_path):
    """
    Input :
    cell_mat_path : path to the .mat file with cell information
    neighboring_cells_path : path to the .txt file with neighboring cell information
    
    Output :
    clusters : list of clusters of neighboring cancer cells (with at least 1 cancer cell)
    total_cancer_cells : total number of cancer cells
    """
    output_time_step = None
    neighboring_cells, total_cells, total_volume, ids_volume_dict = load_files(cell_mat_path, neighboring_cells_path)
    if neighboring_cells is not None :
        clusters = clustering_cells(neighboring_cells, ids_volume_dict)
        output_time_step = [clusters, total_cells, total_volume]
    return output_time_step #output at each time step 

def list_path_folder(folder_path):
    path_folder_list = os.listdir(folder_path)

    # Initialisation du dictionnaire de listes
    files_by_timestep = defaultdict(list)

    # Pattern pour récupérer l'info sur le timestep
    pattern = re.compile(r"(?:output|final|initial)0*([0-9]+)?")

    for filename in path_folder_list:
        # Cas "final" ou "initial"
        if filename.startswith("final"):
            key = "final"
        elif filename.startswith("initial"):
            key = "initial"
        else:
            # Cas "output00001108..."
            match = pattern.match(filename)
            if match:
                timestep = match.group(1)
                key = int(timestep) if timestep else 0
            else:
                continue  # ignorer les fichiers non conformes

        if filename.endswith(("_cells.mat", "_neighbor_graph.txt")):
            files_by_timestep[key].append(filename)

    # Supprimer les entrées avec un nombre de fichiers incorrect, différents de 2 
    keys_to_delete = []
    for key in files_by_timestep.keys(): 
        if len(files_by_timestep[key]) != 2: 
            keys_to_delete.append(key)

    for key in keys_to_delete: 
        del files_by_timestep[key]
    files_by_timestep = dict(files_by_timestep)
    return files_by_timestep

def get_matrix_ids(files_by_timestep, root_path):
    """
    Input :
    files_by_timestep : dict {timestep: [file1, file2]}
    root_path : path to the root folder
    Output :
    result_mat : dict {timestep: [clusters, total_cells, total_volume, timestep]}
    """
    output_path = os.path.join(root_path, "PhysiCell/output")

    result_mat = {} #initialisation du dictionnaire de résultats
    for timestep in files_by_timestep.keys():
        result_mat[timestep]= []
        try:
            file2 = os.path.join(output_path, files_by_timestep[timestep][1]) #neighboring_cells_path
            file1 = os.path.join(output_path, files_by_timestep[timestep][0]) #cell_mat_path
            result_array = cancer_cell_cluster_per_time(file1, file2)  # Input : cell_mat_path, neighboring_cells_path / doit retourner une liste ou array de taille 3 de la forme [liste avec des sets avec les ids de chaque cellule dans chaque cluster, nb d'agents, volume total]
            if result_array is not None :
                result_array.append(timestep) #ajout du timestep à la fin
                result_mat[timestep] = result_array
            else : 
                print(f"Problème à l'étape {timestep} : résultat None")
        except Exception as e:
            print(f"Erreur à l'étape {timestep} : {e}")
            print(f"file2 : {file2}")
            print(f"file1 : {file1}")
    return result_mat

def computation_area_over_time(result_mat): 
    """
    Input :
    result_mat : dictionnaire avec pour chaque pas de temps [clusters, total_cells, total_volume, timestep]
    Output :
    area_over_time : float, aire totale de la tumeur sur le temps
    """
    area_over_time = 0.0
    dt = 2 #heures
    for key in result_mat.keys():
        if result_mat[key] != []:
            for cluster in result_mat[key][0]:  #clusters
                cluster_volume = sum(cluster.values())
                area_over_time += cluster_volume*dt #(volume tumeur / volume total) * dt
    return area_over_time