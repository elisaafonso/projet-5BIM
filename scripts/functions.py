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

#*******************************************************************************

# ############# INTERACTION AVEC LE FICHIER XML (LECTURE, GENERATION) ##########

#*******************************************************************************

def explore_tree(element):
    """
    Objectif : Convertir un élément XML en dictionnaire récursivement.
    Args:
        element (xml.etree.ElementTree.Element): L'élément XML à convertir.
    Returns:
        dict: Dictionnaire représentant l'élément XML avec ses attributs, enfants et texte.
    Note : un document xml est composé de balises et chaque balise peut avoir des attributs, du texte (text_explanation) et des balises enfants. 
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

def get_cell_type(xml_path): 
    """
    Input : xml_path, chemin vers le fichier PhysiCell_settings.xml

    Output : 
    dict_correspondance : dict {cell_type_name: cell_type_ID}, correspondance entre le nom des types cellulaires et leur ID
    dict_corres_microenv : dict {microenv_variable_name: microenv_variable_ID}, correspondance entre le nom des variables du microenvironnement et leur ID
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()
    dict_xml = {root.tag: explore_tree(root)}
    dict_correspondance = {}
    for i in range (len(dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"])): 
        dict_correspondance[dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][i]["attrib"]["name"]] = int(dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][i]["attrib"]["ID"])

    dict_corres_microenv = {} 
    for i in range (len(dict_xml["PhysiCell_settings"]["microenvironment_setup"]["variable"])): 
        dict_corres_microenv[dict_xml["PhysiCell_settings"]["microenvironment_setup"]["variable"][i]["attrib"]["name"]] = int(dict_xml["PhysiCell_settings"]["microenvironment_setup"]["variable"][i]["attrib"]["ID"])
    dict_corres_microenv
    return dict_correspondance, dict_corres_microenv

#*******************************************************************************

############################## LISTER LES FICHIERS (.mat, .txt) ################

#*******************************************************************************

def list_mat_files(folder_path, final_time, interval_time):
    """
    Input :
    folder_path : Chemin vers le dossier contenant les fichiers de sortie ./output
    final_time : temps final de simulation
    interval_time : intervalle de temps entre deux sauvegardes

    Output :
    files_by_timestep : dict {timestep: [file1]} #file1 = fichier cell_mat file pour chaque pas de temps
    """
    path_folder_list = os.listdir(folder_path)

    # Initialisation du dictionnaire de listes
    files_by_timestep = defaultdict(list)

    # Pattern pour récupérer l'info sur le timestep
    pattern = re.compile(r"(?:output)0*([0-9]+)?")

    for filename in path_folder_list:
        match = pattern.match(filename)
        if match:
            timestep = match.group(1)
            key = int(timestep) if timestep else 0
        else:
            continue  # ignorer les fichiers non conformes

        if filename.endswith("_cells.mat"):
            files_by_timestep[key].append(filename)
    print(files_by_timestep)
    
    # Supprimer les entrées avec un nombre de fichiers incorrect, différents de 2 
    keys_to_delete = []
    for key in files_by_timestep.keys(): 
        if len(files_by_timestep[key]) != 1: 
            keys_to_delete.append(key)

    for key in keys_to_delete: 
        del files_by_timestep[key]
    
    nb_timestep = final_time // interval_time
    nb_files = nb_timestep + 1
    files_by_timestep = dict(files_by_timestep)

    if len(files_by_timestep) != nb_files:
        raise ValueError(f"Expected {nb_files} files, but found {len(files_by_timestep)} files.")
    
    return files_by_timestep

def list_path_folder(folder_path, final_time, interval_time):
    """
    Input :
    folder_path : Chemin vers le dossier contenant les fichiers de sortie ./output
    final_time : temps final de simulation
    interval_time : intervalle de temps entre deux sauvegardes
    
    Output :
    files_by_timestep : dict {timestep: [file1]} #file1 = fichier cell_mat file pour chaque pas de temps
    """
    path_folder_list = os.listdir(folder_path)

    # Initialisation du dictionnaire de listes
    files_by_timestep = defaultdict(list)

    # Pattern pour récupérer l'info sur le timestep
    pattern = re.compile(r"(?:output)0*([0-9]+)?")

    for filename in path_folder_list:
        match = pattern.match(filename)
        if match:
            timestep = match.group(1)
            key = int(timestep) if timestep else 0
        else:
            continue  # ignorer les fichiers non conformes

        if filename.endswith("_cells.mat") or filename.endswith("_cell_neighbor_graph.txt"):
            files_by_timestep[key].append(filename)


    # Supprimer les entrées avec un nombre de fichiers incorrect, différents de 2 
    keys_to_delete = []
    for key in files_by_timestep.keys(): 
        if len(files_by_timestep[key]) != 2: 
            keys_to_delete.append(key)

    for key in keys_to_delete: 
        del files_by_timestep[key]
    files_by_timestep = dict(files_by_timestep)

    nb_timestep = final_time // interval_time
    nb_files = nb_timestep + 1
    files_by_timestep = dict(files_by_timestep)

    if len(files_by_timestep) != nb_files:
        raise ValueError(f"Expected {nb_files} files, but found {len(files_by_timestep)} files.")
    
    return files_by_timestep


#*******************************************************************************

############## CALCUL DE LA METRIQUE POUR LA DESCENTE DE LA TUMEUR ##############

#*******************************************************************************

def get_cancer_cell_position(cell_mat_path, neighboring_cells_path, id_cancer_cell, id_cancer_cell_mes, id_connective_tissue, id_fibroblast):
    """
    Input :
    cell_mat_path : Chemin vers le fichier .mat avec les informations des cellules
    neighboring_cells_path : Chemin vers le fichier des cellules voisines
    id_cancer_cell : ID des cellules cancéreuses (non-mésenchymateuses)
    id_cancer_cell_mes : ID des cellules cancéreuses mésenchymateuses

    Output :
    position_y : position y des cellules cancéreuses (None si pas de cellules cancéreuses)
    ratio : temps de descente dans le tissu conjonctif / temps total de simulation
    """
    position_y = None
    is_in_connective_tissue = True 
    try: 
        mat = scipy.io.loadmat(cell_mat_path)
        cells = mat.get("cells")  
        cells_T = cells.T #pour avoir les labels en colonne
        df = pd.DataFrame(cells_T)
        df_filtered = df.iloc[:, 0:6] #récupérer les 6 premières colonnes (id, (x,y,z), volume, type cellulaire)
        filtered_rows = df_filtered[(df_filtered.iloc[:, 5] == id_cancer_cell) | (df_filtered.iloc[:, 5] == id_cancer_cell_mes)]  # cancer cells (mesenchymal or not)

        n = len(filtered_rows)
        if n > 1 : 
            raise ValueError("More than one cancer cell found.")
        elif n == 0 :
            return None, False  # no cancer cell found
        else : 
            position_y = filtered_rows.iloc[0, 2] #Position y de la cellule cancéreuse

            # Charger les voisins
            with open(neighboring_cells_path, "r") as file:
                lines = file.readlines()
                data = [line.strip().split() for line in lines]

            neighbor_ids = []
            for pos in data:
                if len(pos) < 2:
                    continue  # ligne vide ou mal formée
                try:
                    cell_id = int(pos[0].split(":")[0])
                except:
                    continue
                if cell_id == int(filtered_rows.iloc[0,0]):
                    # gérer le cas "953:" sans voisins
                    if ":" in pos[0] and pos[0].endswith(":") and len(pos) == 1:
                        neighbor_ids = []
                    else:
                        neighbor_ids = [int(n) for n in pos[1].split(",") if n.strip().isdigit()]
                    break
            if len(neighbor_ids) == 0 : 
                is_in_connective_tissue = False
            else : 
                for neighbor_id in neighbor_ids : 
                    cell_row = df_filtered[df_filtered.iloc[:,0] == neighbor_id]
                    cell_type = cell_row.iloc[0, 5]
                    if cell_type != id_connective_tissue and cell_type != id_fibroblast: 
                        is_in_connective_tissue = False
                        

    except Exception as e: 
        raise ValueError(f"Error loading file : {cell_mat_path} : {e}")
    return position_y, is_in_connective_tissue

def compute_time_ratio(files_by_timestep, output_path_i, id_cancer_cell, id_cancer_cell_mes, id_connective_tissue, id_fibroblast):
    """
    Input :
    files_by_timestep : dict {timestep: [file1, file2]}, file1 : fichier .mat, file2 : neighboring_cells.txt
    output_path_i : Chemin vers le dossier output_{i}
    id_cancer_cell : ID des cellules cancéreuses (non-mésenchymateuses)
    id_cancer_cell_mes : ID des cellules cancéreuses mésenchymateuses

    Output :
    ratio : temps de descente dans le tissu conjonctif / temps total de simulation
    """
    output_path = output_path_i
    timesteps = files_by_timestep.keys()
    simulation_time = max(timesteps) 
    ratio = 1 #if is_in_connective_tissue is always False, then ratio = 1
    result_mat = {} #initialisation du dictionnaire de résultats
    for timestep in files_by_timestep.keys():
        result_mat[timestep]= []
        try:
            files = files_by_timestep[timestep] #lister les fichiers à ce pas de temps
            cell_file = next(f for f in files if f.endswith("_cells.mat"))
            neighbor_file = next(f for f in files if f.endswith("_cell_neighbor_graph.txt"))

            file1 = os.path.join(output_path, cell_file)
            file2 = os.path.join(output_path, neighbor_file)

            position_y, is_in_connective_tissue = get_cancer_cell_position(file1, file2, id_cancer_cell, id_cancer_cell_mes, id_connective_tissue, id_fibroblast)  # Input : cell_mat_path 
            if position_y is not None : #cela veut dire qu'il y a une cellule cancéreuse
                if is_in_connective_tissue : 
                    time_in_conj_tissue = timestep
                    ratio = time_in_conj_tissue / simulation_time
                    break
            else : 
                ratio = 1
                break #pas de cellule cancéreuse trouvée, tissu conjonctif pas atteint
        except Exception as e:
            raise ValueError(f"Erreur à l'étape {timestep} : {e} \n file1 : {file1}")
    return ratio  

#*******************************************************************************

############## CALCUL DE LA METRIQUE POUR LA PERSISTANCE DE LA TUMEUR ###########

#*******************************************************************************

def get_cancer_volume_per_timestep(cell_mat_path, id_cancer, id_cancer_mes):
    """
    Input :
    cell_mat_path : Chemin vers le fichier .mat avec les informations sur les agents

    Output :
    ids_volume_dict : dict {cancer_cell_id: volume}
    """
    filtered_rows = None
    ids_volume_dict = None
    
    try: 
        # Charger le fichier .mat avec le descriptif (e.g. les labels) de chaque agent 
        mat = scipy.io.loadmat(cell_mat_path)
        cells = mat.get("cells")  
        cells_T = cells.T #pour avoir les labels en colonne
        df = pd.DataFrame(cells_T)
        df_filtered = df.iloc[:, 0:6] #récupérer les 6 premières colonnes (id, (x,y,z), volume, type cellulaire)
        filtered_rows = df_filtered[(df_filtered.iloc[:, 5] == id_cancer) | (df_filtered.iloc[:, 5] == id_cancer_mes)] # cancer cells (mesenchymal or not)

        # IDs et volumes des cellules cancéreuses
        ids_cancer_cells = filtered_rows.iloc[:, 0].astype(int).tolist() #ids  
        ids_cancer_cells_volume = filtered_rows.iloc[:, 4].astype(float).tolist() #volumes
        ids_volume_dict = {ids_cancer_cells[i]: ids_cancer_cells_volume[i] for i in range(len(ids_cancer_cells))} # Dictionnaire id -> volume

    except Exception as e: 
        raise ValueError(f"Error loading file : {cell_mat_path} : {e}")
    return ids_volume_dict

def get_result_mat_persistance(files_by_timestep, output_path_i, id_cancer, id_cancer_mes):
    """
    Input :
    files_by_timestep : dict {timestep: [file1, file2]}
    root_path : Chemin vers le dossier racine

    Output :
    result_mat : dict {timestep: {cancer_cell_ids: volumes}}
    """
    output_path = output_path_i

    result_mat = {} #initialisation du dictionnaire de résultats
    for timestep in files_by_timestep.keys():
        result_mat[timestep]= []
        try:
            file1 = os.path.join(output_path, files_by_timestep[timestep][0]) #cell_mat_path
            result_array = get_cancer_volume_per_timestep(file1, id_cancer, id_cancer_mes)  # Input : cell_mat_path / retourne : dict {cancer_cell_ids: volumes}
            if result_array is not None :
                result_mat[timestep] = result_array
            else : 
                raise ValueError(f"Problème à l'étape {timestep} : résultat None")
        except Exception as e:
            raise ValueError(f"Erreur à l'étape {timestep} : {e} \n file1 : {file1}")
    return result_mat

def computation_area_over_time(result_mat, dt): 
    """
    Input :
    result_mat : dict {timestep: {cancer_cell_ids: volumes}
    dt : intervalle de temps entre deux sauvegardes

    Output :
    area_over_time : float, aire totale de la tumeur sur le temps de simulation
    """
    area_over_time = 0.0
    for key in result_mat.keys():
        if result_mat[key] != {}:
            volume = sum(result_mat[key].values())
            area_over_time += volume*dt #(volume tumeur / volume total) * dt
    return area_over_time

if __name__ == "__main__":
    
    pass
