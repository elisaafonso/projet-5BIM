from scipy.io import loadmat
import numpy as np
from functions import list_path_folder, explore_tree
import xml.etree.ElementTree as ET
import difflib
import os

def compare_txt(file1, file2):
    """Compare deux fichiers texte. 
    Retourne True si identiques, sinon lève une exception avec le diff."""
    with open(file1, encoding="utf-8") as f1, open(file2, encoding="utf-8") as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    if lines1 == lines2:
        return True  # fichiers identiques

    # Si différents, générer le diff
    diff = difflib.unified_diff(
        lines1,
        lines2,
        fromfile=file1,
        tofile=file2,
        lineterm=""
    )
    diff_text = "\n".join(diff)
    raise ValueError(f"Les fichiers sont différents :\n{diff_text}")


if __name__ == "__main__":
    xml_path = "/home/vidium06/src/PhysiCell/config/PhysiCell_settings.xml"

    num_sim = 40

    tree = ET.parse(xml_path)
    root = tree.getroot()
    dict_xml = {root.tag: explore_tree(root)}

    final_time = int(float(dict_xml["PhysiCell_settings"]["overall"]["max_time"]["text_explanation"]))
    interval_time = int(float(dict_xml["PhysiCell_settings"]["save"]["full_data"]["interval"]["text_explanation"]))
	
    for j in range (num_sim): 
        print(f"Comparaison {j}")
        output_path_i_1 = f"/home/vidium06/src/ANA_SENS/Results_PhysiCell/sensitivity_analysis_descent_time_V1_fonctionnel/output_{j}"
        output_path_i_2 = f"/home/vidium06/src/ANA_SENS/Results_PhysiCell/sensitivity_analysis_descent_time_V2_fonctionnel/output_{j}"
        files_by_timestep1 = list_path_folder(output_path_i_1, final_time, interval_time)
        files_by_timestep2 = list_path_folder(output_path_i_2, final_time, interval_time)
        
        if files_by_timestep1 != files_by_timestep2: 
            raise ValueError("Folders do not contain the same amount of files...")

        for i in range (len(files_by_timestep1)): 
            files1 = files_by_timestep1[i] #lister les fichiers à ce pas de temps
            cell_file1 = next(f for f in files1 if f.endswith("_cells.mat"))
            neighbor_file1 = next(f for f in files1 if f.endswith("_cell_neighbor_graph.txt"))
                
            files2 = files_by_timestep2[i] #lister les fichiers à ce pas de temps
            cell_file2 = next(f for f in files2 if f.endswith("_cells.mat"))
            neighbor_file2 = next(f for f in files2 if f.endswith("_cell_neighbor_graph.txt"))

            file11 = os.path.join(output_path_i_1, cell_file1)
            file21 = os.path.join(output_path_i_1, neighbor_file1)

            file12 = os.path.join(output_path_i_1, cell_file1)
            file22 = os.path.join(output_path_i_1, neighbor_file1)

            #COMPARAISON .TXT

            try:
                result = compare_txt(file21, file22)
            except ValueError as e:
                raise ValueError(f"{e}")

            #COMPARAISON FICHIER .MAT

            matrix_name = "cells"   # nom de la matrice dans les fichiers .mat
            mat1 = loadmat(file11)
            mat2 = loadmat(file12)

            A1 = mat1[matrix_name]
            A2 = mat2[matrix_name]

            # 1. Même forme ?
            if A1.shape != A2.shape:
                raise ValueError("Les matrices n'ont pas la même forme.")

            # 2. Égalité stricte
            exact_equal = np.array_equal(A1, A2)
            
            if not exact_equal: 
                raise ValueError(f"Différences présentes dans les fichiers...")
        print(f"Comparaison {j} ok !")

