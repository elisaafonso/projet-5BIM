# A PENSER A MODIFIER : La position du tissu conjonctif dans les conditions initiales et les ids des cellules s'ils ont été modifiés. 

from SALib.analyze.sobol import analyze
from SALib.sample.sobol import sample
import numpy as np
import xml.etree.ElementTree as ET
from functions import explore_tree, dict_to_xml, list_path_folder, get_matrix_ids, computation_area_over_time, list_mat_files, compute_time_ratio, get_result_mat_persistance
import os 
import subprocess
from collections import defaultdict
import re
import scipy.io
import pandas as pd
import re
import shutil

def define_set_param(): 
    print("Defining parameter space and sampling...")

    if analyse_sensibilite == "descent_time" : 
        problem = {
            'num_vars': 5,
            'names': ['attachment_rate', 'cancer_motility_speed', 'transformation_rate_mes', 'cell_adhesion_affinity_LM_EP', 'cell_adhesion_affinity_LB_conj'],
            'bounds': [[720, 2880], [0.0000001, 0.000001]] #TO CHANGE
        }
    elif analyse_sensibilite == "tumor_persistance" :
        problem = {
            'num_vars': 2,
            'names': ['motility_speed_t_cell', 'division_duration_cancer'],
            'bounds': [[720, 2880], [0.0000001, 0.000001]] #TO CHANGE
        }

    param_values = sample(problem, 2) #Génère N*(2+D) jeux de paramètres avec D le nombre de paramètres et N un multiple de 2 fourni en argument
    return problem, param_values

def define_settings(parameters_to_change): 
    print("Defining XML settings for the simulation...")
    tree = ET.parse(os.path.join(root_path, "PhysiCell/config/PhysiCell_settings.xml"))
    root = tree.getroot()
    dict_xml = {root.tag: explore_tree(root)}

    #changement des paramètres 

    ################## DESCENT TIME ####################
    if analyse_sensibilite == "descent_time" :
        # Attachment rate 
        attachment_rate = parameters_to_change[0]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][5]["phenotype"]["mechanics"]["attachment_rate"]["text_explanation"] = attachment_rate

        #Motility speed cancer cells
        cancer_motility_speed = parameters_to_change[1]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][3]["phenotype"]["motility"]["speed"]["text_explanation"] = cancer_motility_speed
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][8]["phenotype"]["motility"]["speed"]["text_explanation"] = cancer_motility_speed

        # Mes transition 
        transformation_rate_mes = parameters_to_change[2]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][3]["phenotype"]["cell_transformations"]["transformation_rates"]["transformation_rate"][8]["text_explanation"] = transformation_rate_mes
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][8]["phenotype"]["cell_transformations"]["transformation_rates"]["transformation_rate"][3]["text_explanation"] = transformation_rate_mes

        #Adhésion lame basale/épi basal
        cell_adhesion_affinity_LM_EP = parameters_to_change[3]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][5]["phenotype"]["mechanics"]["cell_adhesion_affinities"]["cell_adhesion_affinity"][0]["text_explanation"] = cell_adhesion_affinity_LM_EP
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][0]["phenotype"]["mechanics"]["cell_adhesion_affinities"]["cell_adhesion_affinity"][5]["text_explanation"] = cell_adhesion_affinity_LM_EP

        #Adhésion MB/Conj 
        cell_adhesion_affinity_LB_conj = parameters_to_change[4]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][5]["phenotype"]["mechanics"]["cell_adhesion_affinities"]["cell_adhesion_affinity"][6]["text_explanation"] = cell_adhesion_affinity_LB_conj
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][6]["phenotype"]["mechanics"]["cell_adhesion_affinities"]["cell_adhesion_affinity"][5]["text_explanation"] = cell_adhesion_affinity_LB_conj

    ################## TUMOR PERSISTANCE ####################
    elif analyse_sensibilite == "tumor_persistance" :
        # Vitesse motilité T-Cell 
        motility_speed_t_cell = parameters_to_change[0]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][4]["phenotype"]["motility"]["speed"]["text_explanation"] = motility_speed_t_cell

        # Division cellules cancéreuses (mes et non mes)
        division_duration_cancer = parameters_to_change[0]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][3]["phenotype"]["cycle"]["phase_durations"]["duration"]["text_explanation"] = division_duration_cancer
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][8]["phenotype"]["cycle"]["phase_durations"]["duration"]["text_explanation"] = division_duration_cancer

    #dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][0]["phenotype"]["cell_transformations"]["transformation_rates"]["transformation_rate"][3]["text_explanation"] = transformation_rate_epi_basal_cancer
    #dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][3]["phenotype"]["cycle"]["phase_durations"]["duration"]["text_explanation"] = life_cycle_cancer_cell
    #dict_xml["PhysiCell_settings"]["overall"]["max_time"]["text_explanation"] = 240 #modification du temps de simulation
    return dict_xml

def get_physicell_output(param_values): 
    output_storage = np.zeros(param_values.shape[0]) 
    for i in range (len(param_values)): 
        print(f"Running simulation {i+1} with parameters: life_cycle_cancer_cell={param_values[i][0]}, transformation_rate_epi_basal_cancer={param_values[i][1]}")
        #life_cycle_cancer_cell = param_values[i][0]
        #transformation_rate_epi_basal_cancer = param_values[i][1]
        parameters_to_change = param_values[i]
        dict_xml = define_settings(parameters_to_change)

        xml_element = dict_to_xml(dict_xml)
        tree = ET.ElementTree(xml_element)
        tree.write(os.path.join(root_path, "PhysiCell/config/PhysiCell_settings.xml"), encoding="utf-8", xml_declaration=True) #path à changer

        process0 = subprocess.run(
        ["rm", "-rf", "*"],
        capture_output=True,
        text=True,
        cwd=os.path.join(root_path, "PhysiCell/output") #to clean the output folder before each simulation
        )

        print("Configuration file updated. Starting simulation...")
        # Run .exe file 
        process1 = subprocess.run(
        [os.path.join(root_path, "PhysiCell/project.exe")], # !! To change if linux 
        capture_output=True,
        text=True,
        cwd=os.path.join(root_path, "PhysiCell")
        )

        # Print the output and error (if any)
        print("Output:")
        print(process1.stdout) #vérifier le type
        print("Error:")
        print(process1.stderr)

        #get output from output folder 
        output_path = os.path.join(root_path, "PhysiCell/output") 

        if analyse_sensibilite == "tumor_persistance" :  
            pos_conj = 0 # A MODIFIER 
            files_by_timestep = list_mat_files(output_path)
            result_mat = get_result_mat_persistance(files_by_timestep, root_path, pos_conj)
            dt = int(dict_xml["PhysiCell_settings"]["save"]["full_data"]["interval"]["text_explanation"])/60.0  # conversion en heures
            metric = computation_area_over_time(result_mat, dt)
        
        elif analyse_sensibilite == "descent_time" : 
            files_by_timestep = list_mat_files(output_path)
            ratio = compute_time_ratio(files_by_timestep, root_path, position_conj = -200)  
            metric = ratio

        output_storage[i] = metric #computation metric

        #stocker la vidéo 
        process2 = subprocess.run(
        ["make", "jpeg"],
        capture_output=True,
        text=True,
        cwd=os.path.join(root_path, "PhysiCell")
        )
        process3 = subprocess.run(
        ["make", "gif"],
        capture_output=True,
        text=True,
        cwd=os.path.join(root_path, "PhysiCell")
        )
        process4 = subprocess.run(
        ["make", "movie"],
        capture_output=True,
        text=True,
        cwd=os.path.join(root_path, "PhysiCell")
        )

        print("Output:")
        print(f"{process2.stdout}, {process3.stdout}, {process4.stdout}")
        print("Error:")
        print(f"{process2.stderr}, {process3.stderr}, {process4.stderr}")

        src = os.path.join(root_path, "PhysiCell/output/out.mp4")
        dst_folder = os.path.join(root_path, "output_video")
        dst = os.path.join(dst_folder, f"movie_simulation_{i+1}.mp4")
        os.makedirs(dst_folder, exist_ok=True)
        # Move video file and rename it
        shutil.move(src, dst)

    return output_storage

def analyze_sobol():
    problem, param_values = define_set_param()
    Y = get_physicell_output(param_values)
    Si = analyze(problem, Y, print_to_console=True)
    return Si

if __name__ == "__main__":
    analyse_sensibilite = "tumor_persistance"
    root_path = "C://Users/elisa/" # !! A changer
    #get le xml et le copier dans config/ 
    #ouvrir interface graphique Tkinter
    #créer un fichier log
    
    Si = analyze_sobol()
    print(f"Sobol Indices: {Si}")