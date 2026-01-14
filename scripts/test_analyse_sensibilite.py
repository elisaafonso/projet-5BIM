from SALib.analyze.sobol import analyze
from SALib.sample.sobol import sample
import numpy as np
import xml.etree.ElementTree as ET
from functions import explore_tree, dict_to_xml, list_path_folder, get_matrix_ids, computation_area_over_time, list_mat_files, compute_time_ratio, get_result_mat_persistance, get_cell_type
import os 
import subprocess
from collections import defaultdict
import re
import scipy.io
import pandas as pd
import re
import shutil
import sys
import json

def define_set_param(analyse_sensibilite, param_bounds, N): 
    print("Defining parameter space and sampling...")

    if analyse_sensibilite == "descent_time" : 
        problem = {
            'num_vars': 4,
            'names': ['attachment_rate', 'cancer_motility_speed', 'transformation_rate_mes', 'sensitivity_mmp_factor'],
            'bounds': [param_bounds[0], param_bounds[1], param_bounds[2], param_bounds[3]] 
        }
    elif analyse_sensibilite == "tumor_persistance" :
        problem = {
            'num_vars': 2,
            'names': ['motility_speed_t_cell', 'division_duration_cancer'],
            'bounds': [param_bounds[0], param_bounds[1]] 
        }

    param_values = sample(problem, N, seed = 1) #Génère N*(2+D) jeux de paramètres avec D le nombre de paramètres et N un multiple de 2 fourni en argument
    return problem, param_values

def define_settings(parameters_to_change, nb_threads, seed, analyse_sensibilite, xml_path, output_path_i, dict_corresp_name_type, dict_corres_microenv): 
    print("Defining XML settings for the simulation...")
    tree = ET.parse(xml_path)
    root = tree.getroot()
    dict_xml = {root.tag: explore_tree(root)}

    #changement des paramètres 

    ################## DESCENT TIME ####################
    if analyse_sensibilite == "descent_time" :
        # Attachment rate : 0 et 10
        attachment_rate = parameters_to_change[0]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["membrane"]]["phenotype"]["mechanics"]["attachment_rate"]["text_explanation"] = attachment_rate

        #Motility speed cancer cells : 0.01 et 0.5
        cancer_motility_speed = parameters_to_change[1]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer"]]["phenotype"]["motility"]["speed"]["text_explanation"] = cancer_motility_speed
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["motility"]["speed"]["text_explanation"] = cancer_motility_speed

        # Mes transition : 0.1 et 1  
        transformation_rate_mes = parameters_to_change[2]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer"]]["phenotype"]["cell_transformations"]["transformation_rates"]["transformation_rate"][8]["text_explanation"] = transformation_rate_mes
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["cell_transformations"]["transformation_rates"]["transformation_rate"][3]["text_explanation"] = transformation_rate_mes

        # Secretion cancer mes -- 0 et 10
        sensitivity_mmp_factor = parameters_to_change[3]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["secretion"]["substrate"][dict_corres_microenv["mmp_factor"]]["secretion_rate"]["text_explanation"] = sensitivity_mmp_factor

    ################## TUMOR PERSISTANCE ########################################
    #mettre un bord solide tissu conjonctif pour éviter // enlever CAF 
    elif analyse_sensibilite == "tumor_persistance" :
        # Vitesse motilité T-Cell 
        motility_speed_t_cell = parameters_to_change[0]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][4]["phenotype"]["motility"]["speed"]["text_explanation"] = motility_speed_t_cell

        # Division cellules cancéreuses (mes et non mes)
        division_duration_cancer = parameters_to_change[0]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][3]["phenotype"]["cycle"]["phase_durations"]["duration"]["text_explanation"] = division_duration_cancer
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][8]["phenotype"]["cycle"]["phase_durations"]["duration"]["text_explanation"] = division_duration_cancer

        #mort des cellules cancéreuses 

    ################## CHANGER LE NOMBRE DE COEUR ET LA SEED  ####################

    dict_xml["PhysiCell_settings"]["parallel"]["omp_num_threads"]["text_explanation"] = nb_threads
    dict_xml["PhysiCell_settings"]["options"]["random_seed"]["text_explanation"] = seed
    
    #define output_path
    dict_xml["PhysiCell_settings"]["save"]["folder"]["text_explanation"] = output_path_i #output
    return dict_xml

def get_physicell_output(param_values, output_video_folder, nb_threads, seed, root_path, analyse_sensibilite, xml_path, dst_folder, dict_corresp_name_type, dict_corres_microenv): 
    output_storage = np.zeros(param_values.shape[0]) 
    for i in range (len(param_values)): 
        print(f"Running simulation {i+1}...")
        #life_cycle_cancer_cell = param_values[i][0]
        #transformation_rate_epi_basal_cancer = param_values[i][1]
        output_path_i = os.path.normpath(os.path.join(dst_folder, f"output_{i}"))
        os.makedirs(output_path_i, exist_ok=True)

        parameters_to_change = param_values[i]
        dict_xml = define_settings(parameters_to_change, nb_threads, seed, analyse_sensibilite, xml_path, output_path_i, dict_corresp_name_type, dict_corres_microenv)
        xml_element = dict_to_xml(dict_xml)
        tree = ET.ElementTree(xml_element)
        tree.write(os.path.join(root_path, f"PhysiCell/config/PhysiCell_settings_{i}.xml"), encoding="utf-8", xml_declaration=True) #path à changer

        #rajouter make load PROJ=5_bim
        process0 = subprocess.run(
        ["rm", "-rf", "*"],
        capture_output=True,
        text=True,
        cwd=output_path_i #to clean the output folder before each simulation
        )
        #*******************

        print("Configuration file updated. Starting simulation...")
        # Run .exe file 
        exe_path = os.path.join(root_path, "PhysiCell", "project.exe")
        xml_path_i = os.path.join(root_path, "PhysiCell", "config", f"PhysiCell_settings_{i}.xml")

        process1 = subprocess.run(
            [exe_path, xml_path_i],
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
        output_path = output_path_i
        #os.path.join(root_path, "PhysiCell/output") 

        final_time = int(float(dict_xml["PhysiCell_settings"]["overall"]["max_time"]["text_explanation"]))
        interval_time = int(float(dict_xml["PhysiCell_settings"]["save"]["full_data"]["interval"]["text_explanation"]))

        if analyse_sensibilite == "tumor_persistance" :  
            pos_conj = 0 # A MODIFIER 
            files_by_timestep = list_mat_files(output_path, final_time, interval_time)
            result_mat = get_result_mat_persistance(files_by_timestep, root_path, pos_conj)
            dt = int(dict_xml["PhysiCell_settings"]["save"]["full_data"]["interval"]["text_explanation"])/60.0  # conversion en heures
            metric = computation_area_over_time(result_mat, dt)
        
        elif analyse_sensibilite == "descent_time" : 
            id_cancer_cell = dict_corresp_name_type["cancer"]
            id_cancer_cell_mes = dict_corresp_name_type["cancer_mes"]
            files_by_timestep = list_mat_files(output_path, final_time, interval_time)
            position_conj = -215
            ratio = compute_time_ratio(files_by_timestep, output_path_i, position_conj, id_cancer_cell, id_cancer_cell_mes)  
            metric = ratio

        output_storage[i] = metric #computation metric
        
        #stocker la vidéo 
        """process2 = subprocess.run(
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
        dst = os.path.join(output_video_folder, f"movie_simulation_{i+1}.mp4")
        # Move video file and rename it
        shutil.move(src, dst)""" #modifier le fichier output donc plus possible de faire make 

    return output_storage

def analyze_sobol(output_video_folder, nb_threads, seed, root_path, analyse_sensibilite, param_bounds, xml_path, dst_folder, N, dict_corresp_name_type, dict_corres_microenv):
    problem, param_values = define_set_param(analyse_sensibilite, param_bounds, N)
    
    #Save parameters values to a CSV file
    df = pd.DataFrame(param_values, columns=[p['name'] for p in problem['names']]) 
    df.to_csv(os.path.join(dst_folder, "param_values.csv"), index=False)

    Y = get_physicell_output(param_values, output_video_folder, nb_threads, seed, root_path, analyse_sensibilite, xml_path, dst_folder, dict_corresp_name_type, dict_corres_microenv)
    np.savetxt(os.path.join(dst_folder, "output_values.csv"), Y, delimiter=",")

    Si = analyze(problem, Y, print_to_console=True, seed = 1)
    return Si

def run_main_analysis(params):
    print("Lancement de l'analyse de sensibilité...")
    analyse_sensibilite = params["analyse_sensibilite"]
    param_bounds = params["param_bounds"]
    nb_threads = params["nb_threads"]
    seed = params["seed"]
    root_path = params["physicell_path"]
    result_folder_path = params["results_path"]
    xml_path = params["xml_path"]
    N = params["nb_sample_to_generate"] #a power of 2

    print("Analyse :", analyse_sensibilite)
    print("Param bounds :", param_bounds)
    print("Threads :", nb_threads)
    print("Seed :", seed)
    print("PhysiCell path :", root_path)
    print("Results path :", result_folder_path)
    print("XML path :", xml_path)

    dict_corresp_name_type, dict_corres_microenv = get_cell_type(xml_path)

    print("***************************")
    print(dict_corres_microenv)
    print("*************************")

    dst_folder = os.path.join(result_folder_path, "Results_PhysiCell", f"sensitivity_analysis_{analyse_sensibilite}")
    output_video_folder = os.path.join(dst_folder, "output_videos")
    os.makedirs(output_video_folder, exist_ok=True)

        
    Si = analyze_sobol(output_video_folder, nb_threads, seed, root_path, analyse_sensibilite, param_bounds, xml_path, dst_folder, N, dict_corresp_name_type, dict_corres_microenv)
    print(f"Sobol Indices: {Si}")

    #Enregistrer les résultats dans un fichier texte 
    with open(os.path.join(dst_folder, "sobol_indices_result.txt"), "w", encoding="utf-8") as f:
        f.write(f"Sobol Indices for {analyse_sensibilite}:\n {Si}")

if __name__ == "__main__":
    path_param = sys.argv[1]
    print("Parameter path", path_param)
    with open(sys.argv[1], "r") as f:
        params = json.load(f)
    run_main_analysis(params)

    """
    # ***************************************
    # !! A changer avant chaque analyse : 
    analyse_sensibilite = "tumor_persistance"
    root_path = "C://Users/elisa/" 
    result_folder_path =  "C://Users/elisa/results"
    xml_path =  "C://Users/elisa/PhysiCell/config/Physicell_settings.xml"
    nb_threads = 4 
    seed = 19 
    N = 2

    param_bounds = []
    if analyse_sensibilite == "descent_time" : 
        attachment_rate_bound = [0, 1]
        cancer_motility_speed_bound = [0, 1]
        transformation_rate_mes_bound = [0, 1] 
        cell_adhesion_affinity_LM_EP_bound = [0, 1] 
        cell_adhesion_affinity_LB_conj_bound = [0, 1]
        param_bounds = [attachment_rate_bound, cancer_motility_speed_bound, transformation_rate_mes_bound, cell_adhesion_affinity_LM_EP_bound, cell_adhesion_affinity_LB_conj_bound]
    elif analyse_sensibilite == "tumor_persistance" : 
        motility_speed_t_cell_bound = [0, 1] 
        division_duration_cancer_bound = [0, 1]
        param_bounds = [motility_speed_t_cell_bound, division_duration_cancer_bound]
    else : 
        raise ValueError("L'analyse de sensibilité demandée est mal écrite ou n'est pas dans les propositions.")


    # ***************************************

    dst_folder = os.path.join(root_path, "Results_PhysiCell", f"sensitivity_analysis_{analyse_sensibilite}")
    output_video_folder = os.path.join(dst_folder, "output_videos")
    os.makedirs(output_video_folder, exist_ok=True)

        
    Si = analyze_sobol(output_video_folder, nb_threads, seed, root_path, analyse_sensibilite, param_bounds)
    print(f"Sobol Indices: {Si}")

    #Enregistrer les résultats dans un fichier texte 
    with open(os.path.join(dst_folder, "sobol_indices_result.txt"), "w", encoding="utf-8") as f:
        f.write(f"Sobol Indices for {analyse_sensibilite}:\n {Si}") """
