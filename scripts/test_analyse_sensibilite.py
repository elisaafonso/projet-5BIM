from SALib.analyze.sobol import analyze
from SALib.sample.sobol import sample
import numpy as np
import xml.etree.ElementTree as ET
from functions import explore_tree, dict_to_xml, list_path_folder, computation_area_over_time, list_mat_files, compute_time_ratio, get_result_mat_persistance, get_cell_type
import os 
import subprocess
import pandas as pd
import re
import sys
import json

def define_set_param(analyse_sensibilite, param_bounds, N):
    """
    Echantillonage de l'espace des paramètres pour l'analyse de sensibilité Sobol.

    Input : 
    - analyse_sensibilite : str, type d'analyse de sensibilité ("descent_time" ou "tumor_persistance")
    - param_bounds : list of list, bornes des paramètres à échantillonner
    - N : int, nombre d'échantillons à générer (doit être une puissance de 2)

    Output : 
    problem : dict, définition du problème pour l'analyse de sensibilité
    param_values : matrice de paramètres échantillonnés
    """
    print("Defining parameter space and sampling...")

    if analyse_sensibilite == "descent_time1" : 
        problem = {
            'num_vars': 4,
            'names': ['attachment_rate', 'cancer_motility_speed', 'transformation_rate_mes', 'secretion_mmp_factor'],
            'bounds': [param_bounds[0], param_bounds[1], param_bounds[2], param_bounds[3]] 
        }
    elif analyse_sensibilite == "tumor_persistance1" :
        problem = {
            'num_vars': 4,
            'names': ['motility_speed_t_cell', 'division_duration_cancer', 'death_cancer', 'damage_attack_rate'],
            'bounds': [param_bounds[0], param_bounds[1], param_bounds[2], param_bounds[3]]
        }
    elif analyse_sensibilite == "tumor_persistance2" :
        problem = {
            'num_vars': 3,
            'names': ['motility_speed_t_cell', 'division_duration_cancer',  'damage_attack_rate'],
            'bounds': [param_bounds[0], param_bounds[1], param_bounds[2]]
        }
    elif analyse_sensibilite == "descent_time2" : 
        problem = {
            'num_vars': 3,
            'names': ['attachment_rate', 'cancer_motility_speed', 'secretion_mmp_factor'],
            'bounds': [param_bounds[0], param_bounds[1], param_bounds[2]] 
        }

    param_values = sample(problem, N, seed = 1) #Génère N*(2+D) jeux de paramètres avec D le nombre de paramètres et N un multiple de 2 fourni en argument
    return problem, param_values

def define_settings(parameters_to_change, nb_threads, seed, analyse_sensibilite, xml_path, output_path_i, dict_corresp_name_type, dict_corres_microenv):
    """
    Modification des paramètres dans le fichier XML de configuration de PhysiCell.
    Input : 
    - parameters_to_change : list, liste des nouveaux paramètres à insérer dans le fichier XML
    - nb_threads : int, nombre de coeurs à utiliser
    - seed : int
    - analyse_sensibilite : str, type d'analyse de sensibilité ("descent_time" ou "tumor_persistance")
    - xml_path : str, chemin vers le fichier XML de configuration
    - output_path_i : str, chemin vers le dossier de sortie pour cette simulation
    - dict_corresp_name_type : dict, dictionnaire de correspondance entre les noms des types cellulaires et leurs IDs
    - dict_corres_microenv : dict, dictionnaire de correspondance entre des microenvironnements et leurs IDs

    Output : 
    - dict_xml : dict, dictionnaire représentant le fichier XML avec les nouveaux paramètres
    """
    print("Defining XML settings for the simulation...")

    #lecture du fichier XML
    tree = ET.parse(xml_path)
    root = tree.getroot()
    dict_xml = {root.tag: explore_tree(root)}

    #changement des paramètres 

    ################## DESCENT TIME ####################
    if analyse_sensibilite == "descent_time1" :
        # Attachment rate : 0 et 10
        attachment_rate = parameters_to_change[0]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["membrane"]]["phenotype"]["mechanics"]["attachment_rate"]["text_explanation"] = attachment_rate

        #Motility speed cancer cells : 0.01 et 1
        cancer_motility_speed = parameters_to_change[1]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer"]]["phenotype"]["motility"]["speed"]["text_explanation"] = cancer_motility_speed
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["motility"]["speed"]["text_explanation"] = cancer_motility_speed

        # Mes transition : 0.00001 et 0.001
        transformation_rate_mes = parameters_to_change[2]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer"]]["phenotype"]["cell_transformations"]["transformation_rates"]["transformation_rate"][dict_corresp_name_type["cancer_mes"]]["text_explanation"] = transformation_rate_mes
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["cell_transformations"]["transformation_rates"]["transformation_rate"][dict_corresp_name_type["cancer"]]["text_explanation"] = transformation_rate_mes

        # Secretion cancer mes -- 0 et 50
        secretion_mmp_factor = parameters_to_change[3]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["secretion"]["substrate"][dict_corres_microenv["mmp_factor"]]["secretion_rate"]["text_explanation"] = secretion_mmp_factor

    ################## TUMOR PERSISTANCE ########################################
    elif analyse_sensibilite == "tumor_persistance1" :
        # Vitesse motilité T-Cell - 0.01 à 1
        motility_speed_t_cell = parameters_to_change[0]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["TCell"]]["phenotype"]["motility"]["speed"]["text_explanation"] = motility_speed_t_cell

        # Division cellules cancéreuses mes - 1440 à 4320
        division_duration_cancer = parameters_to_change[1]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["cycle"]["phase_durations"]["duration"]["text_explanation"] = division_duration_cancer
        
        #mort des cellules cancéreuses - 0.1 e-5 à 1 e-5
        death_cancer = parameters_to_change[2]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["death"]["model"][0]["death_rate"]["text_explanation"] = death_cancer

        #damage attack rate - 0.2 à 2
        damage_attack_rate = parameters_to_change[3]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["TCell"]]["phenotype"]["cell_interactions"]["attack_damage_rate"]["text_explanation"] = damage_attack_rate


    elif analyse_sensibilite == "tumor_persistance2" :
        # Vitesse motilité T-Cell - 0.01 à 1
        motility_speed_t_cell = parameters_to_change[0]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["TCell"]]["phenotype"]["motility"]["speed"]["text_explanation"] = motility_speed_t_cell

        # Division cellules cancéreuses mes - 1440 à 4320
        division_duration_cancer = parameters_to_change[1]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["cycle"]["phase_durations"]["duration"]["text_explanation"] = division_duration_cancer
        
        #damage attack rate - 0.5 à 2
        damage_attack_rate = parameters_to_change[2]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["TCell"]]["phenotype"]["cell_interactions"]["attack_damage_rate"]["text_explanation"] = damage_attack_rate

        #mort des cellules cancéreuses fixée - 0.31667e-05
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["death"]["model"][0]["death_rate"]["text_explanation"] = 0.31667e-05

        #change simulation time 
        #dict_xml["PhysiCell_settings"]["overall"]["max_time"]["text_explanation"] = 10000 #10000 minutes = ~7 jours

    elif analyse_sensibilite == "descent_time2" :
        # Attachment rate : 0 et 10
        attachment_rate = parameters_to_change[0]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["membrane"]]["phenotype"]["mechanics"]["attachment_rate"]["text_explanation"] = attachment_rate

        #Motility speed cancer cells : 0.01 et 1
        cancer_motility_speed = parameters_to_change[1]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["motility"]["speed"]["text_explanation"] = cancer_motility_speed

        # Secretion cancer mes -- 0 et 50
        secretion_mmp_factor = parameters_to_change[2]
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["secretion"]["substrate"][dict_corres_microenv["mmp_factor"]]["secretion_rate"]["text_explanation"] = secretion_mmp_factor

        # Mes transition fixée à 0
        transformation_rate_mes = 0
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer"]]["phenotype"]["cell_transformations"]["transformation_rates"]["transformation_rate"][dict_corresp_name_type["cancer_mes"]]["text_explanation"] = transformation_rate_mes
        dict_xml["PhysiCell_settings"]["cell_definitions"]["cell_definition"][dict_corresp_name_type["cancer_mes"]]["phenotype"]["cell_transformations"]["transformation_rates"]["transformation_rate"][dict_corresp_name_type["cancer"]]["text_explanation"] = transformation_rate_mes

        #change simulation time 
        #dict_xml["PhysiCell_settings"]["overall"]["max_time"]["text_explanation"] = 10000 #10000 minutes = ~7 jours
        
    ################## CHANGER LE NOMBRE DE COEUR ET LA SEED  ####################
    dict_xml["PhysiCell_settings"]["parallel"]["omp_num_threads"]["text_explanation"] = nb_threads
    dict_xml["PhysiCell_settings"]["options"]["random_seed"]["text_explanation"] = seed
    
    #define output_path
    dict_xml["PhysiCell_settings"]["save"]["folder"]["text_explanation"] = output_path_i #output

    return dict_xml #retourne un dictionnaire XML avec les nouveaux paramètres

def get_physicell_output(param_values, nb_threads, seed, root_path, analyse_sensibilite, xml_path, dst_folder, dict_corresp_name_type, dict_corres_microenv):
    """
    Exécute les simulations PhysiCell pour chaque jeu de paramètres et stocke les métriques.
    Input :
    - param_values : matrice de paramètres échantillonnés
    - nb_threads : int, nombre de coeurs à utiliser
    - seed : int
    - root_path : str, chemin vers le répertoire racine de PhysiCell
    - analyse_sensibilite : str, type d'analyse de sensibilité ("descent_time" ou "tumor_persistance")
    - xml_path : str, chemin vers le fichier XML de configuration
    - dst_folder : str, chemin vers le dossier de sortie pour les résultats
    - dict_corresp_name_type : dict, dictionnaire de correspondance entre les noms des types cellulaires et leurs IDs
    - dict_corres_microenv : dict, dictionnaire de correspondance entre des microenvironnements et leurs IDs

    Output :
    - output_storage : np.array, vecteur des métriques calculées pour chaque simulation
    """
    output_storage = np.zeros(param_values.shape[0]) #initialisation du vecteur de sortie
    for i in range (len(param_values)): 
        print(f"Running simulation {i+1}...")

        #créer le dossier output pour le pas de temps i
        output_path_i = os.path.normpath(os.path.join(dst_folder, f"output_{i}"))
        os.makedirs(output_path_i, exist_ok=True)

        #modification du fichier de configuration XML
        parameters_to_change = param_values[i]
        dict_xml = define_settings(parameters_to_change, nb_threads, seed, analyse_sensibilite, xml_path, output_path_i, dict_corresp_name_type, dict_corres_microenv)
        xml_element = dict_to_xml(dict_xml)
        tree = ET.ElementTree(xml_element)
        tree.write(os.path.join(root_path, f"PhysiCell/config/PhysiCell_settings_{i}.xml"), encoding="utf-8", xml_declaration=True) #path à changer

        #vider le dossier output_{i} avant chaque simulation
        process0 = subprocess.run(
        ["rm", "-rf", "*"],
        capture_output=True,
        text=True,
        cwd=output_path_i 
        )

        print("Configuration file updated. Starting simulation...")
        # Run .exe file 
        exe_path = os.path.join(root_path, "PhysiCell", "project")
        xml_path_i = os.path.join(root_path, "PhysiCell", "config", f"PhysiCell_settings_{i}.xml")

        process1 = subprocess.run(
            [exe_path, xml_path_i],
            capture_output=True,
            text=True,
            cwd=os.path.join(root_path, "PhysiCell")
        )

        # Affichier les outputs et erreurs
        print("Output:")
        print(process1.stdout) 
        print("Error:")
        print(process1.stderr)

        final_time = int(float(dict_xml["PhysiCell_settings"]["overall"]["max_time"]["text_explanation"]))
        interval_time = int(float(dict_xml["PhysiCell_settings"]["save"]["full_data"]["interval"]["text_explanation"]))

        if analyse_sensibilite == "tumor_persistance1" or analyse_sensibilite == "tumor_persistance2":  
            id_cancer = dict_corresp_name_type["cancer"]
            id_cancer_mes = dict_corresp_name_type["cancer_mes"]
            files_by_timestep = list_mat_files(output_path_i, final_time, interval_time)
            result_mat = get_result_mat_persistance(files_by_timestep, output_path_i, id_cancer, id_cancer_mes)
            dt = int(dict_xml["PhysiCell_settings"]["save"]["full_data"]["interval"]["text_explanation"])/60.0  # conversion en heures
            metric = computation_area_over_time(result_mat, dt)
        
        elif analyse_sensibilite == "descent_time1" or analyse_sensibilite == "descent_time2": 
            id_cancer_cell = dict_corresp_name_type["cancer"]
            id_cancer_cell_mes = dict_corresp_name_type["cancer_mes"]
            id_connective_tissue = dict_corresp_name_type["conjonctif"]
            id_fibroblast = dict_corresp_name_type["CAF"]
            files_by_timestep = list_path_folder(output_path_i, final_time, interval_time)
            ratio = compute_time_ratio(files_by_timestep, output_path_i, id_cancer_cell, id_cancer_cell_mes, id_connective_tissue, id_fibroblast)  
            metric = ratio

        output_storage[i] = metric #stockage de la métrique
    return output_storage

def analyze_sobol(nb_threads, seed, root_path, analyse_sensibilite, param_bounds, xml_path, dst_folder, N, dict_corresp_name_type, dict_corres_microenv):
    """
    Effectue l'analyse de sensibilité Sobol en exécutant les simulations PhysiCell et en calculant les indices de Sobol.
    Input :
    - nb_threads : int, nombre de coeurs à utiliser
    - seed : int
    - root_path : str, chemin vers le répertoire racine de PhysiCell
    - analyse_sensibilite : str, type d'analyse de sensibilité ("descent_time" ou "tumor_persistance")
    - param_bounds : list of list, bornes des paramètres à échantillonner
    - xml_path : str, chemin vers le fichier XML de configuration
    - dst_folder : str, chemin vers le dossier de sortie pour les résultats
    - N : int, nombre d'échantillons à générer (doit être une puissance de 2)
    - dict_corresp_name_type : dict, dictionnaire de correspondance entre les noms des types cellulaires et leurs IDs
    - dict_corres_microenv : dict, dictionnaire de correspondance entre des microenvironnements et leurs IDs

    Output : 
    - Si : dict, indices de Sobol calculés
    """
    problem, param_values = define_set_param(analyse_sensibilite, param_bounds, N)
    
    #Save parameters values to a CSV file
    df = pd.DataFrame(param_values) 
    df.to_csv(os.path.join(dst_folder, "param_values.csv"), index=False)

    Y = get_physicell_output(param_values, nb_threads, seed, root_path, analyse_sensibilite, xml_path, dst_folder, dict_corresp_name_type, dict_corres_microenv)
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
    dst_folder = os.path.join(result_folder_path, "Results_PhysiCell", f"sensitivity_analysis_{analyse_sensibilite}")
    os.makedirs(dst_folder, exist_ok=True)

    Si = analyze_sobol(nb_threads, seed, root_path, analyse_sensibilite, param_bounds, xml_path, dst_folder, N, dict_corresp_name_type, dict_corres_microenv)
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

    
