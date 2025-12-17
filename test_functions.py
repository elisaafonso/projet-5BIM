from functions import explore_tree, dict_to_xml, list_path_folder, computation_area_over_time, clustering_cells
import xml.etree.ElementTree as ET
import os 
import re

root_path = "C://Users/elisa/"
projet_path = "_COURS_5A/projet-5BIM/"

def test_xml_reading():
    tree = ET.parse(os.path.join(root_path, projet_path, "test/PhysiCell_settings.xml"))
    root = tree.getroot()
    dict_xml = {root.tag: explore_tree(root)}
    xml_element = dict_to_xml(dict_xml)
    tree = ET.ElementTree(xml_element)
    tree.write(os.path.join(root_path, projet_path, "test/PhysiCell_settings_generated.xml"), encoding="utf-8", xml_declaration=True) 
    tree_generated = ET.parse(os.path.join(root_path, projet_path, "test/PhysiCell_settings_generated.xml"))
    root_gen = tree_generated.getroot()
    dict_xml_to_verify = {root_gen.tag: explore_tree(root_gen)}
    assert dict_xml == dict_xml_to_verify

def test_list_path_folder(): 
    output_path = os.path.join(root_path, "Physicell/output")
    tree = ET.parse(os.path.join(root_path, projet_path, "test/PhysiCell_settings.xml"))
    root = tree.getroot()
    dict_xml = {root.tag: explore_tree(root)}
    final_time = int(dict_xml["PhysiCell_settings"]["overall"]["max_time"]["text_explanation"])
    interval_time = int(dict_xml["PhysiCell_settings"]["save"]["full_data"]["interval"]["text_explanation"])
    nb_timestep = final_time // interval_time
    nb_files = nb_timestep + 1
    files = list_path_folder(output_path)
    assert len(files) == nb_files

def test_computation_area_over_time():
    expected_area = 47600
    dt = 2 #hours
    clusters0 = [{"1": 1000, "2": 2000}, {"3": 1500, "4": 2500}]
    clusters1 = [{"5": 1200, "6": 2200}, {"7": 1800, "8": 2800}]
    clusters2 = [{"9": 1400, "10": 2400}, {"11": 2000, "12": 3000}]
    result_mat = {"0" : [clusters0, 488, 28900], "1" : [clusters1, 490, 28900], "2" : [clusters2, 495, 28900]} #mock data
    area_over_time = computation_area_over_time(result_mat, dt)
    assert area_over_time == expected_area

def test_clustering_cells():
    pass 
