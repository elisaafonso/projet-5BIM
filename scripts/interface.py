import tkinter as tk
from tkinter import messagebox
import xml.etree.ElementTree as ET
from test_analyse_sensibilite import run_main_analysis
import threading
import os

root = tk.Tk()
root.title("Paramètres Physicell")
root.geometry("650x650")

path_valid = False
ana_sens = None

# -----------------------------
# 1. Validation des chemins
# -----------------------------
def validate_path():
    global path_valid

    xml_path = entry_chemin_xml.get()
    physicell_path = entry_chemin_Physicell.get()
    results_path = entry_chemin_results.get()

    if os.path.isfile(xml_path) and os.path.isdir(physicell_path) and os.path.isdir(results_path):
        try:
            ET.parse(xml_path)
        except Exception:
            messagebox.showerror("Erreur", "Le fichier XML est invalide.")
            path_valid = False
            return

        messagebox.showinfo("Succès", "Chemins valides.")
        path_valid = True
        show_analysis_choice()
    else:
        messagebox.showerror("Erreur", "Un ou plusieurs chemins sont invalides.")
        path_valid = False


# -----------------------------
# 2. Choix du type d'analyse
# -----------------------------
def show_analysis_choice():
    frame_choice.pack(pady=10)


def ana_sens_test():
    global ana_sens

    if var_descent_time.get() and var_tumor_persistance.get():
        messagebox.showerror("Erreur", "Veuillez sélectionner un seul type d'analyse.")
        return

    if not var_descent_time.get() and not var_tumor_persistance.get():
        messagebox.showerror("Erreur", "Veuillez sélectionner un type d'analyse.")
        return

    ana_sens = "descent_time" if var_descent_time.get() else "tumor_persistance"
    show_parameter_inputs()


# -----------------------------
# 3. Affichage des paramètres
# -----------------------------
def show_parameter_inputs():
    frame_descent_time.pack_forget()
    frame_tumor_persistance.pack_forget()

    if ana_sens == "descent_time":
        frame_descent_time.pack(pady=10)
    else:
        frame_tumor_persistance.pack(pady=10)


# -----------------------------
# 4. Lecture des paramètres
# -----------------------------
def get_parameters():
    try:
        if ana_sens == "descent_time":
            param_bounds = [
                [float(entry_attachment_rate_min.get()), float(entry_attachment_rate_max.get())],
                [float(entry_cancer_motility_speed_min.get()), float(entry_cancer_motility_speed_max.get())],
                [float(entry_transformation_rate_mes_min.get()), float(entry_transformation_rate_mes_max.get())],
                [float(entry_cell_adhesion_affinity_LM_EP_min.get()), float(entry_cell_adhesion_affinity_LM_EP_max.get())],
                [float(entry_cell_adhesion_affinity_LB_conj_min.get()), float(entry_cell_adhesion_affinity_LB_conj_max.get())]
            ]

            seed = int(entry_seed_dt.get())
            threads = int(entry_threads_dt.get())
            nb_sample = int(entry_nb_sample_to_generate_dt.get())

        else:
            param_bounds = [
                [float(entry_motility_speed_t_cell_min.get()), float(entry_motility_speed_t_cell_max.get())],
                [float(entry_division_duration_cancer_min.get()), float(entry_division_duration_cancer_max.get())]
            ]

            seed = int(entry_seed_tp.get())
            threads = int(entry_threads_tp.get())
            nb_sample = int(entry_nb_sample_to_generate_tp.get())

        # Vérification puissance de 2
        if nb_sample <= 0 or (nb_sample & (nb_sample - 1)) != 0:
            messagebox.showerror("Erreur", "Le nombre d'échantillons doit être une puissance de 2 (ex : 8, 16, 32, 64...).")
            return None

    except ValueError:
        messagebox.showerror("Erreur", "Veuillez entrer des valeurs numériques valides.")
        return None

    return {
        "analyse_sensibilite": ana_sens,
        "param_bounds": param_bounds,
        "seed": seed,
        "nb_threads": threads,
        "xml_path": entry_chemin_xml.get(),
        "physicell_path": entry_chemin_Physicell.get(),
        "results_path": entry_chemin_results.get(),
        "nb_sample_to_generate": nb_sample
    }


# -----------------------------
# 5. Lancer l'analyse
# -----------------------------
def run_analysis():
    params = get_parameters()
    if params is None:
        return

    t = threading.Thread(target=run_main_analysis, args=(params,))
    t.start()

    messagebox.showinfo("Info", "L'analyse démarre.\nLa fenêtre va se fermer dans 5 secondes.")
    root.after(5000, root.destroy)


# -----------------------------
# Widgets : Chemins
# -----------------------------
tk.Label(root, text="Chemin du fichier XML").pack()
entry_chemin_xml = tk.Entry(root, width=50)
entry_chemin_xml.insert(0, "PhysiCell_settings.xml")
entry_chemin_xml.pack()

tk.Label(root, text="Chemin du dossier PhysiCell").pack()
entry_chemin_Physicell = tk.Entry(root, width=50)
entry_chemin_Physicell.insert(0, "C://Users/$USER/")
entry_chemin_Physicell.pack()

tk.Label(root, text="Chemin du dossier résultats").pack()
entry_chemin_results = tk.Entry(root, width=50)
entry_chemin_results.insert(0, "C://Users/$USER/Results/")
entry_chemin_results.pack()

tk.Button(root, text="Valider les chemins", command=validate_path).pack(pady=10)


# -----------------------------
# Choix analyse
# -----------------------------
frame_choice = tk.Frame(root)
tk.Label(frame_choice, text="Choisir le type d'analyse :").pack()

var_descent_time = tk.BooleanVar()
var_tumor_persistance = tk.BooleanVar()

tk.Checkbutton(frame_choice, text="Descent Time", variable=var_descent_time).pack()
tk.Checkbutton(frame_choice, text="Tumor Persistance", variable=var_tumor_persistance).pack()

tk.Button(frame_choice, text="Valider", command=ana_sens_test).pack(pady=10)


# -----------------------------
# Fonction utilitaire
# -----------------------------
def add_param_row(parent, label, default_min, default_max):
    row = tk.Frame(parent)
    tk.Label(row, text=label, width=25, anchor="w").pack(side="left")
    entry_min = tk.Entry(row, width=10)
    entry_min.insert(0, default_min)
    entry_min.pack(side="left", padx=5)
    entry_max = tk.Entry(row, width=10)
    entry_max.insert(0, default_max)
    entry_max.pack(side="left", padx=5)
    row.pack(pady=2)
    return entry_min, entry_max


# -----------------------------
# Paramètres descent_time
# -----------------------------
frame_descent_time = tk.Frame(root)
tk.Label(frame_descent_time, text="Paramètres : descent_time").pack()

entry_attachment_rate_min, entry_attachment_rate_max = add_param_row(frame_descent_time, "Attachment rate", "0.01", "0.1")
entry_cancer_motility_speed_min, entry_cancer_motility_speed_max = add_param_row(frame_descent_time, "Cancer motility speed", "0.1", "2.0")
entry_transformation_rate_mes_min, entry_transformation_rate_mes_max = add_param_row(frame_descent_time, "Transformation rate MES", "0.001", "0.01")
entry_cell_adhesion_affinity_LM_EP_min, entry_cell_adhesion_affinity_LM_EP_max = add_param_row(frame_descent_time, "Adhesion LM-EP", "0.1", "1.0")
entry_cell_adhesion_affinity_LB_conj_min, entry_cell_adhesion_affinity_LB_conj_max = add_param_row(frame_descent_time, "Adhesion LB-conj", "0.1", "1.0")

# Seed / Threads / Samples
frame_extra_dt = tk.Frame(frame_descent_time)
tk.Label(frame_extra_dt, text="Random seed", width=25, anchor="w").pack(side="left")
entry_seed_dt = tk.Entry(frame_extra_dt, width=10)
entry_seed_dt.insert(0, "19")
entry_seed_dt.pack(side="left")
frame_extra_dt.pack()

frame_extra2_dt = tk.Frame(frame_descent_time)
tk.Label(frame_extra2_dt, text="Core number", width=25, anchor="w").pack(side="left")
entry_threads_dt = tk.Entry(frame_extra2_dt, width=10)
entry_threads_dt.insert(0, "4")
entry_threads_dt.pack(side="left")
frame_extra2_dt.pack()

frame_extra3_dt = tk.Frame(frame_descent_time)
tk.Label(frame_extra3_dt, text="Nb samples", width=25, anchor="w").pack(side="left")
entry_nb_sample_to_generate_dt = tk.Entry(frame_extra3_dt, width=10)
entry_nb_sample_to_generate_dt.insert(0, "8")
entry_nb_sample_to_generate_dt.pack(side="left")
frame_extra3_dt.pack()

tk.Button(frame_descent_time, text="Lancer l'analyse", command=run_analysis).pack(pady=10)


# -----------------------------
# Paramètres tumor_persistance
# -----------------------------
frame_tumor_persistance = tk.Frame(root)
tk.Label(frame_tumor_persistance, text="Paramètres : tumor_persistance").pack()

entry_motility_speed_t_cell_min, entry_motility_speed_t_cell_max = add_param_row(frame_tumor_persistance, "T-cell motility speed", "0.1", "2.0")
entry_division_duration_cancer_min, entry_division_duration_cancer_max = add_param_row(frame_tumor_persistance, "Cancer division duration", "600", "2000")

frame_extra_tp = tk.Frame(frame_tumor_persistance)
tk.Label(frame_extra_tp, text="Random seed", width=25, anchor="w").pack(side="left")
entry_seed_tp = tk.Entry(frame_extra_tp, width=10)
entry_seed_tp.insert(0, "19")
entry_seed_tp.pack(side="left")
frame_extra_tp.pack()

frame_extra2_tp = tk.Frame(frame_tumor_persistance)
tk.Label(frame_extra2_tp, text="Core number", width=25, anchor="w").pack(side="left")
entry_threads_tp = tk.Entry(frame_extra2_tp, width=10)
entry_threads_tp.insert(0, "4")
entry_threads_tp.pack(side="left")
frame_extra2_tp.pack()

frame_extra3_tp = tk.Frame(frame_tumor_persistance)
tk.Label(frame_extra3_tp, text="Nb samples", width=25, anchor="w").pack(side="left")
entry_nb_sample_to_generate_tp = tk.Entry(frame_extra3_tp, width=10)
entry_nb_sample_to_generate_tp.insert(0, "8")
entry_nb_sample_to_generate_tp.pack(side="left")
frame_extra3_tp.pack()

tk.Button(frame_tumor_persistance, text="Lancer l'analyse", command=run_analysis).pack(pady=10)


root.mainloop()