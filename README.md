# Modélisation multi‑échelle de la dynamique des cancers de la bouche et de leurs interactions avec le micro‑environnement buccal

## 1. Description de l’analyse de sensibilité

Un modèle agent a été implémenté sur PhysiCell avec un agent par type cellulaire de la bouche :

- Cellules de l’épithélium basal `epi_basal`
- Cellules de l’épithélium intermédiaire `epi_inter`
- Cellules de l’épithélium supérieur `epi_sup`
- Cellules cancéreuses `cancer`
- Cellules cancéreuses mésenchymateuses (perte d’adhésion avec les autres agents) `cancer_mes`
- Cellules T `TCell`
- CAF – Cancer Associated Fibroblast `CAF`
- Membrane basale `membrane`

Deux analyses de sensibilité globale ont été réalisées à l’aide de la méthode de Sobol implémentée dans le package Python SALib : [SALib](https://salib.readthedocs.io/en/latest/api.html).  
L’objectif de ces analyses était d’identifier les facteurs du micro‑environnement tumoral qui favorisent la persistance d’une tumeur et de déterminer quels paramètres ont le plus d’impact.  
Une tumeur est considérée persistante lorsqu’elle descend dans le tissu conjonctif et échappe au phénomène de tapis roulant de l’épithélium.

---

## 1.1 Descente de la tumeur dans le tissu conjonctif

La première analyse avait pour objectif d’évaluer la descente de la tumeur dans le tissu conjonctif.  
Les fichiers d’initialisation de PhysiCell sont disponibles dans le dossier `/initialisation_files/analyse_sens_descent_time/`.

Dans les conditions d’initialisation, une seule cellule cancéreuse (non mésenchymateuse) est présente. Elle ne peut ni se diviser, ni mourir.  
Le ratio du temps de descente de cette cellule dans le tissu conjonctif sur le temps total de simulation est calculé :

\[
descent_time = \frac{t*{\text{descent in conj}}}{t*{\text{tot}}}
\]

- Si le ratio est égal à 1, la cellule a été éjectée par le tapis roulant ou n’a jamais traversé la lame basale.
- Si le ratio est égal à 0, la cellule est passée dans le tissu conjonctif dès le premier pas de temps.

Cette métrique permet d’évaluer la facilité avec laquelle la cellule cancéreuse s’infiltre dans le tissu conjonctif.

Quatre paramètres ont été choisis pour cette analyse de sensibilité :

### • Vitesse de migration des agents `cancer` et `cancer_mes` (\(s\_{mot}\))

Cette vitesse fait partie des paramètres définissant la motilité d’un agent.  
Le biais de migration \(d\_{bias}\) a été fixé à 0.5, ce qui confère à l’agent une motilité **semi‑déterministe** en réponse à un stimulus chimique, ici le `CAF_chemotaxis` qui attire les cellules cancéreuses.

Modifier la vitesse de migration influence donc la rapidité avec laquelle l’agent répond à la chimio‑attraction (Ghaffarizadeh et al., 2018)[1].

La plage de variation choisie est 0.01 à 1 (valeur par défaut).  
À \(s*{mot} = 1\), l’agent traverse la membrane très rapidement, alors qu’à \(s*{mot} = 0.01\), il le fait beaucoup plus rarement.

### • Taux d’attachement des agents de la `membrane` entre eux

Les agents de la membrane sont attachés entre eux.  
Ce paramètre a été fixé arbitrairement entre 0 (pas attachés) et 10 (très attachés).

### • Transition entre `cancer` et `cancer_mes` (et inversement)

Les cellules cancéreuses ont une certaine probabilité par minute de devenir mésenchymateuses (ou de redevenir adhérentes).  
Cette probabilité a été fixée entre 0.00001 et 0.001.  
Au‑delà de 0.001, l’agent changeait trop souvent de type (environ toutes les heures).

### • Sécrétion du facteur MMP par les cellules cancéreuses

Ce paramètre a été fixé arbitrairement entre 0 et 50.  
Le facteur MMP, sécrété par les cellules cancéreuses mésenchymateuses, augmente la probabilité de dégradation des agents de la membrane au contact.

---

## 1.2 Persistance de la tumeur dans le tissu conjonctif

La deuxième analyse de sensibilité visait à évaluer la persistance de la tumeur dans le tissu conjonctif en fonction du micro‑environnement présent.

Une fois la cellule cancéreuse entrée dans le tissu conjonctif, il a été supposé que la tumeur n’interagissait pas avec la matrice extracellulaire (tissu conjonctif) ni avec les CAF.  
Dans les conditions initiales, seules des cellules tumorales et des cellules T ont été ajoutées.  
L’hypothèse est que ce sont principalement ces deux types cellulaires, ainsi que certains paramètres les concernant, qui influencent la persistance tumorale au cours du temps.

La métrique choisie est le volume tumoral cumulé au cours du temps :

\[
volume_over_time = \int_0^T \text{volume des cellules tumorales}(t)\, dt
\]

Quatre paramètres ont également été identifiés pour cette analyse.

### • Vitesse de migration des `TCell`

Les cellules cancéreuses émettent un `cancer_factor` qui attire les cellules T.  
Comme dans la partie 1.1, modifier la vitesse de migration influence la force d’attraction vers le stimulus chimique, la sensibilité au chimiotactisme ayant été fixée.

Il a été observé que modifier la vitesse de migration \(s\_{mot}\) avait un impact plus important que modifier la sensibilité au chimiotactisme.  
La vitesse a donc été fixée entre 0.01 et 1.

### • Division des cellules cancéreuses - 0.000001 à 0.00001

### • Mort des cellules cancéreuses - 0.1 e-5 à 1 e-5

### • Damage attack rate des TCell vers les cellules cancéreuses - 0.2 à 2

## 2. Description du repo GitHub

En fonction de votre système d'exploitation (Windows ou Unix), une modification devra être réalisée. Elle est décrite dans la première partie.

## Dans le fichier `test_analyse_sensibilite.py`

- Les paramètres à mettre dans l'analyse de sensibilité sont à indiquer dans la fonction `define_set_param`.
- Leur emplacement dans le fichier xml est défini dans la fonction `define_settings`.

### Pour un PC Unix

Dans la fonction `get_physicell_output`, changer la commande bash dans process1 :

```python
# Run .exe file
        process1 = subprocess.run(
        [os.path.join(root_path, "PhysiCell/project")], # !! To change
        capture_output=True,
        text=True,
        cwd=os.path.join(root_path, "PhysiCell")
        )
```

---

## Fichier XML

- Le fichier xml est récupéré du dossier `config/` de PhysiCell.

---

## Chemin racine

- Le chemin vers la racine où se trouve le dossier PhysiCell est à indiquer/modifier :
  ```python
  root_path = "C://Users/elisa/"
  ```

---

## Fichier functions.py

Les différentes fonctions sont codées dans le fichier : functions.py.

---

## Fichier requirements.txt

Les différents packages à installer se trouvent dans le fichier : requirements.txt.

---

## Pour lancer le script

- Créer l'environnement virtuel : `python -m venv AS_env` et l'activer
- Installer les dépendances : `python -m pip install -r requirements.txt`
- Lancer l'analyse de sensibilité : `python test_analyse_sensibilite.py` après avoir lancer un environnement virtuel.

rajouter make load PROJ=5_bim

## 3. Bibliography

[1] Ghaffarizadeh, A., Heiland, R., Friedman, S. H., Mumenthaler, S. M., & Macklin, P. (2018). PhysiCell : An open source physics-based cell simulator for 3-D multicellular systems. PLoS Computational Biology, 14(2), e1005991. https://doi.org/10.1371/journal.pcbi.1005991

$$
$$
