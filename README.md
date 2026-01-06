# Analyse de sensibilité

Ce code fonctionne sur un PC Windows, quelques modifications sont à réaliser sur un PC Unix (décrites dans la première partie).

L'analyse de sensibilité est réalisée à l'aide de la méthode de Sobol implémentée dans le package python SALib : [SALib](https://salib.readthedocs.io/en/latest/api.html).

---

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
