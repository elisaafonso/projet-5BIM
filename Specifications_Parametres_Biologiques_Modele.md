# **Spécifications des Paramètres Biologiques du Modèle**

Ce document présente une segmentation du fichier de configuration XML. Chaque bloc de code est accompagné d'une explication technique et d'une justification biologique sourcée, validant le choix des paramètres et des valeurs retenues.

**Notes préliminaires :**

1. Seuls les premiers éléments des parties 3 et 4 sont commentés de manière exhaustive avec un commentaire ligne par ligne sur le code du fichier xml. Ce choix est pris pour éviter de surcharger le document par des répétitions inutiles. Cependant, une signification et justification des autres types est tout de même rapportée.

2. Ce document constitue un support explicatif visant à établir un référentiel de paramètres cohérents avec la réalité biologique. S'il a servi de fondement théorique à la construction de nos modèles, il est important de noter que le fichier XML utilisé pour les simulations finales présente des ajustements par rapport à ce document visant à optimiser la stabilité de la simulation et la pertinence des observations. Cependant, tout paramètre a été finalement choisi afin de rentrer dans des plages de valeurs scientifiquement réalistes.

**Sommaire**

[Partie 1 : Cadre spatio-temporel](#partie-1-:-cadre-spatio-temporel)

[● Domaine de simulation](#domaine-de-simulation)

[● Durée et pas de temps](#durée-et-pas-de-temps)

[Partie 2 : paramètres techniques et globaux](#partie-2-:-paramètres-techniques-et-globaux)

[● Allocation des ressources de calculs (parallel)](<#allocation-des-ressources-de-calculs-(parallel)>)

[● Gestion des données de sorties (save)](<#gestion-des-données-de-sorties-(save)>)

[● Stabilité et reproductibilité](#stabilité-et-reproductibilité)

[Partie 3 : Microenvironnement](#partie-3-:-microenvironnement)

[● Facteur de croissance](#facteur-de-croissance)

[● Facteur cancéreux](#facteur-cancéreux)

[● Death factor](#death-factor)

[● Facteur de dégradation de la membrane basale (mmp_factor)](<#facteur-de-dégradation-de-la-membrane-basale-(mmp_factor)>)

[● Chimiotactisme induit par les CAF (CAF_chemotaxis)](<#chimiotactisme-induit-par-les-caf-(caf_chemotaxis)>)

[● Autres](#autres)

[Partie 4 : types cellulaires](#partie-4-:-types-cellulaires)

[● Epithélium basal](#epithélium-basal)

[● Epithelium intermédiaire](#epithelium-intermédiaire)

[● Epithelium supérieur](#epithelium-supérieur)

[● Cellules cancéreuses](#cellules-cancéreuses)

[● T Cell](#t-cell)

[● Membrane](#membrane)

[● Fibroblasts](#fibroblasts)

[● Cancer mésenchymateux](#cancer-mésenchymateux)

[● CAF](#caf)

[● Conjonctif](#conjonctif)

[Partie 5 : Initialisation et Paramètres](#partie-5-:-initialisation-et-paramètres)

[● Conditions initiales (cells.csv)](<#conditions-initiales-(cells.csv)>)

[● Les règles (cells_rules.csv)](<#les-règles-(cells_rules.csv)>)

# **Partie 1 : Cadre spatio-temporel** {#partie-1-:-cadre-spatio-temporel}

- ## **Domaine de simulation** {#domaine-de-simulation}

      \<domain\>
      \<x\_min\>-300\</x\_min\>
      \<x\_max\>300\</x\_max\>
      \<y\_min\>-300\</y\_min\>
      \<y\_max\>300\</y\_max\>
      \<z\_min\>-10\</z\_min\>
      \<z\_max\>10\</z\_max\>
      \<dx\>20\</dx\>
      \<dy\>20\</dy\>
      \<dz\>20\</dz\>
      \<use\_2D\>true\</use\_2D\>

  \</domain\>

**Signification :** Le domaine définit une fenêtre de simulation carrée de 600×600 microns, découpée en voxels de 20 µm de côté.

**Justification biologique :** Le choix de la fenêtre de simulation permet de contenir l'épaisseur complète d'un épithélium buccal non kératinisé, estimée en moyenne autour de 300µm. Une hauteur totale de 600 µm permet de modéliser confortablement l'épithélium et le chorion pour observer l'invasion. Les cellules de la cavité buccale ayant une dimension de quelques dizaines de microns réduite in vivo par compression, une résolution de 20 µm est le standard optimal.

**Source :** Stasio, Dario Di, et al. “Measurement of Oral Epithelial Thickness by Optical Coherence Tomography.” _Diagnostics_, vol. 9, no. 3, 6 Aug. 2019, p. 90, [https://doi.org/10.3390/diagnostics9030090](https://doi.org/10.3390/diagnostics9030090).

**Citation** : “_The healthy epithelium has a mean thickness of 335.59 ± 150.73 µm”_

- ## **Durée et pas de temps** {#durée-et-pas-de-temps}

\<overall\>  
 \<max_time units="min"\>57600.0\</max_time\>  
 \<time_units\>min\</time_units\>  
 \<space_units\>micron\</space_units\>

        \<dt\_diffusion units="min"\>0.01\</dt\_diffusion\>
        \<dt\_mechanics units="min"\>0.1\</dt\_mechanics\>
        \<dt\_phenotype units="min"\>6\</dt\_phenotype\>
    \</overall\>

**Signification :** Durée totale de 57 600 min (40 jours) avec des mises à jour biologiques toutes les 6 min et physiques toutes les 0.01 min.

**Justification :** Les 40 jours se divisent en une phase obligatoire de mise en place de l’épithélium à partir des cellules basales de 15 jours suivie de 25 jours d'observation de l'invasion tumorale. Ainsi, les 15 jours traduisent aussi le renouvellement de l’épithélium.

**Source** : Brizuela, M., & Winters, R. (2023, 8 mai). _Histology, oral mucosa_. StatPearls \- NCBI Bookshelf. https://www.ncbi.nlm.nih.gov/books/NBK572115/

**Citation** : “_The oral epithelial cells are frequently replaced by cell division, around each 14 to 21 days”_

# **Partie 2 : paramètres techniques et globaux** {#partie-2-:-paramètres-techniques-et-globaux}

- ## **Allocation des ressources de calculs (parallel)** {#allocation-des-ressources-de-calculs-(parallel)}

\<parallel\>  
 \<omp_num_threads\>1\</omp_num_threads\>  
\</parallel\>

**Signification :** Cette balise alloue le nombre de cœurs via la bibliothèque OpenMP pour effectuer des calculs simultanés. L'utilisation de plusieurs threads est indispensable pour accélérer le traitement des milliers d'agents cellulaires simulés, réduisant ainsi considérablement le temps réel de simulation, surtout pour les simulations longues avec beaucoup d’agents.

**Point technique :** Ce paramètre est aussi impliqué dans la reproductibilité des simulations c’est pourquoi en plus de fixer la random_seed il est indispensable de fixer le nombre de threads à 1 pour garantir la reproductibilité. En effet, en mode série (1 cœur), les opérations sont traitées les unes après les autres dans un ordre strict permettant d’avoir des résultats identiques. Cet ordre strict n’est pas assuré par la parallélisation lorsque plusieurs cœurs sont alloués.

**Justification :** Ainsi pour une simulation il est recommandé de se baser sur les capacités de la machine utilisée afin d’optimiser les temps de calcul. Cependant, si c’est la reproductibilité qui est recherchée il est conseillé cette fois-ci de fixer la valeur à 1\.

- ## **Gestion des données de sorties (save)** {#gestion-des-données-de-sorties-(save)}

\<save\>  
 \<folder\>output\</folder\>  
 \<full_data\>  
 \<interval units="min"\>120\</interval\>  
 \<enable\>true\</enable\>  
 \</full_data\>  
 \<SVG\>  
 \<interval units="min"\>120\</interval\>  
 \<enable\>true\</enable\>  
 \<plot_substrate enabled="false" limits="false"\>  
 \<substrate\>growth_factor\</substrate\>  
 \<min_conc /\>  
 \<max_conc /\>  
 \<colormap /\>  
 \</plot_substrate\>  
 \</SVG\>  
 \<legacy_data\>  
 \<enable\>false\</enable\>  
 \</legacy_data\>  
\</save\>

**Signification :** Ce bloc définit la manière d’enregistrement des résultats des simulations. Ainsi, les données sont stockées dans le fichier output avec des images SVG enregistrées toutes les 120 minutes. “plot substrate enabled” permet, si fixé sur True, de visualiser en fond la concentration d’un substrat.

- ## **Stabilité et reproductibilité** {#stabilité-et-reproductibilité}

\<options\>  
 \<legacy_random_points_on_sphere_in_divide\>false\</legacy_random_points_on_sphere_in_divide\>  
 \<virtual_wall_at_domain_edge\>true\</virtual_wall_at_domain_edge\>  
 \<disable_automated_spring_adhesions\>false\</disable_automated_spring_adhesions\>  
 \<random_seed\>19\</random_seed\>  
\</options\>

**Signification :** “virtual wall at domain edge \= true” active une barrière répulsive aux limites du domaine pour empêcher les cellules de sortir et assurer une stabilité de la simulation. “Disable automated spring adhesions=false” permet d’assurer l’intégrité du bloc cellulaire en appliquant les valeurs d'adhésion entre cellules. Enfin, “random seed=19” fixe la graine du générateur de nombres aléatoire. C’est ce qui assure la reproductibilité si couplé à un nombre de cœur qui vaut 1\. (Explication : voir Partie 2 : Allocation des ressources de calculs (parallel))

# **Partie 3 : Microenvironnement** {#partie-3-:-microenvironnement}

- ## **Facteur de croissance** {#facteur-de-croissance}

\<microenvironment_setup\>  
 \<variable name="growth_factor" units="dimensionless" ID="0"\>\#Nom de la variable, du facteur  
 \<physical_parameter_set\>  
 \<diffusion_coefficient units="micron^2/min"\>500.0\</diffusion_coefficient\> \#Vitesse de propagation de la molécule dans le milieu

                \<decay\_rate units="1/min"\>0.0\</decay\_rate\> \#dégradation de la molécule au cours du temps. (0 \= pas de dégradation)
            \</physical\_parameter\_set\>
            \<initial\_condition units="dimensionless"\>0.0\</initial\_condition\> \# concentration à t=0
            \<Dirichlet\_boundary\_condition units="dimensionless" enabled="True"\>0.0\</Dirichlet\_boundary\_condition\> \#Indique si les conditions aux limites changent. Avec True, ce n’est pas le cas donc la source est un puit infini
            \<Dirichlet\_options\>
                \<boundary\_value ID="xmin" enabled="False" /\>
                \<boundary\_value ID="xmax" enabled="False" /\>
                \<boundary\_value ID="ymin" enabled="True"\>0\</boundary\_value\> \#valeur de concentration à la limite inférieure
                \<boundary\_value ID="ymax" enabled="True"\>1\</boundary\_value\> \#valeur de concentration à la limite supérieure
                \<boundary\_value ID="zmin" enabled="False" /\>
                \<boundary\_value ID="zmax" enabled="False" /\>
            \</Dirichlet\_options\>
        \</variable\>

**Signification :** Ce bloc définit le facteur influençant la croissance et la différenciation des cellules. Les paramètres Dirichlet activés sur ymin imposent que les ressources viennent exclusivement du bas du domaine de simulation.

**Justification biologique :** Le growth factor agit comme un morphogène passif. Sa forte concentration au niveau basal maintient le "pool" de cellules souches prolifératives. Sa diminution progressive vers la surface (due à la consommation par les couches cellulaires, ou "effet puits") simule l'éloignement des vaisseaux sanguins. C'est cette baisse de concentration qui signale aux cellules supérieures de sortir du cycle cellulaire et d'entrer en différenciation terminale, reproduisant ainsi la stratification naturelle de l'épithélium.

**Source :** Goodhill, G. J. (1997). Diffusion in Axon Guidance. _European Journal Of Neuroscience_, _9_(7), 1414‑1421. https://doi.org/10.1111/j.1460-9568.1997.tb01496.x  
_Rheinwald, J. G., & Green, H. (1977). Epidermal growth factor and the multiplication of cultured human epidermal keratinocytes. Nature, 265(5593), 421‑424. https://doi.org/10.1038/265421a0_

**Citations** : _Goodhill, G. J. (1997)_ : Cet article établit les coefficients de diffusion à environ 10−7 cm2/s ce qui équivaut à 600 μm2/min.  
_Rheinwald, J. G., & Green, H. (1977)_ : la survie et la prolifération des kératinocytes dépendent de la disponibilité continue de facteurs de croissance. L'utilisation de conditions de Dirichlet simule mathématiquement un réservoir infini de nutriments.

- ## **Facteur cancéreux** {#facteur-cancéreux}

        \<variable name="cancer\_factor" units="dimensionless" ID="1"\>
            \<physical\_parameter\_set\>
                \<diffusion\_coefficient units="micron^2/min"\>500.0\</diffusion\_coefficient\>
                \<decay\_rate units="1/min"\>0.0\</decay\_rate\>
            \</physical\_parameter\_set\>
            \<initial\_condition units="dimensionless"\>0.0\</initial\_condition\>
            \<Dirichlet\_boundary\_condition units="dimensionless" enabled="False"\>0.0\</Dirichlet\_boundary\_condition\>
            \<Dirichlet\_options\>
                \<boundary\_value ID="xmin" enabled="False" /\>
                \<boundary\_value ID="xmax" enabled="False" /\>
                \<boundary\_value ID="ymin" enabled="False" /\>
                \<boundary\_value ID="ymax" enabled="False" /\>
                \<boundary\_value ID="zmin" enabled="False" /\>
                \<boundary\_value ID="zmax" enabled="False" /\>
            \</Dirichlet\_options\>
        \</variable\>


**Signification :** Cette variable représente le "cocktail" de molécules de signalisation (Cytokines, Chimiokines, Facteurs de croissance comme le TGF-β ou l'IL-6) sécrété par les cellules cancéreuses. Son rôle est de diffuser dans le tissu environnant pour : Recruter des cellules immunitaires (Chimiotactisme). Transformer les cellules saines en complices (Fibroblastes → CAF, Macrophages → CAM).

**Justification :** Le coefficient de diffusion (500.0 μm2/min) modélise la mobilité de macromolécules protéiques à travers la matrice extracellulaire dense du stroma. Le taux de dégradation fixé à 0.0 est une hypothèse permettant au signal tumoral de se propager sans perte sur de plus longues distances, assurant un recrutement efficace des cellules immunitaires à travers tout le domaine de simulation

**Sources :** Goodhill, G. J. (1997b). Diffusion in Axon Guidance. _European Journal Of Neuroscience_, _9_(7), 1414‑1421. https://doi.org/10.1111/j.1460-9568.1997.tb01496.x

- ## **Death factor** {#death-factor}

\<variable name="death_factor" units="dimensionless" ID="2"\>  
 \<physical_parameter_set\>  
 \<diffusion_coefficient units="micron^2/min"\>500.0\</diffusion_coefficient\>  
 \<decay_rate units="1/min"\>0.0\</decay_rate\>  
 \</physical_parameter_set\>  
 \<initial_condition units="dimensionless"\>0.0\</initial_condition\>  
 \<Dirichlet_boundary_condition units="dimensionless" enabled="True"\>0.0\</Dirichlet_boundary_condition\>  
 \<Dirichlet_options\>  
 \<boundary_value ID="xmin" enabled="False" /\>  
 \<boundary_value ID="xmax" enabled="False" /\>  
 \<boundary_value ID="ymin" enabled="True"\>1\</boundary_value\>  
 \<boundary_value ID="ymax" enabled="True"\>0\</boundary_value\>  
 \<boundary_value ID="zmin" enabled="False" /\>  
 \<boundary_value ID="zmax" enabled="False" /\>  
 \</Dirichlet_options\>  
 \</variable\>

**Signification :** Cette variable représente un substrat imposé par l’environnement sous la forme d’un gradient qui décroît vers la surface. Comme son nom l'indique, il permet d’éliminer des cellules en fonction de sa valeur en agissant comme un signal déclencheur d’apoptose.

**Justification :** Pour assurer la cohérence des couches du tissu et éviter des dérives et invasions de certains types de cellules, un gradient est mis en place. Ainsi, suivant les plages de valeurs définies, les cellules “ont le droit ou non” de se trouver à certains emplacements. Cela prévient notamment une invasion du tissu conjonctif par simple division ou dérive mécanique de l’épithélium basal en instaurant une zone létale dans les couches profondes.

- ## **Facteur de dégradation de la membrane basale (mmp_factor)** {#facteur-de-dégradation-de-la-membrane-basale-(mmp_factor)}

  \<variable name="mmp_factor" units="dimensionless" ID="3"\>

\<physical_parameter_set\>  
\<diffusion_coefficient units="micron^2/min"\>500.0\</diffusion_coefficient\>  
\<decay_rate units="1/min"\>10\</decay_rate\>  
\</physical_parameter_set\>  
\<initial_condition units="dimensionless"\>0.0\</initial_condition\>  
\<Dirichlet_boundary_condition units="dimensionless" enabled="False"\>0.0\</Dirichlet_boundary_condition\>  
\<Dirichlet_options\>  
\<boundary_value ID="xmin" enabled="False"/\>  
\<boundary_value ID="xmax" enabled="False"/\>  
\<boundary_value ID="ymin" enabled="False"/\>  
\<boundary_value ID="ymax" enabled="False"/\>  
\<boundary_value ID="zmin" enabled="False"/\>  
\<boundary_value ID="zmax" enabled="False"/\>  
\</Dirichlet_options\>  
\</variable\>

**Signification :** Cette variable représente la sécrétion de métalloprotéase matricielle qui représentent, dans cette simulation, le facteur dégradant la membrane basale. Bien que sécrété par de nombreux agents cancéreux, il est principalement relâché par les cellules cancéreuses mésenchymateuses.

**Justification biologique :** Avec une valeur de diffusion moyenne et un decay rate assez important ce facteur agit à faible portée pour assurer une dégradation de la membrane basale au contact des cellules cancéreuses mésenchymateuses.

**Sources :** [LELONGT, Brigitte, RONCO, Pierre et PIEDAGNEL, Rémi, 2002\. Métalloprotéases matricielles : infidélités à la matrice extracellulaire. _médecine/sciences_. mai 2002\. Vol. 18, n° 5, pp. 519‑521. DOI 10.1051/medsci/2002185519.](https://www.zotero.org/google-docs/?JITPxA)

**Citations :** “Chaque protéase cible ds substrats qui lui sont propres, si bien que la famille des MMP est capable de dégrader tous les constituants des membranes basales et des MEC”

- ## **Chimiotactisme induit par les CAF (CAF_chemotaxis)** {#chimiotactisme-induit-par-les-caf-(caf_chemotaxis)}

\<variable name="CAF_chemotaxis" units="dimensionless" ID="4"\>  
 \<physical_parameter_set\>  
 \<diffusion_coefficient units="micron^2/min"\>1000.0\</diffusion_coefficient\>  
 \<decay_rate units="1/min"\>0.00001\</decay_rate\>  
 \</physical_parameter_set\>  
 \<initial_condition units="dimensionless"\>0.0\</initial_condition\>  
 \<Dirichlet_boundary_condition units="dimensionless" enabled="False"\>0.0\</Dirichlet_boundary_condition\>  
 \<Dirichlet_options\>  
 \<boundary_value ID="xmin" enabled="False" /\>  
 \<boundary_value ID="xmax" enabled="False" /\>  
 \<boundary_value ID="ymin" enabled="False" /\>  
 \<boundary_value ID="ymax" enabled="False" /\>  
 \<boundary_value ID="zmin" enabled="False" /\>  
 \<boundary_value ID="zmax" enabled="False" /\>  
 \</Dirichlet_options\>  
 \</variable\>

**Signification :** Cette variable représente un signal de chimiotactisme à longue portée sécrété exclusivement par les fibroblastes activés (CAF). Elle attire activement les cellules cancéreuses du compartiment épithélial vers le tissu conjonctif.

**Justification biologique :** Contrairement aux facteurs paracrines à courte portée, ce facteur possède un coefficient de diffusion très élevé (1000.0) et un taux de dégradation quasi nul (0.00001). Cette configuration physique permet au signal de saturer l'ensemble du domaine et de créer un gradient stable sur de longues distances. Les cellules mésenchymateuses y sont extrêmement sensibles, ce qui fait de ce substrat le véritable moteur de la migration invasive.

**Sources :** Goodhill, G. J. (1997b). Diffusion in Axon Guidance. _European Journal Of Neuroscience_, _9_(7), 1414‑1421. https://doi.org/10.1111/j.1460-9568.1997.tb01496.x  
Kalluri, R. (2016). The biology and function of fibroblasts in cancer. Nature Reviews. Cancer, 16(9), 582‑598. https://doi.org/10.1038/nrc.2016.73  
Sahai, E., Astsaturov, I., Cukierman, E., DeNardo, D. G., Egeblad, M., Evans, R. M., Fearon, D., Greten, F. R., Hingorani, S. R., Hunter, T., Hynes, R. O., Jain, R. K., Janowitz, T., Jorgensen, C., Kimmelman, A. C., Kolonin, M. G., Maki, R. G., Powers, R. S., Puré, E.,. . . Werb, Z. (2020). A framework for advancing our understanding of cancer-associated fibroblasts. Nature Reviews. Cancer, 20(3), 174‑186. https://doi.org/10.1038/s41568-019-0238-1

**Citations :** _Kalluri, R. (2016)_ **:** Détaille comment les CAFs corrompus par la tumeur sécrètent des facteurs pour soutenir et diriger activement la migration des cellules malignes.  
_Sahai et al. (2020)_ : Explique que la communication chimique entre le stroma et la tumeur est essentielle pour créer un environnement permissif à l'invasion et aux métastases.

- ## **Autres** {#autres}

        \<options\>
            \<calculate\_gradients\>true\</calculate\_gradients\>
            \<track\_internalized\_substrates\_in\_each\_agent\>false\</track\_internalized\_substrates\_in\_each\_agent\>
        \</options\>

  \</microenvironment_setup\>

**Signification :** Si activée, cette option demanderait à chaque cellule de tenir un "livre de comptes" stockant la quantité totale de chaque produit chimique qu'elle a stocké depuis sa naissance.

# **Partie 4 : types cellulaires** {#partie-4-:-types-cellulaires}

- ## **Epithélium basal** {#epithélium-basal}

\<cell_definitions\>  
 \<cell_definition name="epi_basal" ID="0"\>  
 \<phenotype\>  
 \<cycle code="5" name="live"\>  
 \<phase_durations units="min"\>  
 \<duration index="0" fixed_duration="true"\>1440\</duration\> \#Cycle cellulaire, ici 1 division toutes les 1440 minutes soit 24 heures  
 \</phase_durations\>  
 \<standard_asymmetric_division enabled="False"\> \#pour définir si une cellule mère d’un type peut se diviser en cellule filles d’autres types ici ce n’est pas le cas  
 \<asymmetric_division_probability name="epi_basal" \#les probabilités de division asymetriques units="dimensionless"\>1.0\</asymmetric_division_probability\>  
 \<asymmetric_division_probability name="epi_inter" units="dimensionless"\>0\</asymmetric_division_probability\>  
 \<asymmetric_division_probability name="epi_sup" units="dimensionless"\>0\</asymmetric_division_probability\>  
 \<asymmetric_division_probability name="cancer" units="dimensionless"\>0\</asymmetric_division_probability\>  
 \<asymmetric_division_probability name="TCell" units="dimensionless"\>0\</asymmetric_division_probability\>  
 \<asymmetric_division_probability name="membrane" units="dimensionless"\>0\</asymmetric_division_probability\>  
 \<asymmetric_division_probability name="conjonctif" units="dimensionless"\>0\</asymmetric_division_probability\>  
 \<asymmetric_division_probability name="cancer_mes" units="dimensionless"\>0\</asymmetric_division_probability\>  
 \<asymmetric_division_probability name="CAF" units="dimensionless"\>0\</asymmetric_division_probability\>  
 \</standard_asymmetric_division\>  
 \</cycle\>  
 \<death\>  
 \<model code="100" name="apoptosis"\>  
 \<death_rate units="1/min"\>9.31667e-05\</death_rate\> \#probabilité par minute de mort  
 \<phase_transition_rates units="1/min"\>  
 \<rate start_index="0" end_index="1" fixed_duration="false"\>0.001938\</rate\>  
 \</phase_transition_rates\>  
 \<parameters\>  
 \<unlysed_fluid_change_rate units="1/min"\>0.05\</unlysed_fluid_change_rate\>  
 \<lysed_fluid_change_rate units="1/min"\>0\</lysed_fluid_change_rate\>  
 \<cytoplasmic_biomass_change_rate units="1/min"\>1.66667e-02\</cytoplasmic_biomass_change_rate\>  
 \<nuclear_biomass_change_rate units="1/min"\>5.83333e-03\</nuclear_biomass_change_rate\>  
 \<calcification_rate units="1/min"\>0\</calcification_rate\>  
 \<relative_rupture_volume units="dimensionless"\>2.0\</relative_rupture_volume\>  
 \</parameters\>  
 \</model\>  
 \<model code="101" name="necrosis"\>  
 \<death_rate units="1/min"\>0.0\</death_rate\> \#probabilité de mort par nécrose  
 \<phase_transition_rates units="1/min"\>  
 \<rate start_index="0" end_index="1" fixed_duration="false"\>9000000000.0\</rate\>  
 \<rate start_index="1" end_index="2" fixed_duration="true"\>1.15741e-05\</rate\>  
 \</phase_transition_rates\>  
 \<parameters\>  
 \<unlysed_fluid_change_rate units="1/min"\>1.11667e-02\</unlysed_fluid_change_rate\>  
 \<lysed_fluid_change_rate units="1/min"\>8.33333e-4\</lysed_fluid_change_rate\>  
 \<cytoplasmic_biomass_change_rate units="1/min"\>5.33333e-05\</cytoplasmic_biomass_change_rate\>  
 \<nuclear_biomass_change_rate units="1/min"\>2.16667e-4\</nuclear_biomass_change_rate\>  
 \<calcification_rate units="1/min"\>7e-05\</calcification_rate\>  
 \<relative_rupture_volume units="dimensionless"\>2.0\</relative_rupture_volume\>  
 \</parameters\>  
 \</model\>  
 \</death\>  
 \<volume\>  
 \<total units="micron^3"\>4913\</total\> \#volume de la cellule en micron cube  
 \<fluid_fraction units="dimensionless"\>0.75\</fluid_fraction\>  
 \<nuclear units="micron^3"\>540\</nuclear\>  
 \<fluid_change_rate units="1/min"\>0.05\</fluid_change_rate\>  
 \<cytoplasmic_biomass_change_rate units="1/min"\>0.0045\</cytoplasmic_biomass_change_rate\>  
 \<nuclear_biomass_change_rate units="1/min"\>0.0055\</nuclear_biomass_change_rate\>  
 \<calcified_fraction units="dimensionless"\>0.0\</calcified_fraction\>  
 \<calcification_rate units="1/min"\>0.0\</calcification_rate\>  
 \<relative_rupture_volume units="dimensionless"\>2\</relative_rupture_volume\>  
 \</volume\>  
 \<mechanics\>  
 \<cell_cell_adhesion_strength units="micron/min"\>0.4\</cell_cell_adhesion_strength\> \#force d’attraction coolant les cellules entre elles  
 \<cell_cell_repulsion_strength units="micron/min"\>100.0\</cell_cell_repulsion_strength\> \#force empêchant les cellules de se chevaucher  
 \<relative_maximum_adhesion_distance units="dimensionless"\>1.25\</relative_maximum_adhesion_distance\> \#distance limite pour maintenir une adhésion entre 2 cellules  
 \<cell_adhesion_affinities\> \#force d’attraction entre ce type cellulaire (basal) et les autres  
 \<cell_adhesion_affinity name="epi_basal"\>1.0\</cell_adhesion_affinity\>  
 \<cell_adhesion_affinity name="epi_inter"\>1.0\</cell_adhesion_affinity\>  
 \<cell_adhesion_affinity name="epi_sup"\>1.0\</cell_adhesion_affinity\>  
 \<cell_adhesion_affinity name="cancer"\>1.0\</cell_adhesion_affinity\>  
 \<cell_adhesion_affinity name="TCell"\>1.0\</cell_adhesion_affinity\>  
 \<cell_adhesion_affinity name="membrane"\>1.0\</cell_adhesion_affinity\>  
 \<cell_adhesion_affinity name="conjonctif"\>0.5\</cell_adhesion_affinity\>  
 \<cell_adhesion_affinity name="cancer_mes"\>1.0\</cell_adhesion_affinity\>  
 \<cell_adhesion_affinity name="CAF"\>1.0\</cell_adhesion_affinity\>  
 \</cell_adhesion_affinities\>  
 \<options\>  
 \<set_relative_equilibrium_distance enabled="false" units="dimensionless"\>1.8\</set_relative_equilibrium_distance\>  
 \<set_absolute_equilibrium_distance enabled="false" units="micron"\>15.12\</set_absolute_equilibrium_distance\>  
 \</options\>  
 \<attachment_elastic_constant units="1/min"\>0.01\</attachment_elastic_constant\> \#raideur des liens entre cellules  
 \<attachment_rate units="1/min"\>10.0\</attachment_rate\> \#vitesse de formation des liens entre cellules  
 \<detachment_rate units="1/min"\>0.0\</detachment_rate\> \#vitesse de rupture  
 \<maximum_number_of_attachments\>12\</maximum_number_of_attachments\> \#nombre de connexions mécaniques simultanées  
 \</mechanics\>  
 \<motility\>  
 \<speed units="micron/min"\>1.0\</speed\> \#vitesse de base de déplacement de la cellule  
 \<persistence_time units="min"\>5.0\</persistence_time\> \#durée pendant laquelle une cellule garde la même direction  
 \<migration_bias units="dimensionless"\>0.5\</migration_bias\> \#influence d’un signal externe sur la direction de migration  
 \<options\>  
 \<enabled\>false\</enabled\> \#interrupteur de mobilité on=true off=false  
 \<use_2D\>true\</use_2D\>  
 \<chemotaxis\>  
 \<enabled\>false\</enabled\> \#interrupteur de chimiotactisme  
 \<substrate\>growth_factor\</substrate\> \#substrat pour lequel s’applique (ou non) le chimiotactisme  
 \<direction\>1\</direction\>  
 \</chemotaxis\>  
 \<advanced_chemotaxis\>  
 \<enabled\>false\</enabled\> \#de même que précédemment  
 \<normalize_each_gradient\>false\</normalize_each_gradient\>  
 \<chemotactic_sensitivities\>  
 \<chemotactic_sensitivity substrate="growth_factor"\>0.0\</chemotactic_sensitivity\> \#sensibilité au type de substrat par chimiotactisme  
 \<chemotactic_sensitivity substrate="cancer_factor"\>0.0\</chemotactic_sensitivity\>  
 \<chemotactic_sensitivity substrate="CAF_factor"\>0.0\</chemotactic_sensitivity\>  
 \<chemotactic_sensitivity substrate="death_factor"\>0.0\</chemotactic_sensitivity\>  
 \<chemotactic_sensitivity substrate="mmp_factor"\>0.0\</chemotactic_sensitivity\>  
 \<chemotactic_sensitivity substrate="CAF_chemotaxis"\>0.0\</chemotactic_sensitivity\>  
 \</chemotactic_sensitivities\>  
 \</advanced_chemotaxis\>  
 \</options\>  
 \</motility\>  
 \<secretion\> \#sécrétion de facteur par le type cellulaire  
 \<substrate name="growth_factor"\>  
 \<secretion_rate units="1/min"\>0.0\</secretion_rate\> \#quantité de substance relâché par le type cellulaire  
 \<secretion_target units="substrate density"\>1.0\</secretion_target\> \#quantité visée par le type cellulaire dans son voisinage au max  
 \<uptake_rate units="1/min"\>0.0\</uptake_rate\> \#vitesse d’absorption des substances du milieu  
 \<net_export_rate units="total substrate/min"\>0.0\</net_export_rate\>  
 \</substrate\>  
 \<substrate name="cancer_factor"\>  
 \<secretion_rate units="1/min"\>0.0\</secretion_rate\>  
 \<secretion_target units="substrate density"\>1.0\</secretion_target\>  
 \<uptake_rate units="1/min"\>0.0\</uptake_rate\>  
 \<net_export_rate units="total substrate/min"\>0.0\</net_export_rate\>  
 \</substrate\>  
 \<substrate name="CAF_factor"\>  
 \<secretion_rate units="1/min"\>0.0\</secretion_rate\>  
 \<secretion_target units="substrate density"\>1.0\</secretion_target\>  
 \<uptake_rate units="1/min"\>0.0\</uptake_rate\>  
 \<net_export_rate units="total substrate/min"\>0.0\</net_export_rate\>  
 \</substrate\>  
 \<substrate name="death_factor"\>  
 \<secretion_rate units="1/min"\>0.0\</secretion_rate\>  
 \<secretion_target units="substrate density"\>1.0\</secretion_target\>  
 \<uptake_rate units="1/min"\>0.0\</uptake_rate\>  
 \<net_export_rate units="total substrate/min"\>0.0\</net_export_rate\>  
 \</substrate\>  
 \<substrate name="mmp_factor"\>  
 \<secretion_rate units="1/min"\>0.0\</secretion_rate\>  
 \<secretion_target units="substrate density"\>1.0\</secretion_target\>  
 \<uptake_rate units="1/min"\>0.0\</uptake_rate\>  
 \<net_export_rate units="total substrate/min"\>0.0\</net_export_rate\>  
 \</substrate\>  
 \<substrate name="CAF_chemotaxis"\>  
 \<secretion_rate units="1/min"\>0.0\</secretion_rate\>  
 \<secretion_target units="substrate density"\>1.0\</secretion_target\>  
 \<uptake_rate units="1/min"\>0.0\</uptake_rate\>  
 \<net_export_rate units="total substrate/min"\>0.0\</net_export_rate\>  
 \</substrate\>  
 \</secretion\>  
 \<cell_interactions\>  
 \<apoptotic_phagocytosis_rate units="1/min"\>0.0\</apoptotic_phagocytosis_rate\>  
 \<necrotic_phagocytosis_rate units="1/min"\>0.0\</necrotic_phagocytosis_rate\>  
 \<other_dead_phagocytosis_rate units="1/min"\>0.0\</other_dead_phagocytosis_rate\> \#phagocytose du type cellulaire sur les autres types  
 \<live_phagocytosis_rates\>  
 \<phagocytosis_rate name="epi_basal" units="1/min"\>0.0\</phagocytosis_rate\>  
 \<phagocytosis_rate name="epi_inter" units="1/min"\>0.0\</phagocytosis_rate\>  
 \<phagocytosis_rate name="epi_sup" units="1/min"\>0.0\</phagocytosis_rate\>  
 \<phagocytosis_rate name="cancer" units="1/min"\>0.0\</phagocytosis_rate\>  
 \<phagocytosis_rate name="TCell" units="1/min"\>0.0\</phagocytosis_rate\>  
 \<phagocytosis_rate name="membrane" units="1/min"\>0.0\</phagocytosis_rate\>  
 \<phagocytosis_rate name="conjonctif" units="1/min"\>0.0\</phagocytosis_rate\>  
 \<phagocytosis_rate name="cancer_mes" units="1/min"\>0.0\</phagocytosis_rate\>  
 \<phagocytosis_rate name="CAF" units="1/min"\>0.0\</phagocytosis_rate\>  
 \</live_phagocytosis_rates\>  
 \<attack_rates\> \#taux d’attaque du type cellulaire sur les autres  
 \<attack_rate name="epi_basal" units="1/min"\>0.0\</attack_rate\>  
 \<attack_rate name="epi_inter" units="1/min"\>0.0\</attack_rate\>  
 \<attack_rate name="epi_sup" units="1/min"\>0.0\</attack_rate\>  
 \<attack_rate name="cancer" units="1/min"\>0.0\</attack_rate\>  
 \<attack_rate name="TCell" units="1/min"\>0.0\</attack_rate\>  
 \<attack_rate name="membrane" units="1/min"\>0.0\</attack_rate\>  
 \<attack_rate name="conjonctif" units="1/min"\>0.0\</attack_rate\>  
 \<attack_rate name="cancer_mes" units="1/min"\>0.0\</attack_rate\>  
 \<attack_rate name="CAF" units="1/min"\>0.0\</attack_rate\>  
 \</attack_rates\>  
 \<attack_damage_rate units="1/min"\>1.0\</attack_damage_rate\>  
 \<attack_duration units="min"\>0.1\</attack_duration\>  
 \<fusion_rates\>  
 \<fusion_rate name="epi_basal" units="1/min"\>0.0\</fusion_rate\>  
 \<fusion_rate name="epi_inter" units="1/min"\>0.0\</fusion_rate\>  
 \<fusion_rate name="epi_sup" units="1/min"\>0.0\</fusion_rate\>  
 \<fusion_rate name="cancer" units="1/min"\>0.0\</fusion_rate\>  
 \<fusion_rate name="TCell" units="1/min"\>0.0\</fusion_rate\>  
 \<fusion_rate name="membrane" units="1/min"\>0.0\</fusion_rate\>  
 \<fusion_rate name="conjonctif" units="1/min"\>0.0\</fusion_rate\>  
 \<fusion_rate name="cancer_mes" units="1/min"\>0.0\</fusion_rate\>  
 \<fusion_rate name="CAF" units="1/min"\>0.0\</fusion_rate\>  
 \</fusion_rates\>  
 \</cell_interactions\>  
 \<cell_transformations\>  
 \<transformation_rates\> \#transformation du type cellulaire en un autre  
 \<transformation_rate name="epi_basal" units="1/min"\>0.0\</transformation_rate\>  
 \<transformation_rate name="epi_inter" units="1/min"\>0.0\</transformation_rate\>  
 \<transformation_rate name="epi_sup" units="1/min"\>0.0\</transformation_rate\>  
 \<transformation_rate name="cancer" units="1/min"\>0.0\</transformation_rate\>  
 \<transformation_rate name="TCell" units="1/min"\>0.0\</transformation_rate\>  
 \<transformation_rate name="membrane" units="1/min"\>0.0\</transformation_rate\>  
 \<transformation_rate name="conjonctif" units="1/min"\>0.0\</transformation_rate\>  
 \<transformation_rate name="cancer_mes" units="1/min"\>0\</transformation_rate\>  
 \<transformation_rate name="CAF" units="1/min"\>0.0\</transformation_rate\>  
 \</transformation_rates\>  
 \</cell_transformations\>  
 \<cell_integrity\>  
 \<damage_rate units="1/min"\>0.0\</damage_rate\>  
 \<damage_repair_rate units="1/min"\>0.0\</damage_repair_rate\>  
 \</cell_integrity\>  
 \</phenotype\>  
 \<custom_data\>  
 \<sample conserved="false" units="dimensionless" description=""\>0.0\</sample\>  
 \</custom_data\>  
 \<initial_parameter_distributions enabled="false"\>  
 \</initial_parameter_distributions\>  
 \</cell_definition\>

**Signification** : Le type epi_basal est l'agent moteur de la simulation, conçu pour agir comme une couche de cellules souches prolifératives et structurelles. Avec un cycle de division de 24 heures et une motilité désactivée, il génère la biomasse du tissu tout en restant physiquement ancré à la base du domaine.

**Justification :** L'absence de sécrétion et de motilité active, combinée à une forte répulsion mécanique (100.0), permet de créer une barrière physique stable qui pousse les cellules filles vers le haut par simple pression de division. Ce paramétrage garantit que l'architecture épithéliale ne s'effondre pas et que la progression des cellules vers les couches supérieures est un phénomène purement convectif, dépendant de la prolifération basale.

**Sources:** Squier, C. A., & Kremer, M. J. (2001). Biology of Oral Mucosa and Esophagus. _JNCI Monographs_, _2001_(29), 7‑15. https://doi.org/10.1093/oxfordjournals.jncimonographs.a003443

**Citations** : Description du comportement de “tapis-roulant” du tissu épithélial

- ## **Epithelium intermédiaire** {#epithelium-intermédiaire}

\<cell_definition name="epi_inter" ID="1"\>  
 \<cycle code="5" name="live"\>  
 \<rate start_index="0" end_index="0" fixed_duration="false"\>0\</rate\>  
 \</cycle\>  
 \<volume\>  
 \<total units="micron^3"\>4913\</total\>  
 \</volume\>  
 \<mechanics\>  
 \<cell_cell_adhesion_strength units="micron/min"\>0.4\</cell_cell_adhesion_strength\>  
 \<cell_cell_repulsion_strength units="micron/min"\>100.0\</cell_cell_repulsion_strength\>  
 \<attachment_rate units="1/min"\>10.0\</attachment_rate\>  
 \<detachment_rate units="1/min"\>0.0\</detachment_rate\>  
 \</mechanics\>  
 \<motility\>  
 \<enabled\>false\</enabled\>  
 \</motility\>  
 \<secretion\>  
 \<substrate name="growth_factor"\>  
 \<secretion_rate units="1/min"\>0.0\</secretion_rate\>  
 \</substrate\>  
 \</secretion\>  
\</cell_definition\>

**Signification :** Ce bloc définit les propriétés biologiques des cellules du _Stratum Spinosum_ ou couche épineuse, qui constituent la structure principale de l'épithélium. Le paramètre de cycle cellulaire réglé sur un taux de transition nul impose un arrêt total de la division, rendant ces cellules quiescentes. Une probabilité infime de 5e-7 min⁻¹ est introduite pour permettre une transformation spontanée en cellule cancéreuse. Le volume est fixé à 4913 µm³, identique à la couche basale pour la stabilité du calcul, et l'adhésion est maintenue à 0.4 pour assurer la cohésion intercellulaire.

**Justification :** La quiescence est justifiée par le fait que les kératinocytes du _stratum spinosum_ sont biologiquement post-mitotiques ; dès qu'ils quittent la niche basale, ils sortent définitivement du cycle cellulaire pour se consacrer à la synthèse de kératines, évitant ainsi une hyperplasie tissulaire pathologique. Le taux de transformation modélise le risque stochastique de "dé-différenciation" ou d'initiation tumorale suite à l'accumulation de dommages à l'ADN, indépendamment des signaux du stroma. L'adhésion modélise la présence dense de desmosomes, ou jonctions fortes, qui donnent son nom à la couche épineuse et assurent la fonction de barrière mécanique du tissu. Enfin, l'immobilité est maintenue car l'ascension de ces cellules vers la surface est un processus passif dû à la pression de prolifération sous-jacente.

**Sources :** Squier, C. A., & Kremer, M. J. (2001). Biology of Oral Mucosa and Esophagus. _JNCI Monographs_, _2001_(29), 7‑15. https://doi.org/10.1093/oxfordjournals.jncimonographs.a003443_._

_Hanahan, D., & Weinberg, R. A. (2011). Hallmarks of Cancer : The Next Generation. Cell, 144(5), 646‑674. https://doi.org/10.1016/j.cell.2011.02.013_

**Citations :** _Squier, C. A. & Kremer, M. J. (2001) : L’épithélium intermédiaire constitue la couche de résistance mécanique majeure de l'épithélium buccal._

_Hanahan, D., & Weinberg, R. A. (2011) : L'instabilité génomique et la mutation spontanée sont les initiateurs fondamentaux du processus tumoral._

- ## **Epithelium supérieur** {#epithelium-supérieur}

\<cell_definition name="epi_sup" ID="2"\>  
 \<cycle code="5" name="live"\>  
 \<phase_transition_rates units="1/min"\>  
 \<rate start_index="0" end_index="0" fixed_duration="false"\>0\</rate\>  
 \</phase_transition_rates\>  
 \</cycle\>  
 \<volume\>  
 \<total units="micron^3"\>729\</total\>  
 \</volume\>  
 \<mechanics\>  
 \<cell_cell_adhesion_strength units="micron/min"\>0.4\</cell_cell_adhesion_strength\>  
 \<cell_cell_repulsion_strength units="micron/min"\>100.0\</cell_cell_repulsion_strength\>  
 \<attachment_rate units="1/min"\>10.0\</attachment_rate\>  
 \<detachment_rate units="1/min"\>0.0\</detachment_rate\>  
 \</mechanics\>  
 \<secretion\>  
 \<substrate name="growth_factor"\>  
 \<secretion_rate units="1/min"\>0.01\</secretion_rate\>  
 \</substrate\>  
 \</secretion\>  
\</cell_definition\>

**Signification :** Ce bloc définit les propriétés des cellules de la couche superficielle (Stratum Superficiale), stade ultime de la vie du kératinocyte avant la mort. Le cycle cellulaire est totalement bloqué avec un taux de transition nul. Le volume cellulaire cible est réduit drastiquement à 729 µm³, contre 4913 µm³ pour les couches inférieures.

**Justification :** La réduction du volume à 729 est le paramètre clé de ce bloc, car elle permet de modéliser physiquement la compaction et la déshydratation des cellules avant leur desquamation. Techniquement, l'arrêt de la division est maintenu par un taux de transition à 0, confirmant l'engagement irréversible dans la différenciation terminale. Contrairement aux autres types sains, l'epi_sup possède un taux de sécrétion de 0.01 pour le growth_factor, agissant comme une source de signal paracrine qui influence le comportement des couches basales situées en dessous. La conservation d'une répulsion maximale (100.0) et d'un attachement fort (10.0) garantit que cette couche superficielle reste une barrière protectrice cohérente malgré sa faible épaisseur.

- ## **Cellules cancéreuses** {#cellules-cancéreuses}

\<cell_definition name="cancer" ID="3"\> \<cycle code="5" name="live"\> \<rate start_index="0" end_index="0" fixed_duration="true"\>0.000002\</rate\> \</cycle\> \<death\> \<death_rate units="1/min"\>0.31667e-05\</death_rate\> \</death\> \<mechanics\> \<cell_cell_adhesion_strength units="micron/min"\>0.0\</cell_cell_adhesion_strength\> \<cell_cell_repulsion_strength units="micron/min"\>10.0\</cell_cell_repulsion_strength\> \</mechanics\> \<motility\> \<speed units="micron/min"\>0.1\</speed\> \<enabled\>true\</enabled\> \<chemotactic_sensitivity substrate="CAF_chemotaxis"\>0.5\</chemotactic_sensitivity\> \</motility\> \<secretion\> \<substrate name="cancer_factor"\> \<secretion_rate units="1/min"\>2.0\</secretion_rate\> \</substrate\> \</secretion\> \<cell_transformations\> \<transformation_rate name="cancer_mes" units="1/min"\>0.00003\</transformation_rate\> \</cell_transformations\> \</cell_definition\>

**Signification :** Le type cancer représente les cellules épithéliales malignes ayant acquis une autonomie de croissance et une capacité d'évasion des signaux de différenciation. Contrairement aux kératinocytes sains, cet agent maintient un cycle de division actif avec un taux de transition de 0.000002 par minute et présente une résistance accrue à l'apoptose avec un taux de mortalité réduit à 0.31667e-05. Il constitue l'initiateur de la pathologie dans le modèle en activant sa motilité propre et en commençant à sécréter des signaux oncogéniques pour manipuler le stroma environnant.

**Justification :** La suppression totale de l'adhésion intercellulaire (0.0) associée à une faible répulsion mécanique (10.0) permet de simuler la perte de cohésion tissulaire nécessaire à l'invasion. L'activation de la motilité avec une vitesse de 0.1 micron par minute et une sensibilité chimiotactique de 0.5 envers le signal des fibroblastes modélise l'attraction exercée par le microenvironnement tumoral sur les cellules malignes. Le taux de sécrétion élevé de cancer_factor (2.0) assure le recrutement paracrine des cellules stromales, tandis que le taux de transformation de 0.00003 vers le type cancer_mes prépare la transition épithélio-mésenchymateuse indispensable à l'invasion profonde du domaine de simulation.

**Sources :** Hanahan, D., & Weinberg, R. A. (2011). Hallmarks of Cancer : The Next Generation. _Cell_, _144_(5), 646‑674. https://doi.org/10.1016/j.cell.2011.02.013

Friedl, P., & Alexander, S. (2011). Cancer Invasion and the Microenvironment : Plasticity and Reciprocity. Cell, 147(5), 992‑1009. https://doi.org/10.1016/j.cell.2011.11.016

Kalluri, R., & Weinberg, R. A. (2009). The basics of epithelial-mesenchymal transition. Journal Of Clinical Investigation, 119(6), 1420‑1428. https://doi.org/10.1172/jci39104

**Citations** : _Hanahan, D. & Weinberg, R. A. (2011)_ : L'autonomie des signaux de croissance est sans doute le trait le plus fondamental des cellules cancéreuses, impliquant la capacité de croître sans stimulation externe.

_Friedl, P. & Alexander, S. (2011)_ : La perte des jonctions cellule-cellule médiées par la E-cadhérine est une condition préalable au détachement des cellules cancéreuses individuelles de la masse tumorale primaire.

_Kalluri, R. & Weinberg, R. A. (2009)_ : Pendant la transition épithélio-mésenchymateuse (EMT), les cellules épithéliales réduisent leurs structures d'adhésion et réorganisent leur cytosquelette pour devenir isolées et mobiles.

- ## **T Cell** {#t-cell}

\<cell_definition name="TCell" ID="4"\>  
 \<cycle code="5" name="live"\>  
 \<phase_transition_rates units="1/min"\>  
 \<rate start_index="0" end_index="0" fixed_duration="true"\>0.0\</rate\>  
 \</phase_transition_rates\>  
 \</cycle\>  
 \<death\>  
 \<model code="100" name="apoptosis"\>  
 \<death_rate units="1/min"\>0.0\</death_rate\>  
 \</model\>  
 \</death\>  
 \<mechanics\>  
 \<cell_cell_adhesion_strength units="micron/min"\>0.0\</cell_cell_adhesion_strength\>  
 \<cell_cell_repulsion_strength units="micron/min"\>10.0\</cell_cell_repulsion_strength\>  
 \</mechanics\>  
 \<motility\>  
 \<speed units="micron/min"\>0.5\</speed\>  
 \<persistence_time units="min"\>5.0\</persistence_time\>  
 \<migration_bias units="dimensionless"\>0.5\</migration_bias\>  
 \<enabled\>true\</enabled\>  
 \<chemotaxis\>  
 \<enabled\>true\</enabled\>  
 \<substrate\>cancer_factor\</substrate\>  
 \<direction\>1\</direction\>  
 \</chemotaxis\>  
 \</motility\>  
\</cell_definition\>

**Signification :** Le type TCell modélise les lymphocytes T infiltrant la tumeur pour assurer la surveillance immunitaire du tissu. Dans la simulation, cet agent est une unité patrouilleuse hautement mobile qui ne se divise pas et ne meurt pas spontanément. Sa fonction principale est de localiser les cellules malignes en suivant les signaux biochimiques émis par ces dernières, agissant ainsi comme un capteur dynamique de la pathologie

**Justification :** La réactivité immunitaire repose sur une motilité active réglée à 0.5 micron par minute associée à un chimiotactisme positif vers le substrat cancer_factor. Cette sensibilité permet aux cellules T de remonter le gradient chimique pour converger précisément vers les foyers tumoraux. L'absence totale d'adhésion et la faible répulsion garantissent que ces cellules peuvent circuler librement entre les kératinocytes et à travers la membrane sans être entravées par les structures physiques du tissu.

**Sources :**

Friedl, P., & Weigelin, B. (2008). Interstitial leukocyte migration and immune function. Nature Immunology, 9(9), 960‑969. https://doi.org/10.1038/ni.f.212

**Citations :**

_Friedl, P. & Weigelin, B. (2008)_ : La migration des leucocytes dans les tissus denses est optimisée par une faible adhésion et une grande plasticité de mouvement dirigé par les gradients chimiques.

- ## **Membrane** {#membrane}

\<cell_definition name="membrane" ID="5"\>  
 \<cycle code="5" name="live"\>  
 \<phase_transition_rates units="1/min"\>  
 \<rate start_index="0" end_index="0" fixed_duration="true"\>0.0\</rate\>  
 \</phase_transition_rates\>  
 \</cycle\>  
 \<death\>  
 \<model code="100" name="apoptosis"\>  
 \<death_rate units="1/min"\>0.0\</death_rate\>  
 \</model\>  
 \</death\>  
 \<mechanics\>  
 \<cell_cell_adhesion_strength units="micron/min"\>0.0\</cell_cell_adhesion_strength\>  
 \<cell_cell_repulsion_strength units="micron/min"\>100.0\</cell_cell_repulsion_strength\>  
 \<attachment_rate units="1/min"\>0.0\</attachment_rate\>  
 \<detachment_rate units="1/min"\>0.0\</detachment_rate\>  
 \</mechanics\>  
 \<motility\>  
 \<enabled\>false\</enabled\>  
 \</motility\>  
\</cell_definition\>

**Signification :** Le type membrane modélise la lame basale, la structure extracellulaire fine et dense qui sépare l'épithélium du tissu conjonctif (stroma) sous-jacent. Dans la simulation, cet agent n'est pas une cellule vivante au sens biologique mais un composant structurel statique. Il est défini par une absence totale de cycle cellulaire et de mortalité, agissant comme une fondation inerte et permanente pour l'architecture du tissu buccal.

**Justification :** Le rôle de barrière physique est traduit techniquement par une force de répulsion maximale de 100.0, ce qui empêche les cellules épithéliales de s'enfoncer dans le stroma en conditions normales. L'absence d'adhésion (0.0) et de motilité (enabled \= false) garantit que la membrane reste parfaitement immobile et ne "colle" pas aux cellules, permettant aux kératinocytes de glisser le long de cette paroi lors de la poussée proliférative. Cette configuration est essentielle pour maintenir la compartimentation du modèle et tester l'invasion tumorale (franchissement de cette barrière).

**Sources :** Squier, C. A., & Kremer, M. J. (2001). Biology of Oral Mucosa and Esophagus. JNCI Monographs, 2001(29), 7‑15. https://doi.org/10.1093/oxfordjournals.jncimonographs.a003443

Yurchenco, P. D. (2010). Basement Membranes : Cell Scaffoldings and Signaling Platforms. Cold Spring Harbor Perspectives In Biology, 3(2), a004911. https://doi.org/10.1101/cshperspect.a004911

**Citations** : _Squier, C. A. & Kremer, M. J. (2001)_ : La membrane basale est une structure complexe qui sépare l'épithélium du tissu conjonctif et sert de support à la couche de cellules basales.

_Yurchenco, P. D. (2011)_ : Les membranes basales sont des matrices extracellulaires denses qui fournissent un soutien structurel aux épithéliums et agissent comme des barrières sélectives au mouvement des cellules.

- ## **Fibroblasts** {#fibroblasts}

\<cell_definition name="fibroblast" ID="6"\>

    \<cycle code="5" name="live"\>

        \<phase\_transition\_rates units="1/min"\>

            \<rate start\_index="0" end\_index="0" fixed\_duration="true"\>0.0\</rate\>

        \</phase\_transition\_rates\>

    \</cycle\>

    \<volume\>

        \<total units="micron^3"\>4913\</total\>

    \</volume\>

    \<mechanics\>

        \<cell\_cell\_adhesion\_strength units="micron/min"\>0.4\</cell\_cell\_adhesion\_strength\>

        \<cell\_cell\_repulsion\_strength units="micron/min"\>100.0\</cell\_cell\_repulsion\_strength\>

        \<attachment\_rate units="1/min"\>0.0\</attachment\_rate\>

        \<detachment\_rate units="1/min"\>0.0\</detachment\_rate\>

    \</mechanics\>

    \<motility\>

        \<speed units="micron/min"\>0.1\</speed\>

        \<enabled\>false\</enabled\>

    \</motility\>

    \<cell\_transformations\>

        \<transformation\_rates\>

            \<transformation\_rate name="CAF" units="1/min"\>0.0\</transformation\_rate\>

        \</transformation\_rates\>

    \</cell\_transformations\>

\</cell_definition\>

**Signification** : Le type fibroblast représente les cellules principales du stroma ou tissu conjonctif situé sous la membrane basale. Dans la simulation, il s'agit d'une cellule structurelle stable qui occupe l'espace dermique avec un volume de 4913 microns cubes, similaire aux cellules de la couche basale. Son rôle est de maintenir l'architecture du tissu de soutien et de servir de réservoir de cellules capables de se transformer en fibroblastes associés au cancer (CAF) dès que le microenvironnement devient tumoral.

**Justification** : Le fibroblaste est configuré dans un état de quiescence absolue avec un taux de division de 0.0 et une motilité désactivée, reflétant son état stationnaire dans un tissu sain. Les paramètres mécaniques de répulsion maximale (100.0) et d'adhésion (0.4) assurent que le stroma reste dense et résistant à la compression exercée par l'épithélium sus-jacent. Bien que le taux de transformation initial vers le type CAF soit nul dans les paramètres de base, cette ligne de code est cruciale car elle permet au modèle d'activer la réaction stromale et la sécrétion de chimiokines dès que la tumeur commence à diffuser ses signaux.

**Sources** : Squier, C. A., & Kremer, M. J. (2001). Biology of Oral Mucosa and Esophagus. JNCI Monographs, 2001(29), 7‑15. https://doi.org/10.1093/oxfordjournals.jncimonographs.a003443

Kalluri, R. (2016). The biology and function of fibroblasts in cancer. Nature Reviews. Cancer, 16(9), 582‑598. https://doi.org/10.1038/nrc.2016.73

Hanahan, D., & Weinberg, R. A. (2011). Hallmarks of Cancer : The Next Generation. Cell, 144(5), 646‑674. https://doi.org/10.1016/j.cell.2011.02.013

**Citations** : _Squier, C. A. & Kremer, M. J. (2001)_ : Les fibroblastes sont les cellules prédominantes du tissu conjonctif, responsables de la production de la matrice extracellulaire qui soutient l'épithélium.

_Kalluri, R. (2016)_ : Les fibroblastes résidents possèdent une plasticité phénotypique leur permettant de répondre aux signaux de stress tissulaire et de se transformer en myofibroblastes ou en fibroblastes associés au cancer.

_Hanahan, D. & Weinberg, R. A. (2011)_ : Le stroma n'est pas simplement un échafaudage structurel mais un partenaire dynamique qui co-évolue avec les cellules cancéreuses durant la progression maligne.

- ## **Cancer mésenchymateux** {#cancer-mésenchymateux}

  \<cell_definition name="cancer_mes" ID="7"\>

      \<cycle code="5" name="live"\>

          \<phase\_transition\_rates units="1/min"\>

              \<rate start\_index="0" end\_index="0" fixed\_duration="true"\>0.000002\</rate\>

          \</phase\_transition\_rates\>

      \</cycle\>

      \<death\>

          \<model code="100" name="apoptosis"\>

              \<death\_rate units="1/min"\>0.31667e-05\</death\_rate\>

          \</model\>

      \</death\>

      \<mechanics\>

          \<cell\_cell\_adhesion\_strength units="micron/min"\>0.0\</cell\_cell\_adhesion\_strength\>

          \<cell\_cell\_repulsion\_strength units="micron/min"\>10.0\</cell\_cell\_repulsion\_strength\>

      \</mechanics\>

      \<motility\>

          \<speed units="micron/min"\>0.7\</speed\>

          \<persistence\_time units="min"\>5.0\</persistence\_time\>

          \<migration\_bias units="dimensionless"\>0.7\</migration\_bias\>

          \<enabled\>true\</enabled\>

          \<chemotaxis\>

              \<enabled\>true\</enabled\>

              \<substrate\>CAF\_chemotaxis\</substrate\>

              \<direction\>1\</direction\>

          \</chemotaxis\>

      \</motility\>

      \<secretion\>

          \<substrate name="mmp\_factor"\>

              \<secretion\_rate units="1/min"\>1.0\</secretion\_rate\>

          \</substrate\>

      \</secretion\>

  \</cell_definition\>

**Signification** : Le type cancer_mes représente les cellules cancéreuses ayant achevé leur transition épithélio-mésenchymateuse (EMT). Cet agent n'est plus une cellule ancrée au tissu mais une entité hautement invasive capable de migrer à travers les structures denses du stroma. Contrairement au type cancer initial, il possède un phénotype migratoire exacerbé et la capacité de dégrader physiquement les obstacles environnants, tels que la membrane basale ou la matrice extracellulaire, pour faciliter l'extension de la tumeur.

**Justification** : Le passage à un état invasif est traduit par une augmentation massive de la vitesse de déplacement, réglée à 0.7 micron par minute, soit sept fois la vitesse du type cancer épithélial. Le modèle renforce également le biais de migration à 0.7, indiquant un mouvement beaucoup plus directionnel en réponse aux gradients chimiques. L'innovation majeure de ce type réside dans la sécrétion du mmp_factor (métalloprotéinases de la matrice) avec un taux de 1.0, permettant la destruction locale des barrières physiques comme la membrane basale. Cette combinaison de motilité rapide et de sécrétion protéolytique modélise avec précision l'invasion métastatique.

**Sources** : Kalluri, R., & Weinberg, R. A. (2009). The basics of epithelial-mesenchymal transition. Journal Of Clinical Investigation, 119(6), 1420‑1428. https://doi.org/10.1172/jci39104

Friedl, P., & Alexander, S. (2011). Cancer Invasion and the Microenvironment : Plasticity and Reciprocity. Cell, 147(5), 992‑1009. https://doi.org/10.1016/j.cell.2011.11.016

Hanahan, D., & Weinberg, R. A. (2011). Hallmarks of Cancer : The Next Generation. Cell, 144(5), 646‑674. https://doi.org/10.1016/j.cell.2011.02.013

**Citations** : _Kalluri, R. & Weinberg, R. A. (2009)_ : La transition épithélio-mésenchymateuse confère aux cellules cancéreuses des traits typiques des cellules souches et une motilité accrue, essentielles pour l'invasion systémique.

_Friedl, P. & Alexander, S. (2011)_ : L'invasion tumorale efficace nécessite une coordination entre la force de propulsion cellulaire et la dégradation de la matrice extracellulaire par des enzymes comme les MMP.

_Hanahan, D. & Weinberg, R. A. (2011)_ : L'activation de l'invasion et des métastases est une caractéristique majeure du cancer, impliquant des changements dans l'adhésion cellulaire et l'acquisition d'un phénotype migratoire.

- ## **CAF** {#caf}

\<cell_definition name="CAF" ID="8"\>  
 \<cycle code="5" name="live"\>  
 \<phase_transition_rates units="1/min"\>  
 \<rate start_index="0" end_index="0" fixed_duration="true"\>0.0\</rate\>  
 \</phase_transition_rates\>  
 \</cycle\>  
 \<mechanics\>  
 \<cell_cell_adhesion_strength units="micron/min"\>0.4\</cell_cell_adhesion_strength\>  
 \<cell_cell_repulsion_strength units="micron/min"\>100.0\</cell_cell_repulsion_strength\>  
 \</mechanics\>  
 \<motility\>  
 \<speed units="micron/min"\>0.1\</speed\>  
 \<enabled\>false\</enabled\>  
 \</motility\>  
 \<secretion\>  
 \<substrate name="CAF_chemotaxis"\>  
 \<secretion_rate units="1/min"\>1.0\</secretion_rate\>  
 \<secretion_target units="substrate density"\>1.0\</secretion_target\>  
 \</substrate\>  
 \</secretion\>  
\</cell_definition\>

**Signification :** Le type CAF représente les fibroblastes associés au cancer, des cellules du stroma qui ont été détournées par la tumeur pour favoriser sa progression. Dans la simulation, ces cellules sont issues de la transformation des fibroblastes sains sous l'influence des signaux tumoraux. Contrairement aux fibroblastes normaux, ils deviennent des centres de signalisation active en sécrétant le substrat CAF_chemotaxis, qui agit comme un phare chimique pour guider les cellules cancéreuses invasives à travers le tissu conjonctif.

**Justification :** Le rôle de facilitateur de l'invasion est traduit techniquement par l'activation de la sécrétion du substrat CAF_chemotaxis avec un taux de 1.0. Ce paramètre est crucial car il crée le gradient chimique que les cellules cancer_mes utilisent pour s'orienter. Bien que leur motilité soit désactivée (enabled \= false) pour simuler leur position fixe dans la matrice extracellulaire, leur influence sur la dynamique globale est majeure. La conservation d'une répulsion mécanique maximale (100.0) garantit que ces cellules continuent de former la charpente physique du stroma tout en modifiant radicalement sa composition chimique.

**Sources :** Kalluri, R. (2016). The biology and function of fibroblasts in cancer. Nature Reviews. Cancer, 16(9), 582‑598. https://doi.org/10.1038/nrc.2016.73

Sahai, E., Astsaturov, I., Cukierman, E., DeNardo, D. G., Egeblad, M., Evans, R. M., Fearon, D., Greten, F. R., Hingorani, S. R., Hunter, T., Hynes, R. O., Jain, R. K., Janowitz, T., Jorgensen, C., Kimmelman, A. C., Kolonin, M. G., Maki, R. G., Powers, R. S., Puré, E.,. . . Werb, Z. (2020). A framework for advancing our understanding of cancer-associated fibroblasts. Nature Reviews. Cancer, 20(3), 174‑186. https://doi.org/10.1038/s41568-019-0238-1

Hanahan, D., & Weinberg, R. A. (2011). Hallmarks of Cancer : The Next Generation. Cell, 144(5), 646‑674. https://doi.org/10.1016/j.cell.2011.02.013

**Citations** : _Kalluri, R. (2016)_ : Les fibroblastes associés au cancer constituent une population hétérogène de cellules qui favorisent activement la croissance tumorale, l'angiogenèse et l'invasion par le remodelage de la matrice et la sécrétion de facteurs de croissance.

_Sahai, E. et al. (2020)_ : La capacité des CAF à sécréter des chimiokines crée des pistes chimiotactiques qui facilitent le recrutement et la migration directionnelle des cellules cancéreuses à travers le stroma.

_Hanahan, D. & Weinberg, R. A. (2011)_ : Le recrutement et l'activation de fibroblastes de soutien en tant que CAF sont des étapes essentielles pour créer un microenvironnement permissif à l'expansion maligne

- ## **Conjonctif** {#conjonctif}

\<cell_definition name="conjonctif" ID="9"\>  
 \<cycle code="5" name="live"\>  
 \<phase_transition_rates units="1/min"\>  
 \<rate start_index="0" end_index="0" fixed_duration="true"\>0.0\</rate\>  
 \</phase_transition_rates\>  
 \</cycle\>  
 \<volume\>  
 \<total units="micron^3"\>4913\</total\>  
 \</volume\>  
 \<mechanics\>  
 \<cell_cell_adhesion_strength units="micron/min"\>0.4\</cell_cell_adhesion_strength\>  
 \<cell_cell_repulsion_strength units="micron/min"\>100.0\</cell_cell_repulsion_strength\>  
 \<attachment_rate units="1/min"\>0.0\</attachment_rate\>  
 \<detachment_rate units="1/min"\>0.0\</detachment_rate\>  
 \</mechanics\>  
 \<motility\>  
 \<enabled\>false\</enabled\>  
 \</motility\>  
\</cell_definition\>

**Signification :** Le type conjonctif représente la matrice structurelle de la lamina propria, le tissu de soutien situé sous l'épithélium buccal. Dans la simulation, ces agents ne sont pas des cellules individuelles actives mais servent à remplir le volume du stroma, créant une densité tissulaire réaliste autour des vaisseaux et des fibres de collagène. Ils assurent la stabilité volumétrique du compartiment conjonctif et servent de support physique sur lequel les fibroblastes et les cellules immunitaires sont positionnés.

**Justification :** Le tissu conjonctif est défini comme un milieu statique avec un taux de division nul et une motilité désactivée, garantissant que le volume dermique reste constant durant toute la simulation. La force de répulsion maximale de 100.0 est cruciale car elle permet au stroma de résister à la pression exercée par la croissance de l'épithélium et de maintenir la séparation nette imposée par la membrane basale. L'adhésion de 0.4 assure une cohésion interne suffisante pour que le stroma ne se désagrège pas lors des interactions mécaniques avec les cellules cancéreuses invasives.

**Sources :** Squier, C. A., & Kremer, M. J. (2001). Biology of Oral Mucosa and Esophagus. JNCI Monographs, 2001(29), 7‑15. https://doi.org/10.1093/oxfordjournals.jncimonographs.a003443

Hanahan, D., & Weinberg, R. A. (2011). Hallmarks of Cancer : The Next Generation. Cell, 144(5), 646‑674. https://doi.org/10.1016/j.cell.2011.02.013

**Citations** _: Squier, C. A. & Kremer, M. J. (2001)_ : Le tissu conjonctif de la muqueuse buccale, ou lamina propria, est un système de soutien mécanique et métabolique essentiel pour l'épithélium sus-jacent.

_Hanahan, D. & Weinberg, R. A. (2011)_ : Le stroma environnant joue un rôle actif dans la progression tumorale, agissant initialement comme une barrière physique avant d'être remodelé par les cellules malignes.

# **Partie 5 : Initialisation et Paramètres** {#partie-5-:-initialisation-et-paramètres}

Ci-dessous la forme est l’explication des fichiers csv utilisés par PhysiCell. Il est possible de les modifier directement mais il est conseillé d’utiliser PhysiCell Studio pour éviter toutes erreurs de syntaxe par exemple.

- ## **Conditions initiales (cells.csv)** {#conditions-initiales-(cells.csv)}

  X | Y | Z | Type

  _Figure 1 : screenshot d’une partie du fichier cells.csv représentant quelques cellules._

Le fichier se présente sous la forme d’un tableau à 4 colonnes (x, y, z, type). **X** représente l’axe horizontal, **Y** l’axe vertical avec donc la surface de l’épithélium qui tend vers des valeurs de 300 et les couches profondes vers \-300. **Z** est fixé à 0 pour toutes les cellules car la simulation est en 2D. Enfin le **type** correspond au type cellulaire concerné.

- ## **Les règles (cells_rules.csv)** {#les-règles-(cells_rules.csv)}

  _Figure 2 : screenshot d’une partie du fichier cells_rules.csv représentant quelques règles._

  De la même manière, il s’agit d’un tableau mais cette fois-ci composé de 8 colonnes. De gauche à droite la signification de chaque colonne est la suivante :

- **Cell type** : le type cellulaire sur lequel est appliquée la règle. Tous les agents de ce type sont concernés. Epi_basal \= toutes les cellules de l’épithélium basal.

- **Signal** : Une substance chimique, un état… qui déclenche la réaction. Growth_factor \= le facteur de croissance

- **Direction** : indique si le comportement réagit à une augmentation (increases) ou une diminution (decreases) du signal.

- **Behavior** : L’action ou transformation qui est appliquée sur le type cellulaire. Transform to epi_inter \= transformation des cellules basales en cellules épithélium intermédiaire.

Ainsi, ces 4 premières colonnes définissent qualitativement la règle. Pour la première ligne de la figure 2 : Les cellules épithéliales basales se transforment en épithélium intermédiaire lorsque le facteur de croissance augmente dans leur voisinage.  
Les trois colonnes suivantes donnent des indications quantitatives se basant sur la fonction de Hill :

- **Saturation value** : la valeur limite vers laquelle tend le comportement. Ici 1 pour la première règle, signifie que quand la courbe a atteint son maximum, la cellule a 100% de chance de se transformer.

- **Half-max** : c’est le seuil 50/50. Toujours pour le même exemple, la valeur vaut 0.3. Ainsi quand le signal atteint 0.3 la cellule à 50% de chance de se transformer.

- **Hill power** : c’est la “pente” de la courbe de Hill. Pour de petites valeurs, la courbe (dé)croit lentement et la transition des cellules est moins déterministe. Pour des valeurs élevées comme ici 400, ça rend le phénomène comparable à un interrupteur on/off. Ce concept est illustré en figure 3 ci dessous.

  _Figure 3 : screenshot d’une partie du fichier cells.csv représentant quelques cellules._

- **Apply to dead** : La dernière colonne est binaire 1 (vrai) ou 0 (faux) pour déterminer si la règle s’applique aux cellules mortes. Pour toutes les règles utilisées la valeur fausse est appliquée.
