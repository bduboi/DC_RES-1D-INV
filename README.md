# DC_RES-1D-INV
Ce dépos a été créé pour les étudiants du stage de terrain en géophysique. Il apporte un code et des explications pour l'inversion de mesures de resistivité 1D avec la configuration Schlumberger.
## Problème inverse
Le principe du problème inverse est de simuler les données obtenues à partir d'un **modèle**, à l'aide d'un **opérateur direct**. Cet opérateur peut être une solution analytique aux équations physiques qui régissent l'expérience, ou une résolution numérique de l'équation. L'opérateur en question simule une grandeur physique (potentiel / tension , champs...) en fonction des **paramètres** de notre domaine, à savoir le sous sol. Ces **paramètres** sont dans notre cas la resistivité / conductivité et l'épaisseur des couches, mais on pourrait très bien faire une autre expérience pour déterminer d'autres paramètres physique.
Dans le cas de la resistivité 1D, le modèle que nous utilisons est une solution analytique pour un modèle de Terre constituée de couches horizontales au dessus d'un **demi espace infini**. Le code a été écrit par Mathias Scheunert (TUBAF). Il prend en entrée M resistivités et M-1 épaisseurs (puisque la dernière couche est infinie), et la liste des N espacements d'électrodes. Avec ces trois informations, il donne en sortie une liste de N resistivités apparentes.
La procédure d'inversion utilise ce modèle pour calculer une liste de resistivité apparente et la comparer avec celle qu'on a vraiment mesuré. Elle se base sur une solution des moindres carrés et à chaque itération calcule une direction vers laquelle les **paramètres** (resistivité, épaisseur) doivent être modifiés pour se rapprocher le mieux que possible du jeu de donnée. 
Il est aussi possible de prendre en compte les erreurs associées aux données en les multipliant par l'inverse de leur erreur associée, ce qui leur donne un poids : les points ayant les plus grosses erreurs seront moins impactants que les autres.
## Liste des fichiers et explication
```
dcfwdf.m           : fonction matlab calculant la réponse à un modèle de resistivité / épaisseurs 1D
Sens_log.m         : fonction matlab utilisée à chaque itération pour calculer la **sensibilité** (ou Jacobienne) et trouver une direction de modification du modèle
Inversion_DC1Dres.m: fonction contenant la boucle d'optimisation du modèle, faisant appel aux deux premières fonctions.
Driver_inversion.m : Code matlab à partir duquel on choisit un modèle de départ, rentre les données et lançons l'inversion.
```
