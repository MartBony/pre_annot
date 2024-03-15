⚠️ Fichier liste_donnees.csv contient la liste des études que nous analysons. La colonne "Done" indique si et qui a analysé le dataset.

# Mini tuto repo
1) Installer GIT (déjà présent sur BIRD)
2) Ouvrir le terminal et se placer dans le répertoire où se trouve le fichier de code avec `cd`
3) Cloner le dossier dans le répertoire avec
```
git clone https://github.com/MartBony/pre_annot.git
cd pre_annot
```

## Modifier et synchroniser les gènes marqueurs
**Modification**
Les gènes marqueurs sont dans le fichier diffGenesBase.csv et sont directement modifiables dans excel.
Les gènes qui ne nous informent pas sur la nature de la cellule sont à mettre dans uselessGenes.csv.

**Synchronisation**
Pour syncroniser le repo, il faut valider les changements en local :
```Bash
git add .
git commit -m "Modification des gènes marqueurs"
```
Puis importer la derniere version depuis le repo :
```Bash
git pull
```
Des problèmes de fusion de fichier peuvent survenir, dans ce cas il faut ouvrir les fichiers problématiques dans excel et les corriger soi-même.
Puis exporter sa version locale vers github :
```Bash
git push origin main
```

# Tuto du code
## Importer les fonctions

Clôner dans un dossier /pre_annot/ dans le répertoire du fichier de code.
⚠️ Ne pas oublier de définir le répertoire de travail au niveau du code :
```R
setwd("répertoire/absolu/de/mon/code")
```

Importer le code de pré-annotation
```R
source("../pre_annot/pre_annot.R", chdir=TRUE)
```

## Prérequis

Pour pré-annoter il faut :
- L'objet seurat de travail `SeurOBJ`
- Le tableau des gènes différentiellement exprimés `diff.expressed.genes` (obtenu après filtration des gènes marqueurs). Il contient les champs :
	- p_val : p-valeur classique.
	- avg_log2FC : Coefficient qui mesure la proportion d'expression différentielle.
	- p_val_adj : p-valeur plus stricte, utilisée pour la pré-annotation.
	- cluster
	- gene

## Fonctions de pré-annotation

### Matrice de comptage
On peut calculer la matrice de comptage des gènes en communs. Cette matrice correspond plus ou moins à la méthode manuelle d'annotation. Elle est biaisée pour les types cellulaires qui possèdent beaucoup de gènes dans le tableau mais elle offre quand même les meilleurs résultats. Il faut quand même jeter un oeil à la matrice d'expression différentielle.

```R
type.annot.matrix <- get_annot_matrix(SeurOBJ, diff.expressed.genes)
```

*PS : Un type cellulaire avec 1 seul gène (ex : LT cd8 pour l'instant) sera 0 ou 1 et apparaîtera prioritéaire sur tous les autres types. Ce cas est à éviter $\to$ ajouter des gènes.*

### Matrice d'expression différentielle
La matrice d'expression différentielle assigne à chaque couple (cluster, type cellulaire) un score qui correspond à la moyenne des coefficients avg_log2FC (log fold-change) des gènes en communs du couple. Voir (Wikipedia : Fold Change)[https://en.wikipedia.org/wiki/Fold_change]

Elle capture à quel point les gènes en commun sont différentiellement exprimés. Elle peut être utilse pour moduler les matrices précédentes.

```R
type.avg.matrix <- get_avg_matrix(SeurOBJ, diff.expressed.genes)
```

Plus le score est grand, plus les coefficients associés aux gènes en commun sont grands.
⚠️ Les résultats sont à croiser avec la matrice de comptage des gènes ci dessus.


### Matrice de comparaison
La matrice de comparaison représente pour chaque cluster et pour chaque type le pourcentage des génes de références du type qui sont différentiellement exprimés dans le cluster. C'est la matrice de comptage corrigée par le nombre de gènes de références qui varie selon type cellulaire. 

Elle est très biaisée pour les types cellulaires qui ont peu de gènes. Par exemple, si un type cellulaire a un seul gène marqueur, il arrivera toujours en premier dans les clusters où ce gène est différentiellement exprimé.

```R
type.corresp.matrix <- get_corresp_matrix(SeurOBJ, diff.expressed.genes)
```

### Visualiser les matrices
```R
display_heatmap(type.annot.matrix)
display_heatmap(type.annot.matrix, "Nombre de gènes") # Possibilité de mettre un titre
```

On peut comparer les deux matrices l'une à côté de l'autre :
```R
display_heatmap(type.annot.matrix) + display_heatmap(type.avg.matrix)
```

### Opérations sur les matrices
On peut croiser les informations de deux matrices en les multipliant élément par élément.
Je pense que celle qui a le plus de potentiel est :
```r
mtx <- type.annot.matrix * type.avg.matrix
```

Elle correspond au % de correspondance aux gènes corrigé par le niveau moyen d'expression différentielle de ces gènes.
Le problème c'est que je ne suis pas sur que c'est ultra rigoureux à utiliser pour la suite.


### Pré-assigner automatiquement

Choisit le type cellulaire le plus probable à partir de la matrice choisie. Utilise une fonction basique pour marquer les choix incertains.
```R
clusters.annot <- pre_labels(type.annot.matrix, seuil = 3) # Seuil optionnel
```

Plus le seuil est haut, plus les assignations serons considérées comme incertaines.

On peut afficher directement les types cellulaires choisis :
```R
clusters.annot
```

### Ajouter les labels à l'objet Seurat et à l'UMAP
```R
names(clusters.annot) <- levels(SeurOBJ)
LabeledSeurOBJ <- RenameIdents(SeurOBJ, clusters.annot)
DimPlot(LabeledSeurOBJ, reduction = "umap", label = TRUE, pt.size = 0.25) + NoLegend()
```


## Fonctions d'aide à l'annotation finale

### Maquer les gènes connus

On peut rajouter une colonne à la table des gènes différentiellement exprimés selon deux critères : 
- Est-ce que le gène est déjà dans les données de pré-assignation ?
- Est-ce que le gène est dans la liste des gènes inexploitables ? (fichier ./pre_annot/uselessGenes.csv)

```R
diff.expressed.genes <- mark_knowns(diff.expressed.genes)
```

Puis on peut continuer l'assignation manuelle en se concentrant sur les gènes pas encore exploités et en s'aidant de la pré-assignation.


## Autres
### Sauvegarder un plot quelconque en PNG
```R
save.plot.png(ma.heatmap, "./mon_nom_de_fichier.png")
```


## Exemple de code final
```R
# Générer les gènes marqueurs avant la pré-annotation
SeurOBJ.markers <- FindAllMarkers(SeurOBJ, only.pos = TRUE)
# Table of most unique genes per cluster
diff.expressed.genes <- SeurOBJ.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) # requires library(dplyr)


# Pré-annotation


source("../pre_annot/pre_annot.R", chdir=TRUE)

diff.expressed.genes <- mark_knowns(diff.expressed.genes) # optionnal 

type.annot.matrix <- get_annot_matrix(SeurOBJ, diff.expressed.genes)
type.avg.matrix <- get_avg_matrix(SeurOBJ, diff.expressed.genes)
type.corresp.matrix <- get_corresp_matrix(SeurOBJ, diff.expressed.genes)
type.modulated.matrix <- type.corresp.matrix * type.avg.matrix
display_heatmap(type.annot.matrix, "Nombre de gènes")
display_heatmap(type.avg.matrix, "Expression différentielle")
display_heatmap(type.corresp.matrix, "% de correspondance")
display_heatmap(type.modulated.matrix, "% de correspondance modulé par l'expression différentielle")

save.plot.png(display_heatmap(type.annot.matrix, "Nombre de gènes"), "./gene_count_matrix.png")



# Display annotations on UMAP
clusters.annot <- pre_labels(type.annot.matrix, seuil = 2)
clusters.annot


# Display annotations on UMAP
clusters.annot.modulated <- pre_labels(type.modulated.matrix, seuil = 2)
clusters.annot.modulated

names(clusters.annot) <- levels(SeurOBJ)
SeurOBJ.labeled <- RenameIdents(SeurOBJ, clusters.annot)
DimPlot(SeurOBJ.labeled, reduction = "umap", label = TRUE, pt.size = 0.25) + 
NoLegend() # Pas sur bird : + labs(title = "Annotation par comparaison")


names(clusters.annot.modulated) <- levels(SeurOBJ)
SeurOBJ.labeled.modulated <- RenameIdents(SeurOBJ, clusters.annot.modulated)
DimPlot(SeurOBJ.labeled.modulated, reduction = "umap", label = TRUE, pt.size = 0.25) + 
  NoLegend() # Pas sur bird : + labs(title = "Annotation par comparaison et expression")

```
