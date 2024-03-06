# Mini tuto repo
1) Installer GIT
2) Ouvrir le terminal et se placer dans le répertoire où se trouve le fichier de code avec `cd`
3) Cloner le dossier dans le répertoire avec `git clone https://github.com/MartBony/pre_annot.git`

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
source("./pre_annot/pre_annot.R")
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

On peut calculer la matrice de comptage des gènes en communs. Cette matrice correspond à la méthode manuelle d'annotation.
```R
type.annot.matrix <- get.annot.matrix(SeurOBJ, diff.expressed.genes)
```

### Matrice d'expression différentielle
La matrice de pré-annotation assigne à chaque couple (cluster, type cellulaire) un score qui correspond à la moyenne des coefficients avg_log2FC des gènes en communs du couple.

Elle capture à quel point les gènes en commun sont différentiellement exprimés. (pertinent ? par la p-valeur est validée de toute facon)

```R
type.avg.matrix <- get.avg.matrix(SeurOBJ, diff.expressed.genes)
```

Plus le score est grand, plus les coefficients associés aux gènes en commun sont grands.
⚠️ Les résultats sont à croiser avec la matrice de comptage des gènes ci dessus.


### Afficher la matrice de pré-annotation
```R
display_heatmap(my.matrix)
```

On peut comparer les deux matrices :
```R
display_heatmap(type.avg.matrix) + display_heatmap(type.alt.matrix)
```


### Pré-assigner automatiquement

Choisit le type cellulaire le plus probable à partir de la moyenne des coefficients avg_log2FC. Utilise une fonction basique pour marquer les choix incertains.
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






## Exemple de code final
```R

source("./pre_annot/pre_annot.R")

diff.expressed.genes <- mark_knowns(diff.expressed.genes)# optionnal 

type.annot.matrix <- get.annot.matrix(SeurOBJ, diff.expressed.genes)
type.avg.matrix <- get.avg.matrix(SeurOBJ, diff.expressed.genes)
display_heatmap(type.annot.matrix) + display_heatmap(type.avg.matrix)


# Display annotations on UMAP
clusters.annot <- pre_labels(type.annot.matrix, seuil = 6)
clusters.annot

names(clusters.annot) <- levels(SeurOBJ)
LabeledSeurOBJ.alt <- RenameIdents(SeurOBJ, clusters.annot)
DimPlot(LabeledSeurOBJ.alt, reduction = "umap", label = TRUE, pt.size = 0.25) + NoLegend()

```
