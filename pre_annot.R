# Requires to set the working directory to parent folder of this with setwd() (outside pre_annot folder)



#Matrice de comptage du nombre de gènes en communs
get.annot.matrix <- function(SeuratObj, diff.expr.genes){
  # Pre-attribution of clusters
  nClusters <- length(levels(SeuratObj)) # Get number of clusters
  
  diffGenesRef <- read.csv("./pre_annot/diffGenesBase.csv", sep = ";")
  cellTypes <- colnames(diffGenesRef)
  
  
  ## Init matrix of association
  type.annot.matrix <- matrix(0, nrow=length(cellTypes), ncol=nClusters)
  colnames(type.annot.matrix) <- 0:(nClusters-1)
  rownames(type.annot.matrix) <- cellTypes
  
  for (i in 1:ncol(diffGenesRef)){ # for each cell type
    genes.list <- c()
    n.significant.genes <- 0 # Pour la moyenne, compter le nombre de gènes statistiquement sig
    irow <- cellTypes[i]
    for(gene in diffGenesRef[,i]){ # for each marker gene
      if(!(gene %in% genes.list)){ # Don't do the same gene twice
        whiches <- which(diff.expr.genes$gene == gene) # find if in diff expressed genes
        for(j in whiches){ # for each row
          if(diff.expr.genes$p_val_adj[j]<0.05 && diff.expr.genes$avg_log2FC[j] > 0){ # use adjusted p-value
            jcol <- diff.expr.genes$cluster[j]
            type.annot.matrix[irow, jcol] <- type.annot.matrix[irow, jcol] + 1
          }
        }
        append(gene, genes.list)
      }
    }
  }
  
  return(type.annot.matrix)
}




get.avg.matrix <- function(SeuratObj, diff.expr.genes){
  # Pre-attribution of clusters
  nClusters <- length(levels(SeuratObj)) # Get number of clusters
  
  diffGenesRef <- read.csv("./pre_annot/diffGenesBase.csv", sep = ";")
  cellTypes <- colnames(diffGenesRef)
  
  
  ## Init matrix of association
  type.avg.matrix <- matrix(0, nrow=length(cellTypes), ncol=nClusters)
  colnames(type.avg.matrix) <- 0:(nClusters-1)
  rownames(type.avg.matrix) <- cellTypes
  
  for (i in 1:ncol(diffGenesRef)){ # for each cell type
    genes.list <- c()
    n.significant.genes <- 0 # Pour la moyenne, compter le nombre de gènes statistiquement sig
    irow <- cellTypes[i]
    for(gene in diffGenesRef[,i]){ # for each marker gene
      if(!(gene %in% genes.list)){ # Don't do the same gene twice
        whiches <- which(diff.expr.genes$gene == gene) # find if in diff expressed genes
        for(j in whiches){ # for each row
          if(diff.expr.genes$p_val_adj[j]<0.05 && diff.expr.genes$avg_log2FC[j] > 0){ # use adjusted p-value + we don't use genes that are underexpressed
            jcol <- diff.expr.genes$cluster[j]
            type.avg.matrix[irow, jcol] <- type.avg.matrix[irow, jcol] + diff.expr.genes$avg_log2FC[j]
            n.significant.genes <- n.significant.genes + 1
          }
        }
        append(gene, genes.list)
      }
    }
    if(n.significant.genes){
      type.avg.matrix[cellTypes[i], ] <- type.avg.matrix[cellTypes[i],] / n.significant.genes # Diviser par le nombre de gènes trouvés pour faire une moyenne
    }
  }
  
  
  # Make the associations more clear by putting max on the diag
  # ⚠️ 1 cluster != 1 cell type : 
  # orderedPredictions <- type.avg.matrix[,HungarianSolver(-type.avg.matrix)$pairs[,2]] # Requires RcppHungarian Package
  
  
  return(type.avg.matrix)
}

mark_knowns <- function(diff.expr.genes){ # Add a column to identify known or useless genes 
  # IE genes that were already analysed and/or that didn't help identify previoulsy
  # Seeks to accelerate the identification of difficult clusters
  diffGenesRef <- read.csv("./pre_annot/diffGenesBase.csv", sep = ";")
  uselessGenes <- read.csv("./pre_annot/uselessGenes.csv", sep = ";")
  diff.expr.genes$known <- FALSE
  
  for(col in diffGenesRef){
    for(gene in col){
      whiches <- which(diff.expr.genes$gene == gene) # find if in diff expressed genes
        for(j in whiches){ # for each row
          diff.expr.genes$known[j] <- TRUE
        }
    }
  }
  
  for(gene in uselessGenes){
    whiches <- which(diff.expr.genes$gene == gene) # find if in diff expressed genes
    for(j in whiches){ # for each row
      diff.expr.genes$known[j] <- TRUE
    }
  }
  
  return(diff.expr.genes)
}


display_heatmap <- function(matrix){ # requires to put - r-reshape and - r-ggplot2 in yml file
  # Old way heatmap(type.avg.matrix, Rowv = TRUE, Colv=NA)
  
  library("ggplot2")
  library("reshape")
  data_melt <- melt(matrix)
  ggp <- ggplot(data_melt, aes(X2, X1)) +                           # Create heatmap with ggplot2
    geom_tile(aes(fill = value)) + scale_x_continuous(label = floor) +
    xlab("Cluster") +
    ylab("Cell type")
  
  return(ggp) 
}

pre_labels <- function(type.avg.matrix, seuil = 3){
  nClusters <- ncol(type.avg.matrix)
  # Assign names to clusters and accuraccy
  ## Accuracy
  is.accurate <- function(score.vector){
    score <-0
    score.vector <- sort(score.vector, decreasing = TRUE)
    if(score.vector[2] > 0){
      score <- score.vector[1]/score.vector[2] # score = highest / second highest
    }
    return(score >= seuil)
  }
  
  clusters.annot <- c()
  
  clusters.annot.plain <- c() # With clusters with same names
  for(i in 1:nClusters){ # For each cluster
    itype <- names(which.max(type.avg.matrix[,i])) # get the predicted type name
    itype.modified <- itype
    
    if(itype %in% clusters.annot.plain){
      # add a number to differentiate similar clusters
      itype.modified <- paste(itype, sum(clusters.annot.plain == itype)+1) 
    }
    clusters.annot.plain <- append(clusters.annot.plain, itype)

    if(!is.accurate(type.avg.matrix[,i])){
      itype.modified <- paste("⚠️", itype.modified)
    }
    clusters.annot <- append(clusters.annot, itype.modified)
  }
  
  
  return(clusters.annot)
}
