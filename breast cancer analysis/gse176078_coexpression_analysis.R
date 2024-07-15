source("C:/Users/andyb/Downloads/FAST_CoExpress_functions.R")


# read in pre-processed data
gene_pairs <- readRDS("/path/to/save/gene_pairs")
processed_counts <- readRDS("/path/to/save/processed_counts")

# establish which gene will play the role of modulator gene
modgene <- "COL6A3"

# remove modgene from y1, y2 candidates
gene_pairs <- gene_pairs[(gene_pairs[,1]!=modgene) & (gene_pairs[,2]!=modgene),]

# set up covariate-dependence of parameters
eq1 <- y1 ~ x1 + x2 + x3
eq2 <- y2 ~ x1 + x2 + x3
eq3 <- ~ x2
eq4 <- ~ x2
eq5 <- ~ x1 + x2 + x3
eqlist <- list(eq1, eq2, eq3, eq4, eq5)

# run analyses
for (i in 1:nrow(gene_pairs)){
  
  gene1 <- gene_pairs[i,1]
  gene2 <- gene_pairs[i,2]
  
  ourdat <- processed_counts[, c(gene1, gene2, modgene, "cell_type")]
  colnames(ourdat) <- c("y1", "y2", "x1", "x2")
  
  ourdat$x2 <- 1*(ourdat$x2 == "CAFs myCAF-like")
  ourdat$x3 <- ourdat$x1*ourdat$x3
  
  ## fit FAST-CoExpress ##
  out_obj <- FAST.CoExpress.nb(formula = eqlist,
                               data = ourdat,
                               copula = "Gaussian")
  
  saveRDS(out_obj$output.mat, paste("/path/to/coefficient_storage/", gene1, "_", gene2, ".rds", sep=""))
  
}

