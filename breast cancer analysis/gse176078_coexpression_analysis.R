source("/path/to/FAST_CoExpress_functions.R")


# read in pre-processed data
gene_pairs <- readRDS("/path/to/genepairs_4471")
processed_counts <- readRDS("/path/to/processed_4471")


paramnames <- c("beta01", "beta11", "beta21", "beta31", "beta02", "beta12", "beta22", "beta32", "alpha01", "alpha11", "alpha02", "alpha12", "tau0", "tau1", "tau2", "tau3", "kappa01", "kappa11", "kappa02", "kappa12")


# set up covariate-dependence of parameters
eq1 <- y1 ~ x1 + x2 + x3
eq2 <- y2 ~ x1 + x2 + x3
eq3 <- ~ x2
eq4 <- ~ x2
eq5 <- ~ x1 + x2 + x3
eq6 <- ~ 1
eq7 <- ~ 1
eqlist <- list(eq1, eq2, eq3, eq4, eq5, eq6, eq7)


# run co-expression analysis on all gene pairs
modgene <- "COL6A3"

gene_pairs <- gene_pairs[(gene_pairs[,1]!=modgene) & (gene_pairs[,2]!=modgene),]


pvalvec <- rep(0, nrow(gene_pairs))
names(pvalvec) <- apply(gene_pairs, 1, function(a) paste(a[1],"_",a[2], sep=""))
for (i in 1:nrow(gene_pairs)){
  
  gene1 <- gene_pairs[i,1]
  gene2 <- gene_pairs[i,2]
  
  print(paste(gene1, gene2))
  
  ourdat <- processed_counts[, c(gene1, gene2, modgene, "cell_type")]
  colnames(ourdat) <- c("y1", "y2", "x1", "x2")
  
  ourdat$x2 <- 1*(ourdat$x2 == "CAFs myCAF-like")
  ourdat$x3 <- ourdat$x1*ourdat$x2
  
  ## fit FAST-CoExpress ##
  out_obj <- FAST.CoExpress.nb(formula = eqlist,
                               data = ourdat,
                               copula = "Gaussian")

  rownames(out_obj$output.mat) <- dimnames(out_obj$varcovmat)[[1]] <- dimnames(out_obj$varcovmat)[[2]] <- paramnames
  
  pvalvec[paste(gene1,"_",gene2,sep="")] <- out_obj$output.mat["tau3", "p.values"]
  
  saveRDS(out_obj, paste("/path/to/realdata results/", gene1, "_", gene2, ".rds", sep=""))

}


# run benjamini hochberg on tau3 p-values
pvals_adj <- p.adjust(pvalvec, method = "BH")

sigpairs <- names(which(pvals_adj < 0.05))



# Now compare myCAFs vs iCAFs in terms of the dependence of their co-expression on the modulator gene 
resmat <- matrix(0, nrow=length(sigpairs), ncol=9)
rownames(resmat) <- sigpairs
colnames(resmat) <- c("lower", "iCAF tau", "upper", "lower", "myCAF tau", "upper", "lower", "iCAF tau-myCAF tau", "upper")
for (i in 1:length(sigpairs)){
  
  genes_i <- strsplit(sigpairs[i], "_")[[1]]
  gene1 <- genes_i[1]
  gene2 <- genes_i[2]
  
  out_obj <- readRDS(paste("/path/to/realdata results/", gene1, "_", gene2, ".rds", sep=""))
  
  varcovmat <- out_obj$varcovmat
  
  iCAFs_taux <- out_obj$output.mat["tau1","coefficients"]
  myCAFs_taux <- sum(out_obj$output.mat[c("tau1","tau3"),"coefficients"])
  diff_taux <- out_obj$output.mat["tau3","coefficients"]
  
  iCAFs_taux_se <- out_obj$output.mat["tau1","std.errors"]
  myCAFs_taux_se <- c(sqrt(t(c(1,1)%*%varcovmat[c("tau1","tau3"),c("tau1","tau3")]%*%c(1,1))))
  diff_taux_se <- out_obj$output.mat["tau1","std.errors"]
  
  iCAF_interval <- iCAFs_taux + c(-1.96, 0, 1.96)*iCAFs_taux_se
  myCAF_interval <- myCAFs_taux + c(-1.96, 0, 1.96)*myCAFs_taux_se
  diff_interval <- -diff_taux + c(-1.96, 0, 1.96)*diff_taux_se

  resmat[i, ] <- c(iCAF_interval, myCAF_interval, diff_interval)
  
  
}


resmat

