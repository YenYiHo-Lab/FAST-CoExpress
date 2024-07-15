library(glmGamPoi)
library(R.utils)
library(sctransform)
library(Seurat)


#gunzip("/path/to/GSM5354525_CID44041.tar.gz")


data_dir <- "/path/to/CID4471/"


matrix_path <- paste(data_dir,"count_matrix_sparse.mtx",sep="")
features_path <- paste(data_dir,"count_matrix_genes.tsv",sep="")
barcodes_path <- paste(data_dir,"count_matrix_barcodes.tsv",sep="")
metadata_path <- paste(data_dir,"metadata.csv",sep="")


# reading the sparse matrix, features, and barcodes
counts <- Matrix::readMM(file = matrix_path) # Requires the Matrix package
features <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)

rownames(counts) <- features$V1
colnames(counts) <- barcodes$V1


# adding metadata
metadat <- read.csv(metadata_path, header=TRUE, row.names=1)
metadat <- metadat[,c("orig.ident", "celltype_minor")]


# subsetting to iCAF myCAF
cafcells <- metadat[metadat$celltype_minor %in% c("CAFs myCAF-like", "CAFs MSC iCAF-like"), ]
cafcells <- rownames(cafcells)


# pre filtering to only focus on genes that are at least moderately expressed
# in myCAF and iCAF cells
subcounts <- counts[, cafcells]

subdf <- as.matrix(subcounts)

rm(subcounts)
gc()

usegene_names <- rownames(subdf)
usecell_names <- colnames(subdf)

subdf <- data.frame(t(unname(subdf)))
colnames(subdf) <- usegene_names
rownames(subdf) <- usecell_names
subdf$cell_type <- metadat[usecell_names,"celltype_minor"]


# gene should be expressed in a certain percentage of iCAF and myCAF cells
icaf_keep <- apply(subdf[subdf$cell_type=="CAFs MSC iCAF-like",1:(ncol(subdf)-1)], 2, function(a) mean(a>0)) > 0.15
icaf_keep <- names(icaf_keep[icaf_keep==TRUE])

mycaf_keep <- apply(subdf[subdf$cell_type=="CAFs myCAF-like",1:(ncol(subdf)-1)], 2, function(a) mean(a>0)) > 0.15
mycaf_keep <- names(mycaf_keep[mycaf_keep==TRUE])



rm(subdf)
gc()


# subset counts down to just these genes
counts <- counts[, cafcells] #counts[intersect(icaf_keep, mycaf_keep),]



# create seurat object
seurat_obj <- CreateSeuratObject(counts = counts, meta.data=metadat, min.cells=3, min.features=200)

rm(counts)
gc()


# filter by mitochondria
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
seurat_obj <- subset(seurat_obj, subset = percent.mt < 20)


# run sctransform
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE, vst.flavor = "v2")


# turn it into a useable dataframe
ourdf <- as.matrix(seurat_obj@assays$SCT$counts)
usegene_names <- rownames(ourdf)
usecell_names <- colnames(ourdf)

ourdf <- data.frame(t(unname(ourdf)))
colnames(ourdf) <- usegene_names
rownames(ourdf) <- usecell_names
ourdf$cell_type <- metadat[usecell_names,"celltype_minor"]


# remove genes which don't have 75th quantile > 3.99
icaf_small <- apply(ourdf[ourdf$cell_type=="CAFs MSC iCAF-like",1:(ncol(ourdf)-1)], 2, function(a) quantile(a, 0.75))<3.99
icaf_small <- names(icaf_small[icaf_small==TRUE])

mycaf_small <- apply(ourdf[ourdf$cell_type=="CAFs myCAF-like",1:(ncol(ourdf)-1)], 2, function(a) quantile(a, 0.75))<3.99
mycaf_small <- names(mycaf_small[icaf_small==TRUE])


ourdf <- ourdf[,!(colnames(ourdf) %in% unique(c(mycaf_small, icaf_small)))]

# remove mitochondria genes
ourdf <- ourdf[,-grep("^MT-", colnames(ourdf))]


# form the gene pairs
gene_pairs <- t(combn(colnames(ourdf)[1:(ncol(ourdf)-1)],2))


# ensure each gene pair has 70% non-overlap with each other for both iCAF and myCAF
small_overlap <- rep(0,nrow(gene_pairs))
for (i in 1:nrow(gene_pairs)){
  
  y1gene <- ourdf[ourdf$cell_type=="CAFs MSC iCAF-like",gene_pairs[i,1]]
  y2gene <- ourdf[ourdf$cell_type=="CAFs MSC iCAF-like",gene_pairs[i,2]]
  
  small_overlap[i] <- 1*(mean((y1gene >0) & (y2gene > 0))<0.3)

  y1gene <- ourdf[ourdf$cell_type=="CAFs myCAF-like",gene_pairs[i,1]]
  y2gene <- ourdf[ourdf$cell_type=="CAFs myCAF-like",gene_pairs[i,2]]  
  
  small_overlap[i] <- small_overlap[i] + 1*(mean((y1gene >0) & (y2gene > 0))<0.3)

}

gene_pairs <- gene_pairs[-which(small_overlap > 0),]



saveRDS(gene_pairs, "/path/to/save/gene_pairs")
saveRDS(ourdf, "/path/to/save/processed_counts")


