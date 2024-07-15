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

counts <- Matrix::readMM(file = matrix_path) # Requires the Matrix package
features <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)

rownames(counts) <- features$V1
colnames(counts) <- barcodes$V1

# add metadata
metadat <- read.csv(metadata_path, header=TRUE, row.names=1)
metadat <- metadat[,c("orig.ident", "celltype_minor")]



# subsetting to iCAF myCAF
cafcells <- metadat[metadat$celltype_minor %in% c("CAFs myCAF-like", "CAFs MSC iCAF-like"), ]
cafcells <- rownames(cafcells)

# initial filtering to only focus on genes that are at least moderately expressed in myCAF and iCAF cells
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

# we want the gene to be expressed in a certain percentage of iCAF and myCAF cells
icaf_keep <- apply(subdf[subdf$cell_type=="CAFs MSC iCAF-like",1:(ncol(subdf)-1)], 2, function(a) mean(a>0)) > 0.1
icaf_keep <- names(icaf_keep[icaf_keep==TRUE])

mycaf_keep <- apply(subdf[subdf$cell_type=="CAFs myCAF-like",1:(ncol(subdf)-1)], 2, function(a) mean(a>0)) > 0.1
mycaf_keep <- names(mycaf_keep[mycaf_keep==TRUE])

rm(subdf)
gc()

# now subsetting counts down to just these genes
counts <- counts[intersect(icaf_keep, mycaf_keep),]

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, meta.data=metadat, min.cells=10, min.features=200)

rm(counts)
gc()

# run sctransform
seurat_obj <- SCTransform(seurat_obj, verbose = TRUE, vst.flavor = "v2")

# turn it into our useable dataframe
ourdf <- as.matrix(seurat_obj@assays$SCT$counts)
usegene_names <- rownames(ourdf)
usecell_names <- colnames(ourdf)

ourdf <- data.frame(t(unname(ourdf)))
colnames(ourdf) <- usegene_names
rownames(ourdf) <- usecell_names
ourdf$cell_type <- metadat[usecell_names,"celltype_minor"]

# now let's restrict ourselves to just iCAFs and myCAFs
ourdf <- ourdf[ourdf$cell_type %in% c("CAFs MSC iCAF-like","CAFs myCAF-like"),]

# remove genes which don't have 75th quantile > 3.99
# we need this to hold for both iCAFs and myCAFs
icaf_small <- apply(ourdf[ourdf$cell_type=="CAFs MSC iCAF-like",1:(ncol(ourdf)-1)], 2, function(a) quantile(a, 0.75))<3.99
icaf_small <- names(icaf_small[icaf_small==TRUE])

mycaf_small <- apply(ourdf[ourdf$cell_type=="CAFs myCAF-like",1:(ncol(ourdf)-1)], 2, function(a) quantile(a, 0.75))<3.99
mycaf_small <- names(mycaf_small[icaf_small==TRUE])

ourdf <- ourdf[,!(colnames(ourdf) %in% unique(c(mycaf_small, icaf_small)))]

                     
# remove genes which have over 70% zeros
icaf_manyzeros <- apply(ourdf[ourdf$cell_type=="CAFs MSC iCAF-like",1:(ncol(ourdf)-1)], 2, function(a) mean(a==0))>0.3
icaf_manyzeros <- names(icaf_small[icaf_small==TRUE])

mycaf_manyzeros <- apply(ourdf[ourdf$cell_type=="CAFs myCAF-like",1:(ncol(ourdf)-1)], 2, function(a) mean(a==0))>0.3
mycaf_manyzeros <- names(mycaf_small[icaf_small==TRUE])

ourdf <- ourdf[,!(colnames(ourdf) %in% unique(c(icaf_manyzeros, mycaf_manyzeros)))]
                         
# remove mitochondrial genes
ourdf <- ourdf[,-grep("^MT-", colnames(ourdf))]

# make the gene pairs
gene_pairs <- t(combn(colnames(ourdf)[1:(ncol(ourdf)-1)],2))

# ensure that the gene pair doesnt have (0,0) for more than 70% of the cells
small_overlap <- rep(0,nrow(gene_pairs))
for (i in 1:nrow(gene_pairs)){
  
  y1gene <- ourdf[ourdf$cell_type=="CAFs MSC iCAF-like",gene_pairs[i,1]]
  y2gene <- ourdf[ourdf$cell_type=="CAFs MSC iCAF-like",gene_pairs[i,2]]
  
  small_overlap[i] <- 1*(mean((y1gene >0) & (y2gene > 0))<0.5)

  y1gene <- ourdf[ourdf$cell_type=="CAFs myCAF-like",gene_pairs[i,1]]
  y2gene <- ourdf[ourdf$cell_type=="CAFs myCAF-like",gene_pairs[i,2]]  
  
  small_overlap[i] <- small_overlap[i] + 1*(mean((y1gene >0) & (y2gene > 0))<0.5)

}

gene_pairs <- gene_pairs[-which(small_overlap > 0),]

saveRDS(gene_pairs, "/path/to/save/gene_pairs")
saveRDS(ourdf, "/path/to/save/processed_counts")

