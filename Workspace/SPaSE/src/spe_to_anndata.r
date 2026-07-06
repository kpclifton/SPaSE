#convert spe to anndata and save as h5ad

###SET UP ###
#Sys.setenv(RETICULATE_PYTHON = "/home/kalen/ENTER/envs/spase/bin/python")
#this isn't working correctly with VScode R interactive so run the following in the terminal before
#export LD_PRELOAD=/home/kalen/ENTER/envs/spase/lib/libstdc++.so.6
#export LD_LIBRARY_PATH=/home/kalen/ENTER/envs/spase/lib:$LD_LIBRARY_PATH
#export RETICULATE_PYTHON=/home/kalen/ENTER/envs/spase/bin/python
#R

#continuosly getting errors from the above spase env so had to make a new spase-clean environment with modules from only conda-forge instead of a mix of conda-forge and pip

library(SEraster) 
library(STcompare)
library(SpatialExperiment)
library(reticulate) 
#use_condaenv("spase", required = TRUE) 
py_config()

anndata <- import("anndata") 
np <- import("numpy")

data("speKidney") 
head(speKidney)

X <- t(as.matrix(assay(speKidney$A, "counts")))
#obs <- as.data.frame(colData(speKidney$A))
#var <- as.data.frame(rowData(speKidney$A))
#ad <- anndata$AnnData(X = X, obs = obs, var = var)
obs <- data.frame(array_row = spatialCoords(speKidney$A)[, 'y'], array_col = spatialCoords(speKidney$A)[, 'x'])
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(speKidney$A))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/speKidney/speKidneyA.h5ad") 
ad$write_h5ad(filepath)


X <- t(as.matrix(assay(speKidney$B, "counts")))
#obs <- as.data.frame(colData(speKidney$B))
#var <- as.data.frame(rowData(speKidney$B))
#ad <- anndata$AnnData(X = X, obs = obs, var = var)
obs <- data.frame(array_row = spatialCoords(speKidney$B)[, 'y'], array_col = spatialCoords(speKidney$B)[, 'x'])
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(speKidney$B))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/speKidney/speKidneyB.h5ad") 
ad$write_h5ad(filepath)

X <- t(as.matrix(assay(speKidney$C, "counts")))
#obs <- as.data.frame(colData(speKidney$C))
#var <- as.data.frame(rowData(speKidney$C))
#ad <- anndata$AnnData(X = X, obs = obs, var = var)
obs <- data.frame(array_row = spatialCoords(speKidney$C)[, 'y'], array_col = spatialCoords(speKidney$C)[, 'x'])
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(speKidney$C))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/speKidney/speKidneyC.h5ad") 
ad$write_h5ad(filepath)

#create 99 genes based on expression of Gene in dataset A
mean_mu <- mean(assay(speKidney$A))
sd_sigma <- sd(assay(speKidney$A))          
target_sd <- sd_sigma / 2

set.seed(42)
n_genes <- 99
n_cells <- length(assay(speKidney$A))

outmatA <- matrix(NA_real_, nrow=n_genes, ncol=n_cells)
for (i in 1:n_genes) {
samp <- rnorm(n_cells, mean=mean_mu, sd=target_sd)
samp[samp < 0.1] <- 0.1
outmatA[i, ] <- round(samp, 6)
}

outmatA
assay(speKidney$A)

# create new gene names
gene_names <- sprintf("Gene%03d", 1:n_genes)
rownames(outmatA) <- gene_names
colnames(outmatA) <- colnames(assay(speKidney$A))

tempAssay <- rbind(assay(speKidney$A), outmatA)

assay(speKidney$A, "counts", withDimnames=FALSE) <- tempAssay



#repeat with rasterization to use the synthensizing healthy sample part of SPaSE which is expecting rows and columns with integer indices
rastKidney <- SEraster::rasterizeGeneExpression(speKidney,
                assay_name = 'counts', resolution = 0.2,
                square = FALSE)

# After rasterization, the output is a SpatialExperiment object, but the spatial units are now raster pixels rather than individual cells or spots. The assay has been converted to pixelval, and additional metadata (num_cell, cellID_list, geometry) records which original cells contributed to each pixel.
head(rastKidney)

X <- t(as.matrix(assay(rastKidney$A, "pixelval")))
obs <- data.frame(array_row = as.integer(spatialCoords(rastKidney$A)[, 'x']*10), array_col = as.integer(spatialCoords(rastKidney$A)[, 'y']*10))
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(rastKidney$A))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/rastKidney/rastKidneyA.h5ad") 
ad$write_h5ad(filepath)


X <- t(as.matrix(assay(rastKidney$B, "pixelval")))
obs <- data.frame(array_row = as.integer(spatialCoords(rastKidney$B)[, 'x']*10), array_col = as.integer(spatialCoords(rastKidney$B)[, 'y']*10))
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(rastKidney$B))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/rastKidney/rastKidneyB.h5ad") 
ad$write_h5ad(filepath)

X <- t(as.matrix(assay(rastKidney$C, "pixelval")))
obs <- data.frame(array_row = as.integer(spatialCoords(rastKidney$C)[, 'x']*10), array_col = as.integer(spatialCoords(rastKidney$C)[, 'y']*10))
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(rastKidney$C))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/rastKidney/rastKidneyC.h5ad") 
ad$write_h5ad(filepath)


simRanPatternRasts[[1]]
simRanPatternRasts[[7]]

X <- t(as.matrix(assay(simRanPatternRasts[[1]])))
#obs <- data.frame(sample_id = colData(simRanPatternRasts[[1]])[,'sample_id'])
#var <- as.data.frame(rowData(simRanPatternRasts[[1]]))
#ad <- anndata$AnnData(X = X, obs = obs, var = var)
obs <- data.frame(array_row = as.integer((spatialCoords(simRanPatternRasts[[1]])[, 'x'] + 3)*10), array_col = as.integer((spatialCoords(simRanPatternRasts[[1]])[, 'y'] + 3) *10))
rownames(obs) <- rownames(spatialCoords(simRanPatternRasts[[1]]))
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(obs)
rownames(coords) <- rownames(spatialCoords(simRanPatternRasts[[1]]))
colnames(coords) <- colnames(spatialCoords(simRanPatternRasts[[1]]))
#coords <- as.matrix(spatialCoords(simRanPatternRasts[[1]]))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/speKidney_simRan/simRanPatternRasts1.h5ad") 
ad$write_h5ad(filepath)

for (i in 2:length(simRanPatternRasts)){
X <- t(as.matrix(assay(simRanPatternRasts[[i]])))
#obs <- data.frame(sample_id = colData(simRanPatternRasts[[i]])[,'sample_id'])
#var <- as.data.frame(rowData(simRanPatternRasts[[i]]))
#ad <- anndata$AnnData(X = X, obs = obs, var = var)
obs <- data.frame(array_row = as.integer((spatialCoords(simRanPatternRasts[[i]])[, 'x'] + 3)*10), array_col = as.integer((spatialCoords(simRanPatternRasts[[i]])[, 'y']+3)*10))
rownames(obs) <- rownames(spatialCoords(simRanPatternRasts[[i]]))
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(obs)
rownames(coords) <- rownames(spatialCoords(simRanPatternRasts[[i]]))
colnames(coords) <- colnames(spatialCoords(simRanPatternRasts[[i]]))
#coords <- as.matrix(spatialCoords(simRanPatternRasts[[i]]))
ad$obsm["spatial"] <- coords
filepath <- path.expand(paste0("~/SPaSE/Data/speKidney_simRan/simRanPatternRasts", i, ".h5ad")) 
ad$write_h5ad(filepath)
}

# Load MERFISH datasets from Zenodo -----------------------------------------------

read_csv_gz_url <- function(x, ...) {
  u  <- url(x, open = "rb")
  gz <- gzcon(u, text = TRUE)
  on.exit({
    try(close(gz), silent = TRUE)
    try(close(u),  silent = TRUE)
  }, add = TRUE)
  
  read.csv(gz, row.names = 1, ...)
}

s2r2_url <- "https://zenodo.org/records/10724029/files/STalign_S2R2.csv.gz?download=1"
s2r3_url <- "https://zenodo.org/records/10724029/files/STalign_S2R3_to_S2R2.csv.gz?download=1"

target <- read_csv_gz_url(s2r2_url)
source <- read_csv_gz_url(s2r3_url)

dim(target)
head(target)
pos_target <- target[, c('x', 'y')]
gene_target <- target[, 3:dim(target)[2]]

## check class of counts matrix
class(gene_target)

## convert classes
target_sparse <- as(t(gene_target), "dgCMatrix")

## format into SpatialExperiment
spe_target <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = target_sparse),
  spatialCoords = as.matrix(pos_target)
)

dim(source)
head(source)
pos_source <- source[, c('STalign_x', 'STalign_y')]
colnames(pos_source) <- c("x", "y")
gene_source <- source[, 7:dim(source)[2]]

## check class of counts matrix
class(gene_source)

## convert classes
source_sparse <- as(t(gene_source), "dgCMatrix")

## format into SpatialExperiment
spe_source <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = source_sparse),
  spatialCoords = as.matrix(pos_source)
)

spe_list <- list(target = spe_target, source = spe_source)
resolution = 200
output <- SEraster::rasterizeGeneExpression(spe_list, resolution = resolution, fun = "mean", BPPARAM = BiocParallel::MulticoreParam())




load(file = "~/ST_compare/data/merfish_data/moransI_20250702.RData")

svg_int <- intersect(svg_merfish_target, svg_merfish_source)
svg_uni <- union(svg_merfish_target, svg_merfish_source)
non_svg <- rownames(output$target)[!(rownames(output$target) %in% svg_int)]

length(svg_int)

output$target[svg_int,]

X <- t(as.matrix(assay(output$target[svg_int,])))
obs <- data.frame(array_row = as.integer(spatialCoords(output$target)[, 'y']), array_col = as.integer(spatialCoords(output$target)[, 'x']))
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(output$target))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/speMERFISH/target_S2R2.h5ad") 
ad$write_h5ad(filepath)


X <- t(as.matrix(assay(output$source[svg_int,])))
obs <- data.frame(array_row = as.integer(spatialCoords(output$source)[, 'y']), array_col = as.integer(spatialCoords(output$source)[, 'x']))
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(output$source))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/speMERFISH/source_S2R3.h5ad") 
ad$write_h5ad(filepath)


#load Visium data

library(STcompare)
library(SpatialExperiment)
#library(MERINGUE)
library(rhdf5)
library(ggplot2)
library(dplyr)
#library(patchwork)
library(Matrix)
library(BiocGenerics)

# The Visium mouse kidney datasets can be found Zenodo: https://doi.org/10.5281/zenodo.17676991
# IL3 - ischemia acute kidney injury dataset 
# NL3 - control dataset 

# AKI counts ###################################################################
zenodo_url <- "https://zenodo.org/records/19074288/files/IL3_filtered_feature_bc_matrix.h5?download=1"
temp_h5 <- tempfile(fileext = ".h5")
download.file(zenodo_url, temp_h5, mode = "wb")

# this is a list of the internal groups and datasets stored in the HDF5 file 
# we are interested in the "matrix/barcodes" dataset 
rhdf5::h5ls(temp_h5)

# Read the CSC-encoded matrix components from the HDF5 file
file_barcodes <- as.character(rhdf5::h5read(temp_h5, "matrix/barcodes"))
barcodes_matrix <- rhdf5::h5read(temp_h5, "matrix") # stores the sparse matrix components 

# convert CSC encoding into an R dgCMatrix
aki_counts <- Matrix::sparseMatrix(
  dims = barcodes_matrix$shape,
  i = as.numeric(barcodes_matrix$indices),
  p = as.numeric(barcodes_matrix$indptr),
  x = as.numeric(barcodes_matrix$data),
  index1 = FALSE
)
colnames(aki_counts) <- file_barcodes
rownames(aki_counts) <- barcodes_matrix[["features"]]$name

unlink(temp_h5)
head(aki_counts)

zenodo_url <- "https://zenodo.org/records/19074288/files/NL3_filtered_feature_bc_matrix.h5?download=1"
temp_h5 <- tempfile(fileext = ".h5")
download.file(zenodo_url, temp_h5, mode = "wb")

# Read the CSC-encoded matrix components from the HDF5 file
file_barcodes <- as.character(rhdf5::h5read(temp_h5, "matrix/barcodes"))
barcodes_matrix <- rhdf5::h5read(temp_h5, "matrix") # stores the sparse matrix components 

# convert CSC encoding into an R dgCMatrix
control_counts <- Matrix::sparseMatrix(
  dims = barcodes_matrix$shape,
  i = as.numeric(barcodes_matrix$indices),
  p = as.numeric(barcodes_matrix$indptr),
  x = as.numeric(barcodes_matrix$data),
  index1 = FALSE
)
colnames(control_counts) <- file_barcodes
rownames(control_counts) <- barcodes_matrix[["features"]]$name

unlink(temp_h5)
head(control_counts)

zenodo_url <- "https://zenodo.org/records/19074288/files/IL3_tissue_positions.csv?download=1"
aki_pos <- read.csv(zenodo_url, header = TRUE, stringsAsFactors = FALSE)

zenodo_url <- "https://zenodo.org/records/19074288/files/NL3_tissue_positions.csv?download=1"
ctrl_pos <- read.csv(zenodo_url, header = TRUE, stringsAsFactors = FALSE)

head(ctrl_pos)
head(aki_pos)

dim(aki_counts)
dim(aki_pos)

dim(control_counts)
dim(ctrl_pos)

# limit the counts to the spots with positions and make sure they are in the same order
aki_counts <- aki_counts[, colnames(aki_counts) %in% aki_pos$barcode]
control_counts <- control_counts[, colnames(control_counts) %in% ctrl_pos$barcode]

# number of columns in counts is the same number of rows in pos
dim(aki_counts)
dim(aki_pos)

dim(control_counts)
dim(ctrl_pos)

# set the barcode column as the rownames
rownames(ctrl_pos) <- ctrl_pos$barcode
ctrl_pos <- ctrl_pos[,-1]
colnames(ctrl_pos) <- c("x", "y")

rownames(aki_pos) <- aki_pos$barcode
aki_pos <- aki_pos[,-1]
colnames(aki_pos) <- c("x", "y")

ctrl_pos$group <- "Control"
aki_pos$group  <- "AKI"

# visualize both control and AKI positions on the same plot
df <- rbind(ctrl_pos, aki_pos)
pos_plot <- ggplot(df, aes(x = x, y = y, color = group)) +
  geom_point() +
  scale_color_manual(values = c("Control" = "blue", "AKI" = "red")) +
  coord_fixed() +
  labs(x = "x", y = "y", color = "Group") +
  theme_classic()
pos_plot

# rotate the tissue by 90-degrees for consistency with the paper 
ctrl_pos_rot <- transform(ctrl_pos, x = y, y = -x + max(ctrl_pos$x))
aki_pos_rot <- transform(aki_pos, x = y, y = -x + max(aki_pos$x))

df_rot <- rbind(ctrl_pos_rot, aki_pos_rot)
pos_rot_plot <- ggplot(df_rot, aes(x = x, y = y, color = group)) +
  geom_point(size=0.5) +
  scale_color_manual(values = c("Control" = "blue", "AKI" = "red")) +
  coord_fixed() +
  labs(x = "x", y = "y", color = "Group") +
  theme_classic()
pos_rot_plot

# read the csv which has the AKI positions aligned to control via STalign
zenodo_url <- "https://zenodo.org/records/19486091/files/aki_region_onehot_STalign_to_ctrl_region_onehot_affine_only.csv.gz?download=1"
aki_STalign <- read_csv_gz_url(zenodo_url)

# store the aligned positions
aki_pos_aligned <- aki_STalign[,c("aligned_x", "aligned_y")]

#rename columns and add group for consistency with ctrl_pos_rot
colnames(aki_pos_aligned) <- c("x", "y")
aki_pos_aligned$group  <- "AKI"
# rotate the tissue by 90-degrees for consistency with the paper 
aki_pos_rot <- transform(aki_pos_aligned, x = y, y = -x + max(aki_pos_aligned$x))
aki_pos_aligned <- aki_pos_rot

df_aligned <- rbind(ctrl_pos_rot, aki_pos_aligned)
pos_aligned_plot <- ggplot(df_aligned, aes(x = x, y = y, color = group)) +
  geom_point(size=0.5) +
  scale_color_manual(values = c("Control" = "blue", "AKI" = "red")) +
  coord_fixed() +
  labs(x = "x", y = "y", color = "Group") +
  theme_classic()
pos_aligned_plot

AKI_ctrl_SE <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = control_counts),
  spatialCoords = as.matrix(ctrl_pos_rot[,1:2]),
)

AKI_aki_SE <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = aki_counts),
  spatialCoords = as.matrix(aki_pos_aligned[,1:2])
)

input <- list(AKI_ctrl = AKI_ctrl_SE, AKI_aki = AKI_aki_SE)
rast <- SEraster::rasterizeGeneExpression(
  input,
  resolution = 5,
  fun        = "sum",
  square     = FALSE,
  assay_name = 'counts'
)

# add a CPM normalization assay to the rasterized SpatialExperiment
assay(rast$AKI_ctrl, "CPM") <- Matrix::t(Matrix::t(assay(rast$AKI_ctrl))/Matrix::colSums(assay(rast$AKI_ctrl)))*1e6
assay(rast$AKI_aki, "CPM") <- Matrix::t(Matrix::t(assay(rast$AKI_aki))/Matrix::colSums(assay(rast$AKI_aki)))*1e6


load(file = "~/ST_compare/data/kidney_data/moransI_svg01per_20250817.RData")

svg_ctrl <- moransI_ctrl %>% filter(minPercentCells > 0.01) %>% rownames()
svg_irl <- moransI_irl %>% filter(minPercentCells > 0.01) %>% rownames()

svg_int_kid <- intersect(svg_ctrl, svg_irl)
svg_uni_kid <- union(svg_ctrl, svg_irl)

# make anndata
X <- t(as.matrix(assay(AKI_ctrl_SE[svg_int_kid,])))
obs <- data.frame(array_row = spatialCoords(AKI_ctrl_SE)[, 'y'], array_col = spatialCoords(AKI_ctrl_SE)[, 'x'])
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(AKI_ctrl_SE))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/AKI_Kidney/AKI_ctrl_SE.h5ad") 
ad$write_h5ad(filepath)


X <- t(as.matrix(assay(AKI_aki_SE[svg_int_kid,])))
obs <- data.frame(array_row = spatialCoords(AKI_aki_SE)[, 'y'], array_col = spatialCoords(AKI_aki_SE)[, 'x'])
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(AKI_aki_SE))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/AKI_Kidney/AKI_aki_SE.h5ad") 
ad$write_h5ad(filepath)

#make anndata rasterized and CPM
X <- t(as.matrix(assay(rast$AKI_ctrl[svg_int_kid,], "CPM")))
obs <- data.frame(array_row = as.integer(spatialCoords(rast$AKI_ctrl)[, 'y']), array_col = as.integer(spatialCoords(rast$AKI_ctrl)[, 'x']))
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(rast$AKI_ctrl))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/AKI_Kidney_rast/AKI_ctrl.h5ad") 
ad$write_h5ad(filepath)


X <- t(as.matrix(assay(rast$AKI_aki[svg_int_kid,], "CPM")))
obs <- data.frame(array_row = as.integer(spatialCoords(rast$AKI_aki)[, 'y']), array_col = as.integer(spatialCoords(rast$AKI_aki)[, 'x']))
ad <- anndata$AnnData(X = X, obs = obs)
coords <- as.matrix(spatialCoords(rast$AKI_aki))
ad$obsm["spatial"] <- coords
filepath <- path.expand("~/SPaSE/Data/AKI_Kidney_rast/AKI_aki.h5ad") 
ad$write_h5ad(filepath)
