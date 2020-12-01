# Title: NMS Analysis
# Purpose: Preliminary and final analysis for triennial and annual NMFS survey NMS ordinations
# Date Created: 10/16/2020

###########################################################################################################################
# Load libraries ----
library(vegan)
library(MASS)
library(dplyr)
library(tidyr)
library(fossil)
library(tibble)
library(plyr)
library(caret)
library(MVA)
library(readxl)

###########################################################################################################################
# Load matrices ----
setwd('C:/Users/howar/Documents/Oregon State/ORnearshore_groundfish/code')
species_matrix_t <- read.csv("../data/NMFS_data/species_matrix_t.csv",
                             header = T,
                             row.names = 1,
                             check.names = F) # Triennial
env_matrix_t <- read.csv("../data/NMFS_data/env_matrix_t.csv",
                         header = T,
                         row.names = 1,
                         check.names = F)
species_matrix_a <- read.csv("../data/NMFS_data/species_matrix_a.csv",
                             header = T,
                             row.names = 1,
                             check.names = F) # Annual
env_matrix_a <- read.csv("../data/NMFS_data/env_matrix_a.csv",
                         header = T,
                         row.names = 1,
                         check.names = F)

###########################################################################################################################
# Split into two data sets ----
# Separate matrices with trawl_id row names and environmental/location columns
# Use the matrices created in NMS_matrices script, these contain log(x+1) transformed CPUE for each species and tow

# Turn matrices into dataframes and make first column the row names
# Annual species
species_matrix_a <- as.data.frame(species_matrix_a)
rownames(species_matrix_a) <- species_matrix_a[, 1]
species_matrix_a <- species_matrix_a[, -1]

# Triennial species
species_matrix_t <- as.data.frame(species_matrix_t)
rownames(species_matrix_t) <- species_matrix_t[, 1]
species_matrix_t <- species_matrix_t[, -1]

# Annual environmental data
env_matrix_a <- as.data.frame(env_matrix_a)
rownames(env_matrix_a) <- env_matrix_a[, 1]
env_matrix_a <- env_matrix_a[, -1]

# Triennial environmental data
env_matrix_t <- as.data.frame(env_matrix_t)
rownames(env_matrix_t) <- env_matrix_t[, 1]
env_matrix_t <- env_matrix_t[, -1]

###########################################################################################################################
# Add richness and Shannon diversity index to environmental matrix from species matrix ----
species_matrix_a1 <- species_matrix_a
species_matrix_t1 <- species_matrix_t
env_matrix_a$richness <- specnumber(species_matrix_a)
env_matrix_t$richness <- specnumber(species_matrix_t)
env_matrix_a$diversity <- diversity(species_matrix_a)
env_matrix_t$diversity <- diversity(species_matrix_t)

###########################################################################################################################
# Annual NMS ----
# ****Preliminary analysis ----
# Randomly sample 50% of matrix
smp_size <- floor(0.5 * nrow(species_matrix_a))
set.seed(361)
annual_ind <- sample(seq_len(nrow(species_matrix_a)), size = smp_size)
practice_a <- species_matrix_a[annual_ind,]
final <- species_matrix_a[-annual_ind,]
smp_env <- floor(0.5 * nrow(env_matrix_a))
set.seed(361)
env_matrix_a2 <- sample(seq_len(nrow(env_matrix_a)), size = smp_env)
practice_env <- env_matrix_a[env_matrix_a2,]
final_env <- env_matrix_a[-env_matrix_a2,]

# Create dissimilarity matrix and check stress
practice_a <- as.data.frame(practice_a)
mds_nulla <- isoMDS(practice_a,         # uses isoMDS engine
                   tol = 1e-7,
                   trace = F)
mds_testa <- metaMDS(practice_a,           # uses monoMDS engine
                 autotransform = F,
                 wascores = T,
                 trace = F)
stressplot(mds_nulla, annual_dis)

# Visualize ordination
windows()
ordiplot(mds_nulla, type = "t")
windows()
ordiplot(mds_testa, type = "t")

#############################################################################################################################
# ****Create final NMS ----
annual_dis <- vegdist(practice_t, "bray") # Create dissimilarity matrix with best distance measure (can skip)

# Run initial NMS with 2 axes, can choose either engine
# Can use either the dissimilarity matrix above or put in the matrix itself
# ****Use monoMDS engine -----
annual_mds_mono <- metaMDS(species_matrix_a,
                           autotransform = F,
                           trymax = 999,
                           trace = 0,
                           k = 2)
# Additional run to further lower stress
annual_mds_mono <- metaMDS(species_matrix_a,
                           autotransform = F,
                           trymax = 999,
                           trace = 0,
                           k = 2,
                           previous.best = annual_mds_mono)

# ****Use isoMDS engine ----
annual_mds_iso <- metaMDS(species_matrix_a,
                          autotransform = F,
                          zerodist = "add",
                          engine = "isoMDS",
                          trymax = 999,
                          trace = 0,
                          k = 2)
# Additional run to further lower stress
annual_mds_iso <- metaMDS(species_matrix_a,
                          autotransform = F,
                          zerodist = "add",
                          engine = "isoMDS",
                          trymax = 999,
                          trace = 0,
                          k = 2,
                          previous.best = annual_mds_iso)

##############################################################################################################################
# Fit environmental data to the ordination ----
annual_env_mds <- envfit(annual_mds_mono,
                         env_matrix_a,
                         permu = 999)
annual_env_mds

save(ef_a, file = "enviro_fit_annual")

# Rotate the ordination to depth
annual_mds_mono <- with(env_matrix_a,
                   MDSrotate(annual_mds_mono,
                             depth_m))

save(annual_mds, file = "annual_mds")

#########################################################################################################################
## plot the final ordination ----
# upload index with species groupings
index_annual <- read.csv("../data/NMFS_data/index_annual.csv", header = T)

# Create an empty plot
# Plot the species labels with ordipointlabel and then change the point colors to the species groupings
windows()
par(family = "serif")
ordiplot(annual_mds_mono,
         type = "n",
         xlim = c(-.2, .5),
         ylim = c(-1.2, 2.2),
         xlab = "Axis 1",
         ylab = "Axis 2",
         main = "Annual Survey")
#ordipointlabel(annual_mds, display = "spec", cex = 0.9, col = "black", add = T)
points(annual_mds_mono,
       display = "spec",
       select = index_annual$group == "other",
       cex = 1.7, col = "darkorchid4",
       bg = "darkorchid4",
       pch = 25)
points(annual_mds_mono,
       display = "spec",
       select = index_annual$group == "roundfish",
       cex = 1.9,
       col = "darkgreen",
       pch = 18)
points(annual_mds_mono,
       display = "spec",
       select = index_annual$group == "elasmobranch",
       cex = 1.7,
       col = "goldenrod4",
       pch = 17)
points(annual_mds_mono,
       display = "spec",
       select = index_annual$group == "rockfish",
       cex = 1.7,
       col = "navy",
       pch = 19)
points(annual_mds_mono,
       display = "spec",
       select = index_annual$group == "flatfish",
       cex = 1.7,
       col = "firebrick4",
       pch = 15)
# with environmental variables
with(env_matrix_a,
     ordisurf(annual_mds_mono,
              depth_m,
              add = T,
              col = "green4",
              cex = 4,
              labcex = .7)) # depth contours
plot(envfit(annual_mds_mono,
            env_matrix_a[ ,8:9]),
     col = "firebrick4",
     cex = 0.9) # diversity and richness
point_colors <- c("firebrick4", "navy", "darkgreen", "goldenrod4", "darkorchid4")
legend("bottomleft",
       legend = c("Flatfish", "Rockfish", "Roundfish", "Elasmobranch", "Other"),
       pch = c(15, 19, 18, 17, 25),
       col = point_colors,
       pt.bg = "darkorchid4",
       bty = "n",
       pt.cex = 1.9,
       cex = 0.9,
       inset = c(0.01, 0.05))

# Move labels

##############################################################################################################################
# Multi-response permutation procedure (MRPP)
# Create categories for climate indices
env_matrix_a$pdo <- 1 * (env_matrix_a$PDO > 0)
env_matrix_a$npgo <- 1 * (env_matrix_a$NPGO > 0)

# PDO
annual_mrpp_pdo <- with(env_matrix_a,
                        mrpp(species_matrix_a,
                             pdo,
                             distance = "bray"))
layout(matrix(1:2, nr = 1))
annual_mrpp_pdo

# NPGO
annual_mrpp_npgo <- with(env_matrix_a,
                         mrpp(species_matrix_a,
                              npgo,
                              distance = "bray"))
layout(matrix(1:2, nr = 1))
annual_mrpp_npgo

# Year
annual_mrpp_year <- with(env_matrix_a,
                         mrpp(species_matrix_a,
                              year,
                              distance = "bray"))
layout(matrix(1:2, nr = 1))
annual_mrpp_year

##############################################################################################################################
##############################################################################################################################

# Triennial NMS
# ****Preliminary analysis ----
# Randomly sample 50% of matrix
smp_size_t <- floor(0.5 * nrow(species_matrix_t))
set.seed(361)
triennial_ind <- sample(seq_len(nrow(species_matrix_t)), size = smp_size_t)
practice_t <- species_matrix_t[triennial_ind, ]
final_t <- species_matrix_t[-triennial_ind, ]
smp_env_t <- floor(0.5 * nrow(env_matrix_t))
set.seed(361)
env_matrix_t2 <- sample(seq_len(nrow(env_matrix_t)), size = smp_env_t)
practice_env_t <- env_matrix_t[env_matrix_t2, ]
final_env_t <- env_matrix_t[-env_matrix_t2, ]

# Create dissimilarity matrix and check stress
# Use null model
practice_t <- as.data.frame(practice_t)
mds_nulla <- isoMDS(practice_t,         # uses isoMDS engine
                    tol = 1e-7,
                    trace = F)
mds_testt <- metaMDS(practice_t,           # uses monoMDS engine
                     autotransform = F,
                     wascores = T,
                     trace = F)
stressplot(mds_nullt, triennial_dis)

# Visualize ordination
windows()
ordiplot(mds_nullt, type = "t")
windows()
ordiplot(mds_testt, type = "t")

#############################################################################################################################
# ****Create final NMS ----
triennial_dis <- vegdist(practice_t, "bray") # Create dissimilarity matrix with best distance measure (can skip)

# Run initial NMS with 2 axes, can choose either engine
# Can use either the dissimilarity matrix above or put in the matrix itself
# ****Use monoMDS engine -----
triennial_mds_mono <- metaMDS(species_matrix_t,
                              autotransform = F,
                              trymax = 999,
                              trace = 0,
                              k = 2)
# Additional run to further lower stress
triennial_mds_mono <- metaMDS(species_matrix_t,
                              autotransform = F,
                              trymax = 999,
                              trace = 0,
                              k = 2,
                              previous.best = triennial_mds_mono)

# ****Use isoMDS engine ----
triennial_mds_iso <- metaMDS(species_matrix_t,
                             autotransform = F,
                             zerodist = "add",
                             engine = "isoMDS",
                             trymax = 999,
                             trace = 0,
                             k = 2)
# Additional run to further lower stress
triennial_mds_iso <- metaMDS(species_matrix_t,
                             autotransform = F,
                             zerodist = "add",
                             engine = "isoMDS",
                             trymax = 999,
                             trace = 0,
                             k = 2,
                             previous.best = triennial_mds_iso)

##############################################################################################################################
# Fit environmental data to the ordination ----
ef_t <- envfit(triennial_mds,
               env_matrix_t,
               permu = 999)
ef_t

save(ef_t, file = "enviro_fit_triennial")

# Rotate the ordination to depth
triennial_mds <- with(env_matrix_t,
                      MDSrotate(triennial_mds,
                                depth_m))

save(triennial_mds, file = "triennial_mds")

##############################################################################################################################
# Plot the final ordination ----
# Import index with species groupings
index_triennial <- read.csv("../data/NMFS_data/index_triennial.csv", header = T)

# Create an empty plot
dev.new(width = 650,
        height = 550,
        unit = "px")
par(family = "serif")
ordiplot(triennial_mds,
         type = "n",
         xlim = c(-1.2, 2.2),
         ylim = c(-1.2, 2.3),
         xlab = "Axis 1",
         ylab = "Axis 2",
         main = "Triennial Survey")
# ordipointlabel(triennial_mds, display = "spec", cex = 0.9, col = "black", add = T)
points(triennial_mds,
       display = "spec",
       select = index_triennial$group == "other",
       cex = 1.7,
       col = "darkorchid4",
       bg = "darkorchid4",
       pch = 25)
points(triennial_mds,
       display = "spec",
       select = index_triennial$group == "roundfish",
       cex = 1.9,
       col = "darkgreen",
       pch = 18)
points(triennial_mds,
       display = "spec",
       select = index_triennial$group == "elasmobranch",
       cex = 1.7,
       col = "goldenrod4",
       pch = 17)
points(triennial_mds,
       display = "spec",
       select = index_triennial$group == "rockfish",
       cex = 1.7,
       col = "navy",
       pch = 19)
points(triennial_mds,
       display = "spec",
       select = index_triennial$group == "flatfish",
       cex = 1.7,
       col = "firebrick4",
       pch = 15)
# with environmental variables
with(env_matrix_t,
     ordisurf(triennial_mds,
              depth_m,
              add = T,
              col = "green4",
              cex = 4,
              labcex = .7)) # depth contours
plot(envfit(triennial_mds, env_matrix_t[, 8:9]),
     col = "firebrick4",
     cex = 0.01,
     labels = F) # diversity and richness
text(locator(1), "richness", cex = 0.9, col = "firebrick4")
text(locator(1), "diversity", cex = 0.9, col = "firebrick4")
point_colors <- c("firebrick4", "navy", "darkgreen", "goldenrod4", "darkorchid4")
#legend("bottomleft", legend = c("Flatfish", "Rockfish", "Roundfish", "Elasmobranch", "Other"), pch = c(15, 19, 18, 17, 25),
#col = point_colors, pt.bg = "darkorchid4", bty = "n", pt.cex = 1, cex = 0.9, inset = c(0.75, 0.75))

##############################################################################################################################
# Multi-response permutation procedure (MRPP)
# Create categories for climate indices
env_matrix_t$pdo <- 1 * (env_matrix_t$PDO > 0)
env_matrix_t$npgo <- 1 * (env_matrix_t$NPGO > 0)

# PDO
triennial_mrpp_pdo <- with(env_matrix_t,
                       mrpp(species_matrix_t,
                            pdo,
                            distance = "bray"))
layout(matrix(1:2, nr = 1))
triennial_mrpp_pdo

# NPGO
triennial_mrpp_npgo <- with(env_matrix_t,
                        mrpp(species_matrix_t,
                             npgo,
                             distance = "bray"))
layout(matrix(1:2, nr = 1))
triennial_mrpp_npgo

# Year
triennial_mrpp_year <- with(env_matrix_t,
                        mrpp(species_matrix_t,
                             year,
                             distance = "bray"))
layout(matrix(1:2, nr = 1))
triennial_mrpp_year
