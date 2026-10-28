
#PACKAGES
library("ape")
library("phylolm")

`%!in%` <- function(x, y) !(`%in%`(x, y))

# 1. load data
dat <- read.delim("/Users/Elena/Documents/irx/inverted_repeats_dtol/meta/TWLASSAQO")
tr  <- read.tree("/Users/Elena/Documents/irx/inverted_repeats_dtol/phylo/otol_grafen.nwk")

# 74 taxa not in the data but in the tree
tree_set <- setdiff(gsub(pattern = "ott", replacement = "", tr$tip.label), dat$ott_id)
tr_drop <- drop.tip(tip = paste0("ott", tree_set), phy = tr)

dat_drop <- setdiff(dat$ott_id, gsub(pattern = "ott", replacement = "", tr_drop$tip.label))
dat_2 <- unique(dat[dat$ott_id %!in% dat_drop, ])

# now we can replace the ott tip labels with names
tr_drop$tip.label <- dat_2$species[match(tr_drop$tip.label, paste0("ott", dat_2$ott_id))]
anyNA(tr_drop$tip.label)
anyDuplicated(tr_drop$tip.label)

#make into pdf so I can zoom
pdf("tr_drop.pdf", width = 20, height = 150)
plot(tr_drop, cex = 0.25, no.margin = TRUE)
dev.off()

#add colours to separate phylums into p,f,a
animal_phylum <- c("Arthropoda", "Chordata", "Mollusca", "Annelida", "Echinodermata", "Bryozoa", "Nemertea", "Nematoda", "Cnidaria", "Porifera")
plant_phylum <- c("Streptophyta", "Chlorophyta")
fungi_phylum <- c("Basidiomycota", "Ascomycota", "Mucoromycota")

phy <- dat_2$phylum[match(tr_drop$tip.label, dat_2$species)]
tipcolours <- ifelse(phy %in% animal_phylum, "red",
              ifelse(phy %in% plant_phylum, "forestgreen",
              ifelse(phy %in% fungi_phylum, "pink", "grey")))
pdf("treedrop_colour.pdf", width = 20, height = 150)
plot(tr_drop, cex = 0.25, tip.color = tipcolours, no.margin = TRUE)
legend("topleft", legend = c("Animal", "Plant", "Fungi", "Other"),
       text.col = c("red", "forestgreen", "pink", "grey"), bty= "n", cex = 3)
dev.off()

# 2. transform predictors
dat_2$log_gs  <- log10(dat_2$genome_size_bp)
dat_2$log_n50 <- log10(dat_2$scaffold_n50)
dat_2$log_nc  <- log10(dat_2$n_contigs)

# 6. model
# MCMCglmm
library(ape)
library(MCMCglmm)
#remove node labels
tr_drop$node.label <- NULL

library(data.table)
setDT(dat_2)

# must have no NAs
m0data <- dat_2[, .(species, ir_density_per_mb, genome_size_bp, scaffold_n50, n_contigs, phylum)]
m0data2 <- m0data[complete.cases(m0data),]

# subset the tree again
na_species_m0data2 <- setdiff(tr_drop$tip.label, m0data2$species)
tr_drop2 <- drop.tip(tr_drop, na_species_m0data2)

inv <- inverseA(tr_drop2, nodes = "TIPS", scale = TRUE)$Ainv
m0 <- MCMCglmm(ir_density_per_mb ~ log10(genome_size_bp) + log10(scaffold_n50) + log10(n_contigs) + phylum,
               random = ~ species,
               #family="exponential",
               ginverse = list(species = inv),
               data = m0data2, 
               nitt = 13000, 
               burnin = 3000, 
               thin = 100)
summary(m0)

# pgls?
m0 <- phylolm(ir_density_per_mb ~ genome_size_bp + scaffold_n50 + n_contigs + phylum, data = dat_2, phy = tr_drop, model = "lambda")

print(summary(m0))
