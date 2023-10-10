# Criacao do objeto phyloseq para plotagem de graficos

setwd('./4 Semestre/IP III/Dados Projeto/')

library(readxl)
library(phylosmith)
library(phyloseq)
library(tibble)
library(ape)

# Importacao dos arquivos

seqtab_data <- read_excel("./16S_seqtab.xlsx")
taxa_data <- read_excel("./16S_taxa_data.xlsx")
metadata <- read_excel("./16S_metadata.xlsx")

# Correspondencia da coluna ASV

correspondencia <- match(seqtab_data$ASV, taxa_data$ASV)

# Correspondencia nos dataframes

seqtab_data <- seqtab_data[correspondencia, ]
taxa_data <- taxa_data[correspondencia, ]

seqtab_data$ASV <- paste0("ASV", 1:nrow(seqtab_data))
taxa_data$ASV <- paste0("ASV", 1:nrow(taxa_data))

seqtab_data <- tibble::column_to_rownames(seqtab_data, "ASV")
taxa_data <- tibble::column_to_rownames(taxa_data, "ASV")
metadata <- tibble::column_to_rownames(metadata, "Samples")

OTU <- otu_table(as.matrix(seqtab_data), taxa_are_rows=TRUE)
TAX <- tax_table(as.matrix(taxa_data))
metada <- as.matrix(metadata)

ps <- phyloseq(OTU,TAX,sample_data(metadata))

# Plot da árvore filogenética

random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
plot(random_tree)

# Juncao do objeto phyloseq com o metada e a arvore filogenetica

physeq1 = merge_phyloseq(ps, metadata, random_tree)
physeq1

# Junção do objeto OTU com o metada,TAX e a árvore filogenética

physeq2 = phyloseq(OTU, sample_data(metadata), TAX, phy_tree(random_tree))
physeq2

# O que significa ser identico?
identical(physeq1, physeq2)

FSr  = transform_sample_counts(physeq1, function(x) x / sum(x) )
FSfr = filter_taxa(FSr, function (x) sum(x) > .1, TRUE)

pdf(file = "D:/CNPEM/LNBR/PB/pictures/Beta_D_M_ITS.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6) # The height of the plot in inches

taxa_abundance_bars(taxa_filter(Fsr, frequency = 0.8), classification = 'Phylum', treatment = c('Plastic'), transformation = 'mean')


ntaxa(ps)
nsamples(ps)
sample_names(ps)
rank_names(ps)

table(tax_table(ps)[,'Phylum'],exclude=NULL)

richness <- estimate_richness(ps)





