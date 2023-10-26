# Criacao do objeto phyloseq para plotagem de graficos

setwd('./4 Semestre/IP III/Dados Projeto/')

library(BiocManager)
library(ggplot2)
library(readxl)
library(phylosmith)
library(phyloseq)
library(tibble)
library(ape)
library(dendextend)
library(tidyr)
library(ggtree)
library(vegan)
library(phangorn)
library(phytools)
library(geiger)
library(readr)

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

write.csv(OTU, file = "16S_otu.csv", row.names = TRUE)
write.csv(TAX, file = "16S_tax.csv", row.names = TRUE)
write_tsv(data.frame(TAX), '16S_tax.txt')

ps <- phyloseq(OTU,TAX,sample_data(metadata))

# Plot da árvore filogenética

random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
plot(random_tree)

# Juncao do objeto phyloseq com o metada e a arvore filogenetica

physeq1 = merge_phyloseq(ps, metadata, random_tree)
physeq1

# Juncao do objeto OTU com o metada,TAX e a árvore filogenética

physeq2 = phyloseq(OTU, sample_data(metadata), TAX, phy_tree(random_tree))
physeq2

# O que significa ser identico?
identical(physeq1, physeq2)

FSr  = transform_sample_counts(physeq1, function(x) x / sum(x) )
FSfr = filter_taxa(physeq1, function (x) sum(x) > .1, TRUE)

# Informacoes importantes extraidas

ntaxa(ps)
nsamples(ps)
sample_names(ps)
rank_names(ps)

table(tax_table(ps)[,'Phylum'],exclude=NULL)

# Riqueza encontrada nos dados

richness <- estimate_richness(ps)

anova.sh = aov(richness$Shannon ~ sample_data(ps)$Plastic)

TukeyHSD(anova.sh)

# Perfil filogenético

phylogeny_profile(FSfr, classification = 'Phylum', treatment = c('Plastic'), merge = TRUE, relative_abundance = TRUE)

# Perfil filogenético para o gênero

FSfr_copy <- FSfr
tax_table_nova <- FSfr_copy@tax_table
abund <- taxa_data$Genus
abund_count <- aggregate(data.frame(count = abund), list(value = abund), length)
nova_coluna

f1  <- filterfun_sample(topk(20))
wh1 <- genefilter_sample(OTU, f1, A=1)
OTU_trim <- prune_taxa(wh1, OTU)

ps_trim <- phyloseq(OTU_trim, sample_data(metadata), TAX)

phylogeny_profile(ps_trim, classification = 'Genus', treatment = c('Plastic'), merge = TRUE, relative_abundance = TRUE)

# Plot da diversidade alpha

alpha_meas = c("Observed", "Shannon", "Simpson")

p <- plot_richness(ps, "Plastic", "Plastic", measures=alpha_meas)
p + geom_boxplot(data=p$data, aes(x=Plastic, y=value),color ="black", fill="darkblue", outlier.fill = "red", alpha=0.1) + theme_bw() + ggtitle("Alpha Diversity") + theme(text = element_text(size = 25)) + geom_point(aes(colour=Plastic), size=4) + geom_point(shape = 21,size = 4,colour = "black")
dev.off()

alpha_diversity_graph(ps, index = 'shannon', treatment = c('Plastic'), subset = NULL, colors = 'default')

# Plot da diversidade beta

ordinate(ps, "PCoA", "bray") %>%
plot_ordination(ps, ., color = "Plastic", title = "Bray-Curtis") + geom_point(aes(size = 10)) + theme_bw() + ggtitle("Beta Diversity Bray-Curtis - Considered all ASVs") + theme(text = element_text(size = 20))
dev.off()

# Grafico de barras da abundancia

taxa_abundance_bars(FSfr, classification = 'Phylum', treatment = 'Plastic')

# Heatmap de abundancia

abundance_heatmap(FSfr, classification = 'Phylum', treatment = c('Plastic'), transformation = 'log10')

# Em linhas nao esta legal, mas o heatmap é mais explicativo

abundance_lines(FSfr, classification = 'Phylum', treatment = c('Plastic'), relative_abundance = TRUE)


nmds_phyloseq(FSfr, c('Sample', 'Plastic'), circle = TRUE, verbose = FALSE)

# PCOA

pcoa_phyloseq(FSfr, c('Sample','Plastic'), circle = TRUE)

tsne_phyloseq(FSfr, treatment = c('Sample','Plastic'), perplexity = 2)

# Dendograma

dendrogram_phyloseq(FSfr, treatment = 'Plastic', method = 'bray', colors = 'default')

# Plot taxa core

taxa_core_graph(FSfr, abundance_thresholds = seq(0.01, 0.25, 0.01))


# Rede de co ocorrência

filtered_obj <- conglomerate_taxa(ps, "Phylum")
co_occurrence_network(filtered_obj, treatment = 'Plastic', 
                      classification = 'Phylum')

# Rede de correlacao

filtered_obj <- taxa_filter(ps, frequency = 0.65)
variable_correlation_network(filtered_obj, variables = 'Sample', classification = 'Phylum', method = 'spearman')

# Arvore filogenetica

random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
random_tree_fs = rtree(ntaxa(FSfr), rooted=TRUE, tip.label=taxa_names(FSfr))

plotTree(random_tree_fs,ftype="i",fsize=0.6,lwd=1)
plotTree(random_tree_fs,type="fan",fsize=0.1,lwd=1, ftype="i")


ggplot(random_tree, branch.length="none", aes(x, y, label = NULL)) + geom_tree() + theme_tree()
ggplot(random_tree, aes(x, y, label = NULL)) + geom_tree() + theme_tree() + layout_inward_circular(xlim=15)


ggtree(random_tree) + geom_tiplab() + theme_tree()

ggtree(random_tree, layout = "circular") + geom_tiplab(size=1, aes(angle=angle)) + theme_tree()
ggtree(random_tree, layout = "ellipse") + geom_tiplab(size=1, aes(angle=angle)) + theme_tree()

# Exportação para fasta

tax_to_fasta <- read_excel("./16S_taxa_data.xlsx")
dna <- Biostrings::DNAStringSet(tax_to_fasta$ASV)
Biostrings::writeXStringSet(dna,'asv.fasta')

