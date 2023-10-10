# O presente script analisa os dados de sequenciamento por amplicon na regiao V4 do 16S de rRNA de um gene

# Importacao da biblioteca do dada2

library(dada2); packageVersion("dada2")
library(tibble)
library(phyloseq)


# Importacao dos arquivos de saida do sequeciamento

path <- "./4 Semestre/IP III/Dados Projeto/dados/16S"
list.files(path)


# Inserir em "pattern" a forma com que os dados estao configurados

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


# Extracao dos nomes das amostras supondo que o formato seja "SAMPLENAME_X.fastq"
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# Plot da qualidade dos reads

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


# Adaptador 16S
FWD <- "NNNGTGYCAGCMGCCGCGGTAA"
REV <- "NNNGGACTACNVGGGTWTCTAAT"

allOrients <- function(primer)  {
  require(Biostrings)
  dna <- DNAString(primer) 
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  
} 

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients


# Filtragem e corte dos reads

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# Configurar a multithread como falsa caso seja rodado no Windows

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)


# Aprendizado dos erros das amostras
# p.e, A2C caso a leitura da base A tenha sido C

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


# Inferencia da amostra e aplicacao do algoritmo dada

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


# Mescla dos pares Forward e Reverse

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)


# Tabela de sequencia

seqtab <- makeSequenceTable(mergers)


# Remocao de chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))


# Acompanhamento de reads

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# Atribuicao de taxonomia

setwd("./4 Semestre/IP III/Dados Projeto")
taxa <- dada2::assignTaxonomy(seqtab.nochim, "./tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

# Atribuicao de taxonomia exata
# Atribuicao de taxonomia exata

taxa_data <- dada2::addSpecies(taxa, "./tax/silva_species_assignment_v138.1.fa.gz")

# Analisando as atribuicoes de taxonomia

taxa.print <- taxa # Removendo nome de sequências nas linhas para display
rownames(taxa.print) <- NULL
head(taxa.print)

taxa_data.print <- taxa_data
rownames(taxa_data.print) <- NULL
head(taxa_data.print)

df <- as.data.frame(t(seqtab.nochim))
seqtab_data <- tibble::rownames_to_column(df, "ASV")

library(writexl)

library(tibble)
taxa_data <- rownames_to_column(as.data.frame(taxa_data), var = "ASV")

write_xlsx(seqtab_data,'16S_seqtab.xlsx')
write_xlsx(as.data.frame(taxa_data), '16S_taxa_data.xlsx')
