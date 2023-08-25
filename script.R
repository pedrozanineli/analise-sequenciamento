# O presente script analisa os dados de sequenciamento por amplicon na região V4 do 16S de rRNA de um gene

# Importação da biblioteca do dada2

library(dada2); packageVersion("dada2")

# Importação dos arquivos de saída do sequeciamento

path <- "~/"
list.files(path)

# Inserir em "pattern" a forma com que os dados estão configurados

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extração dos nomes das amostras supondo que o formato seja "SAMPLENAME_X.fastq"
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Plot da qualidade dos reads

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

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

# Inferência da amostra e aplicação do algoritmo dada

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Mescla dos pares Forward e Reverse

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Tabela de sequência

seqtab <- makeSequenceTable(mergers)

# Remoção de chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
