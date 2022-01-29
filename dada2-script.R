library(dada2); packageVersion("dada2")
path <- "rawseq" # directory containing the fastq files.
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#序列质量
plotQualityProfile(fnFs[1:2])# quality profiles of the forward reads
plotQualityProfile(fnRs[1:2])# quality profile of the reverse reads

#序列过滤和裁剪
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

#计算错误率
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#基于错误模型进一步质控Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]#1指代第一个样品，样品1的数据中得到52个sequence variants

#序列拼接
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#生成asv表格
seqtab <- makeSequenceTable(mergers)
dim(seqtab) #生成5504ASV
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#去除嵌合体Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)#结果可知嵌合体占了17%

#统计上述分析步骤
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#序列物种注释Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

#另存物种注释变量，去除序列名，只显示物种信息
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#文件导出
setwd(path)
seqtable.taxa.plus <- cbind('#seq'=rownames(taxa), t(seqtab.nochim), taxa)
# ASV表格导出
write.csv(seqtab.nochim, "dada2_counts.csv", sep="\t", quote=F, row.names = T)
# 带注释文件的ASV表格导出
write.csv(seqtable.taxa.plus , "dada2_counts.taxon.species.csv", quote=F, row.names=F)
# track文件保存
write.table(track , "dada2_track.txt", sep="\t", quote=F, row.names = F)
