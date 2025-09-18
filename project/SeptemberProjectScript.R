# Step 1: Gene Annotation

library(GenomicRanges)
library(GenomeInfoDb)
library(S4Vectors)
library(IRanges)
library(rtracklayer)

# Creazione cartella data se non esiste
if(!dir.exists("data")) dir.create("data")

# Scarico GTF da Ensembl Release 114
gtf_file <- "data/Homo_sapiens.GRCh38.114.gtf.gz"
if(!file.exists(gtf_file)) {
  message("⬇️ Scarico GTF da Ensembl...")
  download.file(
    url = "https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz",
    destfile = gtf_file
  )
  message("✅ Download completato.")
}

# Leggo GTF in GRanges
gtf <- import(gtf_file)
genes <- gtf[gtf$type == "gene"]

# Mantengo solo geni con gene_id valido
genes <- genes[!is.na(mcols(genes)$gene_id) & mcols(genes)$gene_id != ""]

# Rinomino i geni coding come ensemblID:geneSymbol
gene_biotype <- mcols(genes)$gene_biotype
gene_id     <- mcols(genes)$gene_id
gene_name   <- mcols(genes)$gene_name

# Aggiorno la colonna gene_id in GRanges
mcols(genes)$gene_id_updated <- ifelse(
  gene_biotype == "protein_coding",
  paste0(gene_id, ":", gene_name),  # coding: ensemblID:geneSymbol
  gene_id                           # noncoding: solo ensemblID
)

# Sostituisco la colonna gene_id con quella aggiornata
mcols(genes)$gene_id <- mcols(genes)$gene_id_updated
mcols(genes)$gene_id_updated <- NULL

# Salvo GRanges completo
gr_all <- genes
save(gr_all, file = "data/gr_all.RData")

# Divido in coding e noncoding
gr_coding    <- gr_all[gene_biotype == "protein_coding"]
gr_noncoding <- gr_all[gene_biotype != "protein_coding"]

# Stampa info
cat("✅ Totale geni:", length(gr_all), "\n")
cat("🧬 Geni coding:", length(gr_coding), "\n")
cat("🌀 Geni noncoding:", length(gr_noncoding), "\n")

# Salvo i subset
save(gr_coding, file = "data/gr_coding.RData")
save(gr_noncoding, file = "data/gr_noncoding.RData")
cat("💾 File salvati in data/: gr_all.RData, gr_coding.RData, gr_noncoding.RData\n")



# Step 2: PCA

# Carica librerie
library(dplyr)
library(ggplot2)
library(Biobase)

# Percorsi ai file TSV nella sottocartella data
file1 <- "data/GSE244486_raw_counts.tsv"
file2 <- "data/GSE244485_raw_counts.tsv"

# Leggi i dataset
expr1 <- read.table(file1, header = TRUE, sep = "\t", row.names = 1)
expr2 <- read.table(file2, header = TRUE, sep = "\t", row.names = 1)

cat("Dimensioni expr1:", dim(expr1), "\n")
cat("Dimensioni expr2:", dim(expr2), "\n")

# Trova geni comuni
common_genes <- intersect(rownames(expr1), rownames(expr2))
cat("Geni comuni:", length(common_genes), "\n")

if(length(common_genes) == 0){
  stop("❌ Nessun gene in comune tra i due dataset. Controlla piattaforma/annotazione.")
}

# Mantieni solo geni comuni
expr1_common <- expr1[common_genes, ]
expr2_common <- expr2[common_genes, ]

# Combina i dataset
combined_expr <- cbind(expr1_common, expr2_common)

# Trasponi per PCA (campioni come righe)
combined_expr_t <- t(combined_expr)

# PCA
pca_res <- prcomp(combined_expr_t, scale. = TRUE)

# Crea dataframe per ggplot
pca_df <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Dataset = rep(c("GSE244486", "GSE244485"),
                times = c(ncol(expr1_common), ncol(expr2_common)))
)

# Grafico PCA
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Dataset)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA dei dataset GEO", x = "PC1", y = "PC2")

ggsave("results/PCA_plot.png", plot = p, width = 6, height = 5)  # salva su file PNG

# Calcola i centroidi dei due dataset
centroids <- pca_df %>%
  group_by(Dataset) %>%
  summarise(PC1_mean = mean(PC1), PC2_mean = mean(PC2))

# Distanza euclidea tra i centroidi
dist_centroids <- sqrt((centroids$PC1_mean[1] - centroids$PC1_mean[2])^2 +
                         (centroids$PC2_mean[1] - centroids$PC2_mean[2])^2)

cat("Distanza tra i centroidi dei dataset (batch effect proxy):", dist_centroids, "\n")

if(dist_centroids > 5){
  cat("⚠️ Probabile batch effect rilevante.\n")
} else {
  cat("✅ Batch effect minimo.\n")
}


# Step 3.1: If a batch effect is detected (prepare metadata)

library(dplyr)

# File di input
file1 <- "data/GSE244486_raw_counts.tsv"
file2 <- "data/GSE244485_raw_counts.tsv"

# Leggi i due dataset
expr1 <- read.table(file1, header = TRUE, sep = "\t", row.names = 1)
expr2 <- read.table(file2, header = TRUE, sep = "\t", row.names = 1)

# Rinomina colonne con suffisso batch
colnames(expr1) <- paste0(colnames(expr1), "_486")
colnames(expr2) <- paste0(colnames(expr2), "_485")

# Salva versioni rinominate
write.table(expr1, "data/GSE244486_raw_counts_renamed.tsv",
            sep = "\t", quote = FALSE, col.names = NA)
write.table(expr2, "data/GSE244485_raw_counts_renamed.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# Costruisci metadata
samples <- c(colnames(expr1), colnames(expr2))
conditions <- gsub("_.*", "", samples)  # prende la parte prima di "_"

metadata <- data.frame(
  sample = samples,
  condition = conditions
)

# Salva metadata
write.table(metadata, "data/sample_metadata.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ Creati:\n - data/GSE244486_raw_counts_renamed.tsv\n - data/GSE244485_raw_counts_renamed.tsv\n - data/sample_metadata.tsv\n")

# Step 3.1: If a batch effect is detected 

library(data.table)
library(dplyr)
library(sva)
library(DESeq2)

# ---  1: Caricamento counts ---
counts_file_486 <- "data/GSE244486_raw_counts_renamed.tsv"
counts_file_485 <- "data/GSE244485_raw_counts_renamed.tsv"

counts_486 <- fread(counts_file_486, data.table = FALSE)
rownames(counts_486) <- counts_486[[1]]
counts_486 <- counts_486[,-1]

counts_485 <- fread(counts_file_485, data.table = FALSE)
rownames(counts_485) <- counts_485[[1]]
counts_485 <- counts_485[,-1]

counts_combined <- cbind(counts_486, counts_485)
colnames(counts_combined) <- make.unique(colnames(counts_combined))
cat("Dimensioni counts combinati:", dim(counts_combined), "\n")

# ---  1b: Caricamento metadata ---
metadata_file <- "data/sample_metadata.tsv"
metadata <- fread(metadata_file, data.table = FALSE)
rownames(metadata) <- metadata$sample
cat("Dimensioni metadata:", dim(metadata), "\n")

# Mantieni solo campioni in comune
common_samples <- intersect(colnames(counts_combined), rownames(metadata))
counts_combined <- counts_combined[, common_samples]
metadata <- metadata[common_samples, , drop=FALSE]
cat("Campioni in comune:", length(common_samples), "\n")

# --- 2: Correzione batch con ComBat-seq ---
metadata$batch <- ifelse(grepl("_486", rownames(metadata)), "486", "485")
combat_counts <- ComBat_seq(as.matrix(counts_combined), batch=metadata$batch, group=metadata$condition)
cat("Dimensioni combat_counts:", dim(combat_counts), "\n")

# --- 3: DESeq2 su dati COMBAT-corrected ---
dds_combat <- DESeqDataSetFromMatrix(countData = round(combat_counts),
                                     colData = metadata,
                                     design = ~ condition)
dds_combat <- DESeq(dds_combat)

res_infected_combat <- results(dds_combat, contrast=c("condition","infected","mock"))
res_bystander_combat <- results(dds_combat, contrast=c("condition","bystander","mock"))

write.csv(as.data.frame(res_infected_combat), file="results/results_DE_infected_vs_mock_COMBAT.csv", row.names=TRUE)
write.csv(as.data.frame(res_bystander_combat), file="results/results_DE_bystander_vs_mock_COMBAT.csv", row.names=TRUE)
cat("DE results COMBAT salvati in CSV\n")

# --- 4: DESeq2 su dati RAW includendo batch come covariata ---
dds_raw <- DESeqDataSetFromMatrix(countData = round(counts_combined),
                                  colData = metadata,
                                  design = ~ batch + condition)
dds_raw <- DESeq(dds_raw)

res_infected_raw <- results(dds_raw, contrast=c("condition","infected","mock"))
res_bystander_raw <- results(dds_raw, contrast=c("condition","bystander","mock"))

write.csv(as.data.frame(res_infected_raw), file="results/results_DE_infected_vs_mock_RAW.csv", row.names=TRUE)
write.csv(as.data.frame(res_bystander_raw), file="results/results_DE_bystander_vs_mock_RAW.csv", row.names=TRUE)
cat("DE results RAW+batch salvati in CSV\n")

# --- 5: Confronto COMBAT vs RAW+batch ---
# Trasformiamo in data.frame e aggiungiamo il nome dei geni
df_combat_infected <- as.data.frame(res_infected_combat)
df_combat_infected$gene <- rownames(df_combat_infected)

df_raw_infected <- as.data.frame(res_infected_raw)
df_raw_infected$gene <- rownames(df_raw_infected)

compare_infected <- merge(df_combat_infected[, c("gene", "log2FoldChange", "padj")],
                          df_raw_infected[, c("gene", "log2FoldChange", "padj")],
                          by="gene",
                          suffixes=c("_combat", "_raw"))

write.csv(compare_infected, file="results/compare_DE_infected_vs_mock.csv", row.names=FALSE)
cat("Confronto DE infected vs mock COMBAT vs RAW salvato\n")

df_combat_bystander <- as.data.frame(res_bystander_combat)
df_combat_bystander$gene <- rownames(df_combat_bystander)

df_raw_bystander <- as.data.frame(res_bystander_raw)
df_raw_bystander$gene <- rownames(df_raw_bystander)

compare_bystander <- merge(df_combat_bystander[, c("gene", "log2FoldChange", "padj")],
                           df_raw_bystander[, c("gene", "log2FoldChange", "padj")],
                           by="gene",
                           suffixes=c("_combat", "_raw"))

write.csv(compare_bystander, file="results/compare_DE_bystander_vs_mock.csv", row.names=FALSE)
cat("Confronto DE bystander vs mock COMBAT vs RAW salvato\n")
# --- Salvataggio oggetti DESeq2 ---
save(dds_combat, file="results/dds_combat.RData")
save(dds_raw, file="results/dds_raw.RData")
cat("Oggetti DESeq2 salvati: dds_combat e dds_raw\n")


# Step 4: Heat Map and Clustering

library(DESeq2)
library(pheatmap)
library(matrixStats)
library(grid)
library(RColorBrewer)

# --- Carica dds_combat ---
load("results/dds_combat.RData")  # contiene dds_combat

# --- Variance stabilizing transformation ---
vsd <- vst(dds_combat, blind = FALSE)

# --- Selezione top 100 geni più variabili ---
rv <- matrixStats::rowVars(assay(vsd), useNames = TRUE)
select <- order(rv, decreasing = TRUE)[1:100]

mat <- assay(vsd)[select, ]
mat <- mat - rowMeans(mat)

# --- Metadata per annotazioni ---
anno <- as.data.frame(colData(dds_combat)[, c("condition","batch")])

# --- Palette colori più leggibile ---
col_fun <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

# --- Crea heatmap con maggiore spazio per dendrogrammi ---
hm <- pheatmap(mat,
               annotation_col = anno,
               show_rownames = FALSE,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               color = col_fun,
               main = "Top 100 variable genes - COMBAT corrected",
               silent = TRUE,
               treeheight_row = 50,   # altezza dendrogramma righe
               treeheight_col = 50,   # altezza dendrogramma colonne
               fontsize_row = 6,      # dimensione testo righe
               fontsize_col = 8       # dimensione testo colonne
)

# --- Salva su PDF e PNG ---
pdf("results/heatmap_top100.pdf", width=12, height=10)  # ingrandito
grid::grid.newpage()
grid::grid.draw(hm$gtable)
dev.off()

png("results/heatmap_top100.png", width=1500, height=1200, res=150)
grid::grid.newpage()
grid::grid.draw(hm$gtable)
dev.off()

cat("✅ Heatmap salvata in results/heatmap_top100.pdf\n")


# Step 5: Coding vs Noncoding Gene Relationships

library(DESeq2)
library(GenomicRanges)
library(BiocGenerics)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)

# ----- 1️⃣ Caricamento DE results da CSV -----
cat("✅ Caricamento DE results...\n")
de_infected <- read.csv("results/results_DE_infected_vs_mock_COMBAT.csv", row.names = 1)
de_bystander <- read.csv("results/results_DE_bystander_vs_mock_COMBAT.csv", row.names = 1)
cat("✅ DE results caricati correttamente\n")

# ----- 2️⃣ Caricamento GenomicRanges -----
cat("✅ Caricamento GenomicRanges...\n")
load("data/gr_coding.RData")      # coding genes
load("data/gr_noncoding.RData")   # noncoding genes
cat("✅ GenomicRanges caricati correttamente\n")

# ----- 3️⃣ Caricamento oggetto dds_combat -----
cat("✅ Caricamento dds_combat...\n")
load("results/dds_combat.RData")  # contiene dds_combat
cat("✅ dds_combat caricato correttamente\n")

# ----- 4️⃣ VST transformation -----
cat("✅ Calcolo VST...\n")
vst_obj <- vst(dds_combat, blind = FALSE)
cat("✅ VST calcolata\n")

# ----- 5️⃣ Estrazione top 5 geni DE significativi -----
cat("✅ Pulizia gene IDs e selezione top 5...\n")
# Filtra DE significativi
de_genes_infected <- rownames(de_infected)[de_infected$padj < 0.05]
de_genes_bystander <- rownames(de_bystander)[de_bystander$padj < 0.05]

# Intersect con geni presenti in VST
de_genes_infected <- intersect(de_genes_infected, rownames(vst_obj))
de_genes_bystander <- intersect(de_genes_bystander, rownames(vst_obj))

# Ordina per log2FoldChange assoluto e prendi top 5
top5_infected <- head(de_genes_infected[order(abs(de_infected[de_genes_infected, "log2FoldChange"]), decreasing = TRUE)], 5)
top5_bystander <- head(de_genes_bystander[order(abs(de_bystander[de_genes_bystander, "log2FoldChange"]), decreasing = TRUE)], 5)
cat("✅ Top 5 geni selezionati\n")

# ----- 6️⃣ Funzione per violin plot -----
plot_violin <- function(genes, vst_obj, group_col = "condition", title = "Violin Plot") {
  mat <- assay(vst_obj)[genes, , drop = FALSE]
  df <- as.data.frame(t(mat))
  df$sample <- rownames(df)
  
  df_long <- pivot_longer(df, cols = -sample, names_to = "gene", values_to = "expression")
  df_long$group <- colData(vst_obj)[df_long$sample, group_col]
  
  p <- ggplot(df_long, aes(x = gene, y = expression, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, y = "VST expression", x = "Gene")
  
  return(p)
}

# ----- 7️⃣ Generazione dei plot per top 5 -----
cat("✅ Generazione violin plots top 5...\n")
p_infected <- plot_violin(top5_infected, vst_obj, group_col = "condition", title = "Top 5 Infected vs Mock")
p_bystander <- plot_violin(top5_bystander, vst_obj, group_col = "condition", title = "Top 5 Bystander vs Mock")

# Salvataggio
ggsave("results/violin_top5_infected_vs_mock.png", p_infected, width = 10, height = 6)
ggsave("results/violin_top5_bystander_vs_mock.png", p_bystander, width = 10, height = 6)
cat("✅ Violin plots top 5 salvati in results/\n")


