library(dplyr)

splicing <- read.table("../geo/GSE60630_PDC3_vs_DMSO_splicing.diff", header = T)
gene <- read.table("../geo/GSE60630_PDC3_vs_DMSO_gene_exp.diff", header = T)

up_genes <- gene$gene_id[gene$q_value < 0.01 & gene$log2.fold_change. > 1]
dn_genes <- gene$gene_id[gene$q_value < 0.01 & gene$log2.fold_change. < -1]

sig_splicing <- splicing %>%
  filter(q_value < 0.05)

sum(sig_splicing$gene_id %in% up_genes)
sum(sig_splicing$gene_id %in% dn_genes)

iso <- read.table("../geo/GSE60630_PDC3_vs_DMSO_isoform_exp.diff", header = T)
names(iso)[10] <- "fc"
names(iso)[1] <- "acc"

get_ex1_end <- function(exon_start, exon_end, strand) {
  if(strand == "+") {
    ex1_end <- as.numeric(strsplit(exon_end,",")[[1]][1])
  } else if(strand == "-") {
    ex_split <- strsplit(exon_start,",")
    ex1_end <- as.numeric(ex_split[[1]][length(ex_split[[1]])])
  }
  ex1_end
}

refgene <- read.table("../refseq/refGene.txt")
names(refgene) <- c("bin","acc","chr","strand","tx_start","tx_end","cds_start","cds_end","exon_count","exon_start","exon_end",
                    "score","gene_id","cds_start_stat","cds_end_stat","exon_frames")

refgene <- refgene %>%
  rowwise() %>%
  mutate(ex1_end = get_ex1_end(exon_start, exon_end, strand))

tss <- refgene %>%
  select(acc, tx_start)

iso_diff <- iso %>%
  left_join(tss) %>%
  filter(status == "OK") %>%
  group_by(gene_id, tx_start) %>%
  summarize(n_up = sum(q_value < 0.01 & fc > 1),
            n_dn = sum(q_value < 0.01 & fc < -1),
            n_nc = sum(q_value >= 0.01 & abs(fc) < 1),
            n_iso = n()) %>%
  rowwise() %>%
  mutate(n_de = sum(n_up, n_dn))

gt1_iso <- iso_diff %>%
  filter(n_iso > 1)

gt1_iso_diff <- gt1_iso %>%
  # at least 1 differentially expressed isoform
  filter(n_de > 0) %>%
  # either differential in two directions or one isoform changes
  # while others don't
  filter((n_up > 0 & n_dn > 0) | (n_de > 0 & n_nc > 0))

sum(gt1_iso_diff$n_up)
sum(gt1_iso_diff$n_dn)
sum(gt1_iso_diff$n_nc)

up_iso_diff <- gt1_iso_diff %>%
  filter(n_up > 0 & n_nc > 0)

dn_iso_diff <- gt1_iso_diff %>%
  filter(n_dn > 0 & n_nc > 0)

up_set <- iso %>%
  filter(status == "OK") %>%
  filter(gene_id %in% up_iso_diff$gene_id) %>%
  mutate(group = ifelse(q_value >= 0.01,"nc","sig"),
         group = ifelse(q_value < 0.01 & fc > 1, "up", group))

dn_set <- iso %>%
  filter(status == "OK") %>%
  filter(gene_id %in% dn_iso_diff$gene_id) %>%
  mutate(group = ifelse(q_value >= 0.01,"nc","sig"),
         group = ifelse(q_value < 0.01 & fc < -1, "dn", group))

write.table(up_set,"../cuffdiff/isoform_up.txt", row.names = F, quote = F)
write.table(dn_set,"../cuffdiff/isoform_dn.txt", row.names = F, quote = F)

diff_set <- iso %>%
  filter(status == "OK") %>%
  filter(gene_id %in% gt1_iso_diff$gene_id) %>%
  mutate(group = ifelse(q_value >= 0.01,"nc","sig"),
         group = ifelse(q_value < 0.01 & fc > 1, "up", group),
         group = ifelse(q_value < 0.01 & fc < -1, "dn", group))

length(unique(diff_set$gene_id[diff_set$group != "nc"]))

write.table(diff_set,"../cuffdiff/isoform_set.txt",row.names = F, quote = F)

