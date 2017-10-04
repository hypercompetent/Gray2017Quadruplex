##########
#Settings#
##########


#output directory
outdir <- "exp_quintile_pdc3_plot/"

#gene sets
gene_set_file <- "../geo/GSE60630_PDC3_vs_DMSO_gene_exp.diff"
gene_sets_raw <- read.table(gene_set_file,sep="\t",stringsAsFactors=F,header = T)
gene_sets <- gene_sets_raw %>% filter(status == "OK")

vals <- gene_sets$value_2
vals <- vals[order(vals)]
quantiles <- c(vals[1],vals[1:5*length(vals)/5])

#how to interface with refseq file. can be: "acc" or "name"
get_ref_by <- "name"

#accession and gene name columns
if (get_ref_by == "acc") {
  acc_col <- 1
  names(gene_sets)[acc_col] <- "acc"
}
if (get_ref_by == "name") {
  name_col <- 1
  names(gene_sets)[name_col] <- "name"
}

#group definitions
groups <- list(gene_sets$value_2 >= quantiles[1] & gene_sets$value_2 <= quantiles[2],
               gene_sets$value_2 > quantiles[2] & gene_sets$value_2 <= quantiles[3],
               gene_sets$value_2 > quantiles[3] & gene_sets$value_2 <= quantiles[4],
               gene_sets$value_2 > quantiles[4] & gene_sets$value_2 <= quantiles[5],
               gene_sets$value_2 > quantiles[5] & gene_sets$value_2 <= quantiles[6])
group_colors <- rev(c("#FF0000","#FA931E","#F0C943","#008000","#3B8CFB"))
group_names <- rev(c("quint1","quint2","quint3","quint4","quint5"))

#refseq file
ref_file <- "../refseq/refGene.txt"
ref <- read.table(ref_file,sep="\t",stringsAsFactors=F)
names(ref)[2] <- "acc"
names(ref)[13] <- "name"
names(ref)[3] <- "chr"
names(ref)[4] <- "strand"
names(ref)[9] <- "exons"
names(ref)[10] <- "exst"
names(ref)[11] <- "exen"

#refseq filtering
# filter_file <- "allgenes.txt"
# filter <- read.table(filter_file,sep="\t",stringsAsFactors=F)
filter <- gene_sets
# names(filter)[5] <- "name"
# names(filter)[6] <- "acc"

if (get_ref_by == "acc") {
  ref <- ref[ref$acc %in% filter$acc,]
}
if (get_ref_by == "name") {
  ref <- ref[ref$name %in% filter$name,]
}


#remove annotations on non-standard chromosomes, such as chr1_gl000192_random?
remove_ns_chr <- TRUE
#remove annotations for highly repetitive genes (hrgs), such as UTY (77 locations)
remove_hrgs = FALSE
if(remove_hrgs) {
  hrg_cutoff = 10
}

#window sizes
full_window <- 500
half_window <- round(full_window/2,0)

#regions can be any of (ex1,int1,ex2,int2,ex3,int3,end)
regions <- c("ex1","int1")

#bigwig file
bw_file_plus <- "../bw/g4-12_plus.bw"
bw_file_minus <- "../bw/g4-12_minus.bw"

bw_plus_gr <- import(bw_file_plus, format = "bw")
bw_minus_gr <- import(bw_file_minus, format = "bw")


#random control settings
use_rc <- TRUE
n_random <- 2000
random_color <- "gray50"

##############################
#Shouldn't need to edit below#
##############################

############################
#Get regions for each group#
############################
source("functions.R")

# Make an output directory
dir.create(outdir)

if(remove_ns_chr) {
  ref <- ref[grep("_",ref$chr,invert=T),]
}

if(remove_hrgs) {
  gene_counts <- data.frame(table(ref[,get_ref_by]))
  bad_genes <- gene_counts$Var1[gene_counts$freq >= hrg_cutoff]
  ref <- ref[!ref[,get_ref_by] %in% bad_genes,]
}

# get the identifiers (either acc or name) for each group
# this yields a list of vectors containing the appropriate id columns
group_ids <- list()
for(i in 1:length(groups)) {
  group_ids <- c(group_ids,list(gene_sets[groups[[i]],get_ref_by]))
}

if(use_rc) {
  group_names <- c("random",group_names)
  group_ids <- c(list(sample(ref[,get_ref_by],n_random,replace=F)),group_ids)
  group_colors <- c(random_color,group_colors)
}

# use the ids to retrieve locations from the ref table, and use the half-window size to build
# window locations based on each region that we're looking at (ex1, int1, etc)
# this yields a list of lists per region containing a list of data frames per group.
# thus, one can retrieve a data frame for ex1 windows for a group called "group1" using group_regions$ex1$group1
group_regions <- list()
for(i in 1:length(regions)) {
  reg <- regions[i]
  group_regions <- c(group_regions,list(get_regions(reg,half_window,get_ref_by,group_ids)))
  names(group_regions)[[i]] <- reg
}

#################
#Perform pileups#
#################
require(ggplot2)

max_val <- 0
min_val <- 0

pileups <- list()

for(reg in regions) {
  region_pile <- data.frame(pos=-half_window:(half_window))
  
  for(group in group_names) {
    bedname <- paste(outdir,reg,"_",group,".bed",sep="")
    write.table(group_regions[[reg]][group],bedname,quote=F,row.names=F,col.names=F,sep="\t")
    
    pilename_nt <- paste(outdir,reg,"_",group,"_nt.pile",sep="")
    pilename_ts <- paste(outdir,reg,"_",group,"_ts.pile",sep="")
    
    group_bed <- group_regions[[reg]][[group]]
    
    group_piles <- pile_stranded(group_bed, bw_plus_gr, bw_minus_gr, window_size = 250, norm = TRUE)
    
    region_pile <- cbind(region_pile, group_piles$pile_nt$val)
    names(region_pile)[ncol(region_pile)] <- paste(group,"nt",sep="")
    region_pile <- cbind(region_pile, group_piles$pile_t$val*-1)
    names(region_pile)[ncol(region_pile)] <- paste(group,"ts",sep="")
    
  }
  reg_max_val <- max(region_pile[2:ncol(region_pile)])
  if(reg_max_val > max_val) { max_val <- reg_max_val }
  reg_min_val <- min(region_pile[2:ncol(region_pile)])
  if(reg_min_val < min_val) { min_val <- reg_min_val }
  
  pileups <- c(pileups,list(region_pile))
  names(pileups)[[length(pileups)]] <- reg
}

##############
#Plot pileups#
##############

max_val <- max_val*1.1
if(min_val > 0) { min_val <- 0 } else { min_val <- min_val*1.1 }

for(reg in regions) {
  region_pile <- pileups[[reg]]
  
  region_plot <- ggplot(data=region_pile) + 
    theme_bw() + 
    ylim(min_val,max_val) +
    scale_x_continuous(expand = c(0,0),limits = c(-100,250)) +
    geom_ribbon(aes_string(x="pos",ymin = "randomts",ymax = "randomnt"), color = group_colors[1], alpha = 0.2)
  for(i in 2:length(group_names)) {
    name_nt <- paste(group_names[i],"nt",sep="")
    name_ts <- paste(group_names[i],"ts",sep="")
    pos <- "pos"
    setcolor <- group_colors[i]
    region_plot = region_plot + geom_line(aes_string(x=pos,y=name_nt),color=setcolor)
    region_plot = region_plot + geom_line(aes_string(x=pos,y=name_ts),color=setcolor)
  }

  plotname <- paste(outdir,reg,"_plot.pdf",sep="")
  region_plot
  ggsave(plotname,width=6,height=6)
  
}

