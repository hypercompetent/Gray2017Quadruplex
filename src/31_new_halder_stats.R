library(dplyr)
library(rtracklayer)

##########
#Settings#
##########


#output directory
outdir <- "new_halder_g4_stats/"

#gene sets
gene_set_file <- "../treatments/halder_phendc3.txt"
gene_sets <- read.table(gene_set_file,sep="\t",stringsAsFactors=F,header = T)
gene_sets <- gene_sets[-1] %>% unique()

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
groups <- list(gene_sets[2] == "up",
               gene_sets[2] == "dn")
group_colors <- c("orangered","dodgerblue")
group_names <- c("pdc3up","pdc3dn")


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
filter_file <- "../treatments/all_halder_genes.txt"
filter <- read.table(filter_file,sep="\t",stringsAsFactors=F,header = T)
names(filter)[7] <- "name"

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
n_random <- 1000

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

group_names <- c("all",group_names)
group_ids <- c(list(ref[,get_ref_by]),group_ids)

# use the ids to retrieve locations from the ref table, and use the half-window size to build
# window locations based on each region that we're looking at (ex1, int1, etc)
# this yields a list of lists per region containing a list of data frames per group.
# thus, one can retrieve a data frame for ex1 windows for a group called "group1" using group_regions$ex1$group1
group_regions <- list()
for(i in 1:length(regions)) {
  reg <- regions[i]
  group_regions <- c(group_regions,list(get_half_regions(reg,half_window,get_ref_by,group_ids)))
  names(group_regions)[[i]] <- reg
}

#########################
#Calculate G4 Statistics#
#########################


g4_results <- group_regions

for(reg in regions) {
  for(half in c("up","dn")) {
    for(group in group_names) {
  
      bedname <- paste(outdir,reg,"_",half,"_",group,".bed",sep="")
      write.table(group_regions[[reg]][[half]][group],bedname,quote=F,row.names=F,col.names=F,sep="\t")
      
      resultname <- paste(outdir,reg,"_",half,"_",group,".txt",sep="")
      
      group_bed <- group_regions[[reg]][[half]][[group]]
      
      g4_table <- bigwig_overlap_stranded(group_bed, bw_plus_gr, bw_minus_gr, out_file = resultname, write_out = TRUE)
      
      g4_results[[reg]][[half]][[group]] <- g4_table[,c("t.overlap","nt.overlap","sum.overlap")]
      
    }
  }
}


#Calculate G4 statistics
dir.create(paste0(outdir,"t"))
dir.create(paste0(outdir,"nt"))
dir.create(paste0(outdir,"sum"))

exp_groups <- group_names[2:length(group_names)]
g4_statistics <- g4_results
for(reg in regions) {
  for(half in c("up","dn")) {
    g4_statistics[[reg]][[half]]$all <- NULL
    for(group in exp_groups) {
      
      all.results <- g4_results[[reg]][[half]]$all
      exp.results <- g4_results[[reg]][[half]][[group]]
      
      group_result <- list()
      
      for(i in 1:ncol(exp.results)) {
        res <- exp.results[[i]]
        strand <- sub("\\..+","",names(exp.results)[i])
        exp <- data.frame(sum=sum(res),mean=mean(res),nz=sum(res != 0),nz.perc=sum(res != 0)/length(res),nz.mean=mean(res[res != 0]))
        
        rand <- data.frame(sum=0,mean=0,nz=0,nz.perc=0,nz.mean=0)
        for(j in 1:n_random) {
          randset <- sample(all.results[[i]],length(res),replace=F)
          rand[j,] <- c(sum(randset),mean(randset),sum(randset != 0),sum(randset != 0)/length(randset),mean(randset[randset != 0]))
        }
        
        rand.mean <- c(mean(rand$sum),mean(rand$mean),mean(rand$nz),mean(rand$nz.perc),mean(rand$nz.mean))
        rand.sd <- c(sd(rand$sum),sd(rand$mean),sd(rand$nz),sd(rand$nz.perc),sd(rand$nz.mean))
        fdr.gt <- c(sum(rand$sum >= exp$sum)/n_random,
                    sum(rand$mean >= exp$mean)/n_random,
                    sum(rand$nz >= exp$nz)/n_random,
                    sum(rand$nz.perc >= exp$nz.perc)/n_random,
                    sum(rand$nz.mean >= exp$nz.mean)/n_random)
        fdr.lt <- c(sum(rand$sum <= exp$sum)/n_random,
                    sum(rand$mean <= exp$mean)/n_random,
                    sum(rand$nz <= exp$nz)/n_random,
                    sum(rand$nz.perc <= exp$nz.perc)/n_random,
                    sum(rand$nz.mean <= exp$nz.mean)/n_random)
        
        result <- rbind(exp,rand.mean,rand.sd,fdr.gt,fdr.lt)
        row.names(result) <- c(group,"rand.mean","rand.sd","fdr.gt","fdr.lt")
        
        group_result <- c(group_result,list(result))
        
        statfile <- paste(outdir,strand,"/",reg,"_",half,"_",group,"_stats.txt",sep="")
        write.table(result,statfile,quote=F,sep="\t")
        
        
      }
      
      names(group_result) <- c("t","nt","sum")
      
      g4_statistics[[reg]][[half]][[group]] <- group_result
      
    }
  }
}

g4_combined <- data.frame()

blank_line <- rep("",13)

strands <- c("nt","t","sum")

for(strand in strands) {
  for(group in exp_groups){
    for(region in regions) {
      regionline <- c(strand,region,group,rep("",10))
      header <- c("up","sum","mean","nz","nz.perc","nz.mean","",
                  "dn","sum","mean","nz","nz.perc","nz.mean")
      
      up <- g4_statistics[[region]]$up[[group]][[strand]]
      up <- round(up,3)
      up <- cbind(rownames(up),up)
      up <- data.frame(lapply(up, as.character), stringsAsFactors=FALSE)
      
      dn <- g4_statistics[[region]]$dn[[group]][[strand]]
      dn <- round(dn,3)
      dn <- cbind(rownames(dn),dn)
      dn <- data.frame(lapply(dn, as.character), stringsAsFactors=FALSE)
      
      comb <- cbind(up,"",dn)
      names(comb) <- names(g4_combined)
      
      g4_combined <- rbind(g4_combined,regionline,header,comb,blank_line)
      
    }
  }
}

comb_file <- paste(outdir,"g4_combined.txt",sep="/")
write.table(g4_combined,comb_file,quote=F,row.names=F,col.names=F,sep="\t")

library(xlsx)

xlsx_file <- paste(outdir,"g4_combined.xlsx",sep="/")
write.xlsx(g4_combined,xlsx_file,col.names=F,row.names=F)

simpler_g4_statistics <- g4_combined[,c(1,2,4,7,8,9,11)]
write.table(simpler_g4_statistics,paste(outdir,"g4_sum_nz.txt",sep="/"),quote=F,row.names=F,col.names=F,sep="\t")
write.xlsx(simpler_g4_statistics,paste(outdir,"g4_sum_nz.xls",sep="/"),row.names=F,col.names=F)



dead_simple <- data.frame(group = rep(exp_groups,2),
                          strand = rep(c("nt","t"),each = length(exp_groups)),
                          tss_250_up = c(extract_results(g4_statistics, "ex1","up",exp_groups,"nt"),
                                      extract_results(g4_statistics, "ex1","up",exp_groups,"t")),
                          tss_250_dn = c(extract_results(g4_statistics, "ex1","dn",exp_groups,"nt"),
                                      extract_results(g4_statistics, "ex1","dn",exp_groups,"t")),
                          int_250 = c(extract_results(g4_statistics, "int1","dn",exp_groups,"nt"),
                                      extract_results(g4_statistics, "int1","dn",exp_groups,"t"))) %>%
  arrange(strand, group)
write.xlsx(dead_simple, paste(outdir,"g4_simple_summary.xlsx",sep = "/"), row.names = F)
