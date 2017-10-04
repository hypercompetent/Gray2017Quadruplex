pile_stranded <- function(group_bed, bw_plus_gr, bw_minus_gr, window_size = 500, norm = TRUE) {
  
  group_gr <- GRanges(seqnames = group_bed$chr,
                      IRanges(group_bed$start,
                              group_bed$end),
                      strand = group_bed$strand)
  
  group_plus <- group_gr[strand(group_gr) == "+"]
  group_minus <- group_gr[strand(group_gr) == "-"]
  
  print("Finding Overlaps")
  
  ol_p_nt <- suppressWarnings(as.data.frame(findOverlaps(group_plus, bw_plus_gr)))
  names(ol_p_nt) <- c("group_hit","g4_hit")
  
  ol_p_t <- suppressWarnings(as.data.frame(findOverlaps(group_plus, bw_minus_gr)))
  names(ol_p_t) <- c("group_hit","g4_hit")
  
  ol_m_nt <- suppressWarnings(as.data.frame(findOverlaps(group_minus, bw_minus_gr)))
  names(ol_m_nt) <- c("group_hit","g4_hit")
  
  ol_m_t <- suppressWarnings(as.data.frame(findOverlaps(group_minus, bw_plus_gr)))
  names(ol_m_t) <- c("group_hit","g4_hit")
  
  pile_nt <- data.frame(pos = -window_size:window_size, val = 0)
  pile_t <- data.frame(pos = -window_size:window_size, val = 0)
  
  window_end <- 2*window_size + 1
  
  print("Calculating NT pileup")
  
  targets <- group_plus
  ol <- ol_p_nt
  g4 <- bw_plus_gr
  group_start <- start(targets)[ol$group_hit]
  group_end <- end(targets)[ol$group_hit]
  group_strand <- as.character(strand(targets)[ol$group_hit])
  g4_start <- start(g4)[ol$g4_hit]
  g4_end <- end(g4)[ol$g4_hit]
  g4_width <- width(g4)[ol$g4_hit]
  
  for(i in 1:nrow(ol)) {
    
    if(group_strand[i] == "+") {
      hit_start <- g4_start[i] - group_start[i] + 1
    } else if (group_strand[i] == "-") {
      hit_start <- group_end[i] - g4_end[i] + 1
    }
    hit_end <- hit_start + g4_width[i]
    
    if(hit_start < 1) { hit_start <- 1 }
    if(hit_end > window_end) { hit_end <- window_end }
    
    if( hit_start <= window_end & hit_end >= 1) {
      pile_nt$val[hit_start:hit_end] <- pile_nt$val[hit_start:hit_end] + 1
    }
    
  }
  
  targets <- group_minus
  ol <- ol_m_nt
  g4 <- bw_minus_gr
  group_start <- start(targets)[ol$group_hit]
  group_end <- end(targets)[ol$group_hit]
  group_strand <- as.character(strand(targets)[ol$group_hit])
  g4_start <- start(g4)[ol$g4_hit]
  g4_end <- end(g4)[ol$g4_hit]
  g4_width <- width(g4)[ol$g4_hit]
  
  for(i in 1:nrow(ol)) {
    
    if(group_strand[i] == "+") {
      hit_start <- g4_start[i] - group_start[i] + 1
    } else if (group_strand[i] == "-") {
      hit_start <- group_end[i] - g4_end[i] + 1
    }
    hit_end <- hit_start + g4_width[i]
    
    if(hit_start < 1) { hit_start <- 1 }
    if(hit_end > window_end) { hit_end <- window_end }
    
    if( hit_start <= window_end & hit_end >= 1) {
      pile_nt$val[hit_start:hit_end] <- pile_nt$val[hit_start:hit_end] + 1
    }
    
  }
  
  if(norm) {
    pile_nt$val <- pile_nt$val/(length(group_plus) + length(group_minus))
  }
  
  
  print("Calculating TS pileup")
  
  targets <- group_plus
  ol <- ol_p_t
  g4 <- bw_minus_gr
  group_start <- start(targets)[ol$group_hit]
  group_end <- end(targets)[ol$group_hit]
  group_strand <- as.character(strand(targets)[ol$group_hit])
  g4_start <- start(g4)[ol$g4_hit]
  g4_end <- end(g4)[ol$g4_hit]
  g4_width <- width(g4)[ol$g4_hit]
  
  for(i in 1:nrow(ol)) {
    
    if(group_strand[i] == "+") {
      hit_start <- g4_start[i] - group_start[i] + 1
    } else if (group_strand[i] == "-") {
      hit_start <- group_end[i] - g4_end[i] + 1
    }
    hit_end <- hit_start + g4_width[i]
    
    if(hit_start < 1) { hit_start <- 1 }
    if(hit_end > window_end) { hit_end <- window_end }
    
    if( hit_start <= window_end & hit_end >= 1) {
      pile_t$val[hit_start:hit_end] <- pile_t$val[hit_start:hit_end] + 1
    }
    
  }
  
  targets <- group_minus
  ol <- ol_m_t
  g4 <- bw_plus_gr
  group_start <- start(targets)[ol$group_hit]
  group_end <- end(targets)[ol$group_hit]
  group_strand <- as.character(strand(targets)[ol$group_hit])
  g4_start <- start(g4)[ol$g4_hit]
  g4_end <- end(g4)[ol$g4_hit]
  g4_width <- width(g4)[ol$g4_hit]
  
  for(i in 1:nrow(ol)) {
    
    if(group_strand[i] == "+") {
      hit_start <- g4_start[i] - group_start[i] + 1
    } else if (group_strand[i] == "-") {
      hit_start <- group_end[i] - g4_end[i] + 1
    }
    hit_end <- hit_start + g4_width[i]
    
    if(hit_start < 1) { hit_start <- 1 }
    if(hit_end > window_end) { hit_end <- window_end }
    
    if( hit_start <= window_end & hit_end >= 1) {
      pile_t$val[hit_start:hit_end] <- pile_t$val[hit_start:hit_end] + 1
    }
    
  }
  
  if(norm) {
    pile_t$val <- pile_t$val/(length(group_plus) + length(group_minus))
  }
  
  list(pile_nt = pile_nt,
       pile_t = pile_t)
  
}


extract_results <- function(g4_statistics, region = "ex1", side = "dn", groups, strand = "nt", cutoff = 0.05) {
  
  results <- character()
  
  for(group in groups) {
    lt.results <- g4_statistics[[region]][[side]][[group]][[strand]]["fdr.lt","sum"]
    gt.results <- g4_statistics[[region]][[side]][[group]][[strand]]["fdr.gt","sum"]
    for(i in 1:length(lt.results)) {
      if(lt.results[i] < gt.results[i] & lt.results[i] < cutoff) {
        results <- c(results, paste0("Depleted (",lt.results[i],")"))
      } else if(gt.results[i] < lt.results[i] & gt.results[i] < cutoff) {
        results <- c(results, paste0("Enriched (",gt.results[i],")"))
      } else {
        results <- c(results, "ns")
      }
    }
  }
  
  return(results)
  
}


bigwig_overlap_stranded <- function(group_bed, bw_plus_gr, bw_minus_gr, out_file = NULL, write_out = FALSE) {
  
  group_gr <- GRanges(seqnames = group_bed$chr,
                      IRanges(group_bed$start,
                              group_bed$end),
                      strand = group_bed$strand)
  
  plus_ol <- countOverlaps(group_gr, bw_plus_gr, ignore.strand = TRUE)
  minus_ol <- countOverlaps(group_gr, bw_minus_gr, ignore.strand = TRUE)
  
  results <- group_bed %>%
    mutate(plus_ol = plus_ol,
           minus_ol = minus_ol) %>%
    rowwise() %>%
    mutate(t.overlap = ifelse(strand == "+", minus_ol, plus_ol),
           nt.overlap = ifelse(strand == "+", plus_ol, minus_ol),
           sum.overlap = minus_ol + plus_ol) %>%
    select(-plus_ol, -minus_ol)
  
  if(write_out) {
    write.table(results, out_file, row.names = F, quote = F)
  }
  
  return(results)
  
}

get_list_items <- function(list,i,flip) {
  items <- vector()
  if (flip == FALSE) {
    for(x in 1:length(list)) {
      sublist <- list[[x]]
      items <- c(items,as.numeric(sublist[i]))
    }
  } else if (flip == TRUE) {
    for(x in 1:length(list)) {
      sublist <- list[[x]]
      items <- c(items,as.numeric(rev(sublist)[i]))
    }
  }
  return(items)
}

sort_bed <- function(df) {
  sorted <- df[with(df,order(chr,start)),]
  return(sorted)
}

get_regions <- function(region,half_window,get_ref_by,ids) {
  #define the group region output list
  group_region <- list()
  
  #for each set of ids for each group
  for(i in 1:length(ids)) {
    
    #get just the ref entries that match the ids
    if(get_ref_by == "acc") {
      sub_ref <- ref[ref$acc %in% ids[[i]], ]
    } else if(get_ref_by == "name") {
      sub_ref <- ref[ref$name %in% ids[[i]], ]
    }
    

    #set up which entries in the ref table will be the center of
    #the windows
    if (region == "ex1") {
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exst,1,flip=F)
      minus_center <- get_list_items(exen,1,flip=T)
    
      } else if (region == "ex2") {
        
      sub_ref <- sub_ref[sub_ref$exons >= 2,]

      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exst,2,flip=F)
      minus_center <- get_list_items(exen,2,flip=T)
      
    } else if (region == "ex3") {
      
      sub_ref <- sub_ref[sub_ref$exons >= 3,]
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exst,3,flip=F)
      minus_center <- get_list_items(exen,3,flip=T)
      
    } else if (region == "int1") {
      
      sub_ref <- sub_ref[sub_ref$exons >= 2,]
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exen,1,flip=F)
      minus_center <- get_list_items(exst,1,flip=T)
      
    } else if (region == "int2") {
      
      sub_ref <- sub_ref[sub_ref$exons >= 3,]
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exen,2,flip=F)
      minus_center <- get_list_items(exst,2,flip=T)
      
    } else if (region == "int3") {
      
      sub_ref <- sub_ref[sub_ref$exons >= 4,]
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exen,3,flip=F)
      minus_center <- get_list_items(exst3,flip=T)
      
    } else if (region == "end") {
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exen,1,flip=T)
      minus_center <- get_list_items(exst,1,flip=F)
      
    }
    
    rgst <- numeric()
    rgen <- numeric()
    for(i in 1:nrow(sub_ref)) {
      strand <- sub_ref$strand[i]
      if(strand == "+") {
        rgst[i] <- plus_center[i] - half_window
        rgen[i] <- plus_center[i] + half_window - 1
      } else if (strand == "-") {
        rgst[i] <- minus_center[i] - half_window + 1
        rgen[i] <- minus_center[i] + half_window
      }
    }
    
    if(get_ref_by == "acc") {
      names <- sub_ref$acc
    } else if (get_ref_by == "name") {
      names <- sub_ref$name
    }
    
    regions <- unique(data.frame(chr=sub_ref$chr,start=rgst,end=rgen,name=names,score=0,strand=sub_ref$strand))
    regions <- sort_bed(regions)
    
    group_region <- c(group_region,list(regions))

  }
  names(group_region) <- group_names
  return(group_region)
  
}

get_half_regions <- function(region,half_window,get_ref_by,ids) {
  #define the group region output list
  group_region_up <- list()
  group_region_dn <- list()
  
  #for each set of ids for each group
  for(i in 1:length(ids)) {
    
    #get just the ref entries that match the ids
    if(get_ref_by == "acc") {
      sub_ref <- ref[ref$acc %in% ids[[i]], ]
    } else if(get_ref_by == "name") {
      sub_ref <- ref[ref$name %in% ids[[i]], ]
    }
    

    #set up which entries in the ref table will be the center of
    #the windows
    if (region == "ex1") {
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exst,1,flip=F)
      minus_center <- get_list_items(exen,1,flip=T)
    
      } else if (region == "ex2") {
        
      sub_ref <- sub_ref[sub_ref$exons >= 2,]

      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exst,2,flip=F)
      minus_center <- get_list_items(exen,2,flip=T)
      
    } else if (region == "ex3") {
      
      sub_ref <- sub_ref[sub_ref$exons >= 3,]
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exst,3,flip=F)
      minus_center <- get_list_items(exen,3,flip=T)
      
    } else if (region == "int1") {
      
      sub_ref <- sub_ref[sub_ref$exons >= 2,]
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exen,1,flip=F)
      minus_center <- get_list_items(exst,1,flip=T)
      
    } else if (region == "int2") {
      
      sub_ref <- sub_ref[sub_ref$exons >= 3,]
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exen,2,flip=F)
      minus_center <- get_list_items(exst,2,flip=T)
      
    } else if (region == "int3") {
      
      sub_ref <- sub_ref[sub_ref$exons >= 4,]
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exen,3,flip=F)
      minus_center <- get_list_items(exst3,flip=T)
      
    } else if (region == "end") {
      
      exst <- strsplit(sub_ref$exst,",")
      exen <- strsplit(sub_ref$exen,",")
      
      plus_center <- get_list_items(exen,1,flip=T)
      minus_center <- get_list_items(exst,1,flip=F)
      
    }
    
    uprgst <- numeric()
    uprgen <- numeric()

    dnrgst <- numeric()
    dnrgen <- numeric()

    for(i in 1:nrow(sub_ref)) {
      strand <- sub_ref$strand[i]
      if(strand == "+") {
        uprgst[i] <- plus_center[i] - half_window
        uprgen[i] <- plus_center[i]

        dnrgst[i] <- plus_center[i]
        dnrgen[i] <- plus_center[i] + half_window - 1
      } else if (strand == "-") {
        uprgst[i] <- minus_center[i]
        uprgen[i] <- minus_center[i] + half_window

        dnrgst[i] <- minus_center[i] - half_window + 1
        dnrgen[i] <- minus_center[i]
      }
    }
    
    if(get_ref_by == "acc") {
      names <- sub_ref$acc
    } else if (get_ref_by == "name") {
      names <- sub_ref$name
    }
    
    upregions <- unique(data.frame(chr=sub_ref$chr,start=uprgst,end=uprgen,name=names,score=0,strand=sub_ref$strand))
    upregions <- sort_bed(upregions)

    dnregions <- unique(data.frame(chr=sub_ref$chr,start=dnrgst,end=dnrgen,name=names,score=0,strand=sub_ref$strand))
    dnregions <- sort_bed(dnregions)
    
    group_region_up <- c(group_region_up,list(upregions))
    group_region_dn <- c(group_region_dn,list(dnregions))
    
  }
  names(group_region_up) <- group_names
  names(group_region_dn) <- group_names
  
  results <- c(list(group_region_up),list(group_region_dn))
  
  split_names <- c("up","dn")
  names(results) <- split_names
  
  return(results)
  
}

