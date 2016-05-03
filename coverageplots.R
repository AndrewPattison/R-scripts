library(ShortRead)
library(rtracklayer)
library(IRanges)
library(ggplot2)
library(plotly)
library(dplyr)
library(ggbio)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(Rsamtools)
library(GenomicFeatures)
library(gridExtra)
#import("/data/home/apattison/Bioinformatics/work_outsde_my_project/dp_ipsc/genome_locations/Mus_musculus.GRCm38.84.gtf")
#ens_gff <- import("/data/home/apattison/Bioinformatics/work_outsde_my_project/dp_ipsc/genome_locations/Mus_musculus.GRCm38.84.gtf")
#ens_gff  <- makeTxDbFromBiomart(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")



setwd("/data/home/apattison/Bioinformatics/work_outsde_my_project/work_for_minni_anko/coverage_plots_r_code/")

# Get the to 20 genes from this gff file in terms of fold change
filter_gff_for_rows  <- function(gff){
  indexes <- data.frame(index = numeric(), fc = numeric())
  for (i in 1:nrow(gff)){
    sub_1 <- strsplit(gff[i,9],split = ";")[[1]][[3]]
    sub_2 <- gsub(pattern = "Change_in_distal_usage=", replacement =  "",x= sub_1)
    ind_fc <- data.frame(i,as.numeric(sub_2), stringsAsFactors = F)
    indexes <- rbind(indexes,ind_fc)
    
  }
  ind <- indexes
  ordered <- ind[order(abs(ind$as.numeric.sub_2.), decreasing = T),]
  indexes <- ordered[1:20,1]
  top_hits <- gff[indexes,]
  names <- character()
  for (i in 1:nrow(top_hits)){
    sub_1 <- strsplit(top_hits[i,9],split = ";")[[1]][[6]]
    sub_2 <- gsub(pattern = "parent=", replacement =  "",x= sub_1)
    sub_3 <- strsplit(top_hits[i,9],split = ";")[[1]][[5]]
    sub_4 <- gsub(pattern = "name=", replacement =  "",x= sub_3)
    sub_5 <- strsplit(sub_4, split = "_")[[1]][[1]]
    sub_6 <- strsplit(sub_2, split = "_")[[1]][[1]]
    names <- c(names, paste0(sub_5, " (",sub_6, ")"))
  }
  top_hits$name <- names
  return(top_hits)
}

grab_regions <- function (gff_file){
  
  gff <- read.delim(gff_file, header=F, comment.char = "", stringsAsFactors=F)
  regions <- filter_gff_for_rows(gff)
  return(regions)
  
}
get_wig_output <- function (gff_regions,wig_file){  
  final_df <- data.frame(number_list = numeric(), Freq = numeric(), name = character())
  track_df <-  data.frame(start = numeric(), end = numeric(), strand = character(), name = character())
  for (line in seq(1, nrow(gff_regions),2)){
    chromosome <-gsub ("chr", "", gff_regions[line,1])
    strand <- gff_regions[line, 7]
    if (strand == "+"){
      if (grepl("_for", wig_file)){
        ori <- F
        start <- gff_regions[line,4] - 3000
        end <- gff_regions [line +1 ,5]  + 3000          
      }
      else{
        next
      }   
    }
    else{
      if (grepl("_rev", wig_file)){
        ori <- T
        start <- gff_regions[line +1 ,4] - 3000
        end <- gff_regions[line ,5] + 3000
      }
      else{
        next
      }
      
    }
    
    
    cat(chromosome, ori, start,end)
    bw <- BigWigFile(wig_file)
    suma <- summary(bw, GRanges (paste0("chr", chromosome),IRanges(start, end)), type ="max", size = end-start, defaultValue =0)
    
    cov_frame <- data.frame (suma)
    pos_score <- cov_frame[, c("start","score")]   
    colnames(pos_score) <- c("number_list", "Freq")
    
    which <- GRanges (paste0(chromosome),strand = strand ,IRanges(start, end))
    exon_features <- ens_gff[mcols(ens_gff)$type == "exon"]
    seqs <- subsetByOverlaps(exon_features,which)
    starts <- start(seqs)
    ends <- end(seqs)
    strands <- strand (seqs)
    
    track_df_part <- data.frame(starts, ends, strands)    
    
    name <- gff_regions[line,10]
    pos_score$name <- paste0 (name) #min(plot_df$number_list), " to ", max(plot_df$number_list))
    track_df_part$name <- paste0 (name)
    final_df <- rbind (final_df, pos_score)  
    track_df <- rbind(track_df, track_df_part)
  }
  return(list(final_df, track_df))
}

run_program <- function(gff_regions, wig_list, name_list){
  
  wig_files <- wig_list
  final_frame <- data.frame(number_list = numeric(), Freq = numeric(), name = character(), wig_file = character())
  tracks_frame <- data.frame(start = numeric(), end = numeric(), strand = character(), name = character(), wig_file = character())
  
  for (i in 1:length(wig_files)){
    print (wig_files[[i]])
    split <- strsplit(wig_files[[i]], split = "/")
    wig_name <- name_list[[i]] 
    
    wig_file_out <- get_wig_output(gff_regions, wig_files[[i]])
    if (nrow(wig_file_out[[1]]) != 0){
      wig_file_out[[1]]$wig_file <- wig_name
      final_frame <- rbind(final_frame, wig_file_out[[1]]) 
      wig_file_out[[2]]$wig_file <- wig_name
      tracks_frame <- rbind(tracks_frame, wig_file_out[[2]]) 
    }    
  }
  dp_points <- gff_regions[c(1,3,4,5,10)]
  ff1 <- group_by (final_frame, name, wig_file)
  # Mean is col mean of the group as defined by groupings. 
  ff2 <- mutate(ff1,Freq = Freq/max(Freq))  
  ff3 <- ungroup(ff2)
  multi_df <- split(data.frame(ff3),f = ff3$name)
  for (i in 1:length(multi_df)){
    df <- multi_df[[i]]
    df <- df[df$Freq >0,]
    name <- multi_df[[i]][1,3]
    tf <- tracks_frame[tracks_frame$name == name,]
    dpp <- dp_points[dp_points$name == name,]
    dpp$mean <- mean(rowMeans(dpp[3:4]))
    breaks <- seq(0, 1, by = 0.1)
    chr_anno <- dpp[1,]
    dpp$V3 <- gsub(x = dpp$V3, pattern = "_", replacement = " ")
    plot <- ggplot(data = multi_df[[i]], aes(x = number_list , y =Freq, colour = wig_file))+
      geom_line(position = position_jitter(w =0.02, h=0))+
      geom_segment(data = tf, aes(x = starts, xend = ends, y = -0.2, yend= -0.2), colour = "blue")+
      geom_segment(data = dpp, aes(x = V4, xend = V5+20, y = -0.4, yend= -0.4), size = 5, colour = "red")+
      geom_text(data = dpp, aes(x = V4,  y = -0.6, label = V3),  colour = c("red"))+
      #facet_wrap(~name, scales = "free")+
      xlab(paste(chr_anno[1,1]))+
      ylab("Coverage (normalised to max in this window)")+
      ggtitle(name)+
      theme_classic()+
      geom_hline(yintercept = 0, size = 0.5, alpha = 0.5)+
      scale_colour_manual(values = c("red", "orange", "purple", "blue"), name = "Time point")
      if (tf[1,"strands"] == "-"){
        scale_x_reverse()        
      }
    ggsave(plot = plot, paste0("plots/",name,".png"), width = 15)
    print(plot)    
  }
  
  
}



d0_1_for <- "/data/home/apattison/Bioinformatics/work_outsde_my_project/dp_ipsc/bigwigs/test/D170e2P2-d0_ATCACG_sorted_for.bw"
d0_2_for <- "/data/home/apattison/Bioinformatics/work_outsde_my_project/dp_ipsc/bigwigs/test/D170e3P2-d0_CTTGTA_sorted_for.bw"
d0_1_rev <- "/data/home/apattison/Bioinformatics/work_outsde_my_project/dp_ipsc/bigwigs/test/D170e2P2-d0_ATCACG_sorted_rev.bw"
d0_2_rev <- "/data/home/apattison/Bioinformatics/work_outsde_my_project/dp_ipsc/bigwigs/test/D170e3P2-d0_CTTGTA_sorted_rev.bw"
ips8_1_for <- "/data/home/apattison/Bioinformatics/work_outsde_my_project/dp_ipsc/bigwigs/test/D170e2P2-iPSp8_CAGATC_sorted_for.bw"
ips8_2_for <- "/data/home/apattison/Bioinformatics/work_outsde_my_project/dp_ipsc/bigwigs/test/D170e3P2-iPSp8_GTGAAA_sorted_for.bw"
ips8_1_rev <- "/data/home/apattison/Bioinformatics/work_outsde_my_project/dp_ipsc/bigwigs/test/D170e2P2-iPSp8_CAGATC_sorted_rev.bw"
ips8_2_rev <- "/data/home/apattison/Bioinformatics/work_outsde_my_project/dp_ipsc/bigwigs/test/D170e3P2-iPSp8_GTGAAA_sorted_rev.bw"

bam_list <- list(d0_1_for, d0_2_for, d0_1_rev,d0_2_rev ,ips8_1_for,ips8_2_for, ips8_1_rev, ips8_2_rev)
name_list <- list ("Day 0 forward strand 1", "Day 0 forward strand 2", "Day 0 reverse strand 1", "Day 0 reverse strand 2",
              "Day 8 forward strand 1", "Day 8 forward strand 2", "Day 8 reverse strand 1", "Day 8 reverse strand 2")

test_gff <- "ipsp8_pass_filter.gff"
bam_dir <- "/data/home/apattison/Bioinformatics/work_outsde_my_project/dp_ipsc/bigwigs/test"
test_out <- grab_regions(test_gff)
run_program(test_out, bam_list, name_list)


# Show chr fix legend change colours by replicate
# Pile together a bunch of ggplots in R
# grid extra, grid.arrange will multi plot