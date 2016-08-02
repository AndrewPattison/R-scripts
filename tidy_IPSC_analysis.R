library ("dplyr")
library("tidyr")
library("readr")
library("ggplot2")

# The goal of this "tidy data" analysis was to first tidy, and the compare and plot 
# the output of Minni Anko's IPSC dataset from the DaPars APA calling tool.  

# Make sure I am in the correct working directory for this analysis. 
setwd("/home/bigpatto2/dp_ipsc/tidy_data_analysis")

# Set a consistent theme for plotting
theme <- theme(panel.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
         text = element_text(size=10), axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"),
         strip.background = element_blank(), panel.border = element_rect(colour = "black",fill = NA))

# Function to combine all the csvs of the same type into a single data frame
combine_data_frames <-function(file_path_list,read_file_list){
  big_frame <- data.frame()
  for (i in 1:length(read_file_list)){
    name <- basename(file_path_list[[i]])
    named_frame <- mutate(read_file_list [[i]], fname = name)
    if (i == 1){
      big_frame <- named_frame
    }
    else{
      big_frame <- rbind(big_frame, named_frame)
    }
  }
  # Split up the names column based on the dot "." to remove file extensions
  big_frame <- separate(big_frame, fname, into =  c("name", "extension"), 
                        sep = "\\.",
                        remove =T) %>%
    select(-extension)
  return(big_frame)
}

# Read in all the file groups into lists.
end_shift_csvs <- list.files(path = "End_shift_csvs/", 
                             pattern = "*.csv", 
                             full.names = T)

end_shift_read <- lapply(end_shift_csvs, read_csv)

raw_dapars_outs <- list.files(path = "raw_dapars_outputs/", 
                             pattern = "*.txt", 
                             full.names = T)

raw_dapars_read <- lapply(raw_dapars_outs, read_delim, "\t")


validated_dapars_outs <- list.files(path = "manually_validated_calls/", 
                             pattern = "*.csv", 
                             full.names = T)

validated_dapars_read <- lapply(validated_dapars_outs, read_csv)

# Combine and clean up each data set
end_shift_messy <- combine_data_frames (end_shift_csvs, end_shift_read)

raw_dapars_messy <- combine_data_frames (raw_dapars_outs, raw_dapars_read) %>%
  separate(Gene, into =  c("transcript_id", "gene", "chromosome", "strand"), 
           sep = "\\|", remove =T) %>%
  select(-gene, -chromosome, -strand) %>%
  separate(Loci, into =  c("chromosome", "start",  "end"), 
           sep = ":|-", remove =T)%>%
  select(-chromosome) %>%
  separate(transcript_id, into =  c("transcript_id", "extension"), 
           sep = "\\.", remove =T) %>%
  select(-extension) %>%
  mutate(Pass_Filter = factor(Pass_Filter, levels = c("Y", "N"))) %>%
  mutate(colours_for_plot = as.character(PDUI_Group_diff >0)) %>%
  unite(colours_for_plot, colours_for_plot , Pass_Filter, sep = "_") %>%
  mutate(name = factor(name, levels = 
                         c("Day_3", "Day_6", "Day_9", "Day_12", "IPSC_1","IPSC_8"))) 

validated_dapars_messy <- combine_data_frames (validated_dapars_outs, validated_dapars_read) %>%
  rename(transcript_id = Transcript_id) %>%
  separate(transcript_id, into =  c("transcript_id", "extension"), 
           sep = "\\.", remove =T) %>%
  select(-extension, -Gene_name) %>%
  mutate(Real_no_d = ifelse (Real == "D", "Y",Real), direction = Change_in_distal_usage > 0) %>%
  mutate(name = factor(name, levels = 
                         c("Day_3", "Day_6", "Day_9", "Day_12", "IPSC_1","IPSC_8"))) 

# Join the DaPars and end shift datasets 
dp_v_end_shift <- left_join(end_shift_messy, raw_dapars_messy, by = c("name", "transcript_id")) %>%
  left_join(validated_dapars_messy, by = c("name", "transcript_id"))

# Join the raw and validated DaPars datasets
dp_v_validated <- left_join(raw_dapars_messy, validated_dapars_messy, by = c("name", "transcript_id")) %>%
  group_by(grp = Real_no_d== "Y") %>%
  unite(col = united, grp,direction,sep="_") %>%
 # filter(adjusted.P_val < 0.05) %>%
  arrange( united)

# Get counts of true and false variables. 
count(x = dp_v_validated,name, PDUI_Group_diff > 0)

# Plot total and validated DaPars results.
ggplot(data = dp_v_validated, aes(x = Group_A_Mean_PDUI *100, y = Group_B_Mean_PDUI *100, colour = united))+
  facet_wrap(~name,nrow = 2)+
  geom_point(size = 0.7)+
  geom_abline()+
  theme+
  scale_colour_manual(values = c("lightgrey", "lightgrey","lightgrey", "green", "blue"))+
  #scale_alpha_manual(values = c(0.1, 0.1, 0.1, 1, 1))+
  guides(colour = F, alpha = F)+
  coord_fixed()+  
  ggsave("new_points_plot_eps.eps")

# Make a stacked bar plot based on my manual validations of the DaPars data
validated_dapars_for_ggplot <- filter(validated_dapars_messy, Real != "D") %>%
  mutate(Real = factor(Real, levels = c("Y", "M", "N"))) %>%
  group_by(name, Real, direction) %>%
  # n() means sumber of rows in group
  summarise(count = n()) %>%
  mutate(direction = ifelse (direction, yes = "down", no = "up")) %>%
  spread(key = direction, value = count) %>%
  mutate(down = down *-1)

# Make a stacked barplot of validated DaPars results. 
ggplot(validated_dapars_for_ggplot, aes(x= name, fill = Real, y = up))+ geom_bar(stat= "identity")+
  geom_bar(aes(x= name, fill = Real, y = down), stat= "identity")+
  geom_abline(slope = 0)+
  ylab("Number of significnat APA events")+  
  xlab("Day (vs MEF)")+
  ggtitle("APA directional trend")+
  #guides(fill = F)+
  theme_bw()+
  theme

# Plot the common direction of APA change. 
ggplot(validated_dapars_for_ggplot, aes(x= name, fill = Real == "Y", y = up))+ geom_bar(stat= "identity")+
  geom_bar(aes(x= name,  fill = Real == "Y", y = down), stat= "identity")+
  geom_abline(slope = 0)+
  ylab("Number of significnat APA events")+  
  labs("Day (vs MEF)")+
  scale_fill_manual (values = c("green", "blue"))+
  ggtitle("APA directional trend")+
  guides(fill = F)+
  theme_bw()+
  theme +
  ggsave("apa_direction_no_leg.pdf")
  

# Grab a random selection of transcript ids 
random_ids <- sample(unique(validated_dapars_messy$transcript_id),size = 40)

#  
validated_dapars_for_heatmappy_typey_thingy <- filter(validated_dapars_messy, Real != "D") %>%
  filter(transcript_id %in% random_ids) %>%
  mutate(direction = ifelse (direction, yes = "down", no = "up")) 

ggplot(data = validated_dapars_for_heatmappy_typey_thingy, aes(x= name , y= transcript_id, colour = direction))+
  geom_point()

# Make a tidy data frame of DaPars calls suported by PAT-Seq
validated_dapars_ps_sup <- filter(validated_dapars_messy, Real != "D") %>%
  mutate(`PAT-Seq supported` = factor(`PAT-Seq supported`, levels = c("Y", "FR", "prox", "dist", "N"))) %>%
  group_by(name, `PAT-Seq supported`, direction) %>%
  # n() means sumber of rows in group
  summarise(count = n()) %>%
  mutate(direction = ifelse (direction, yes = "down", no = "up")) %>%
  spread(key = direction, value = count)

# Make a stacked barplot of DaPars reads supported by PAT-Seq
ggplot(validated_dapars_ps_sup, aes(x= name, fill = `PAT-Seq supported`, y = up))+ 
  geom_bar(stat= "identity")+
  geom_bar(aes(x= name, fill = `PAT-Seq supported`, y = down * -1), stat= "identity")+
  geom_abline(slope = 0)+
  ylab("Number of PAT-Seq supported APA events")+
  xlab("Day (vs MEF)")+
  #scale_fill_manual (values = c("blue", "lightblue", "grey", "lightgreen", "green"))+
  ggtitle("PAT-Seq Support")+
  guides(fill = F)+
  theme_bw()+
  theme+
  ggsave("PAT-Seq_support_no_leg.pdf")

# Exploratory correlation of DaPars calls and Paul Harrison's end-shift method.
ggplot(data= dp_v_end_shift, aes(x= PDUI_Group_diff *-1,  y= r)) + geom_point()

cor.test(dp_v_end_shift$r,  dp_v_end_shift$PDUI_Group_diff*-1)
