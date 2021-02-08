# Authors: Nayela Zeba, Timothy D. Berry, Thea L. Whitman
# Modified from: An open-source, automated, gas sampling peripheral 
# for laboratory incubation experiments; Authors: Timothy D. Berry, Thea L. Whitman

#Load necessary packages
library (tidyverse)
library (zoo)
library (stats)
library (RColorBrewer)
library (broom)

#The working directory will need to be set by the user to point to where their
#data is stored

#setwd("C:/Users/nayel/Box Sync/WhitmanLab/Projects/Biochar_ageing_NZ/incubations_feb2020")

#Changes overflow options to suppress automatic conversion to scientific notation
options(scipen=99)

################################################################################
#Parameter setup
################################################################################

#When smoothing data, how many points wide should the sliding window be? Larger values = smoother
smooth_window = 50

#When checking for outliers, how many points should be analyzed (nslice) and what
#threshold constitutes an outlier (outliercutoff) smaller values will find more 
#outliers, but increase chances of false positives. It is recommended these values
#not be changed unless you have very choppy data even after smoothing
nslice = 10
outliercutoff = 8

#Designates the valves used for flushing manifolds (flush_positions)
#and for moving gas through system without sampling jars (idle_positions).
#Anything not listed here will be treated like a sample!
flush_positions <-  c(1,17,33,49)
idle_positions <-  c(16,32,48,63)

#Fraction of system volume that comes from the sample jar. In a system with no-deadspace 
#this would be 1. In reality, this depends on jar volume and
#must be determined for each specific sample container that is used
fj = 0.85

#Parameters specific to incubation experimental design:
#Jar volume (L)
jarvol_L = 0.12

################################################################################
#Data preprocessing
################################################################################

#Assembles a list of log files and Picarro data files from the chosen directory and combines them
logfile_list <- list.files(path="./raw_data", pattern = 'log_*', recursive = T, full.names = TRUE)
logfile_list_read = lapply(logfile_list, read.csv, header = TRUE)
relay_log = do.call(rbind, logfile_list_read)

picarro_files <- list.files(path = "./raw_data", pattern='.*[.]dat$', recursive = T, full.names = TRUE)
picarro_files

all.pic.data <- as.data.frame(c())
total.pic.files <- length(picarro_files)

for(i in 1:total.pic.files){
  temp <- read_table(picarro_files[i])
  all.pic.data <- rbind.data.frame(temp, all.pic.data)
}

rm(i, picarro_files)

#Keep only data columns that are needed for this analysis
columns_to_keep <- c("TIME", "EPOCH_TIME", "12CO2_dry")

cleaned.pic.data <- all.pic.data[,columns_to_keep]
rm(all.pic.data, temp)

#Removes unneeded columns and times from data that are from 
#before or after the relay was active (since it's not experiment data)
#Also renames the CO2 concentration column to not start with
#a number so that GG plot doesn't have a fit later
relay_log= relay_log[order(relay_log$Epoch_time,decreasing = FALSE),]
pic.trimmed <-  cleaned.pic.data[cleaned.pic.data$EPOCH_TIME > relay_log$Epoch_time[1],]
pic.trimmed= pic.trimmed[order(pic.trimmed$EPOCH_TIME),]
pic.trimmed <- pic.trimmed[pic.trimmed$EPOCH_TIME < relay_log$Epoch_time[length(relay_log$Epoch_time)],]
names(pic.trimmed)[3] <- "dry_12CO2"


#Coverts epoch times into integers and removes duplicate timestamps caused by
#the Picarro and relay board not measuring at 1 second intervals. 
#This results in a loss of ~30% of the data points, but there should still be plenty 
#left to redraw the curves and do our calculations with
pic.trimmed <- dplyr::mutate(pic.trimmed, EPOCH_TIME = as.integer(EPOCH_TIME))
pic.trimmed.deduped <- dplyr::distinct(pic.trimmed, EPOCH_TIME, .keep_all = TRUE)

#Renames the Epoch_time in the relay file to EPOCH_TIME,
#removes duplicate time stamps, and creates a data frame for the
#relay logs wherein the epoch time matches with a picarro time
names(relay_log)[1] <- "EPOCH_TIME"
relay.deduped <- distinct(relay_log, EPOCH_TIME, .keep_all = TRUE)
relay.matched <- dplyr::filter(relay.deduped, EPOCH_TIME %in% pic.trimmed.deduped$EPOCH_TIME)

merged <-  merge.data.frame(relay.matched, pic.trimmed.deduped, by.x = "EPOCH_TIME")

#Cleanup to salvage some memory
rm(cleaned.pic.data, pic.trimmed, relay.deduped, relay_log, logfile_list_read, logfile_list,
   pic.trimmed.deduped, relay.matched)

#saveRDS(merged, "r_objects/merged")

#########################################
### Start here if using previously merged file
#########################################
# merged=readRDS("r_objects/merged")

#Smoothing using the ZOO package. Takes a moving average (determined by the smooth_window parameter, set above). 
# For this dataset, 300 points roughly corresponds to a 5 minute moving average;
# 50 points corresponds to a 50 second moving average.
CO2_smooth <-rollapply(data = merged$dry_12CO2, FUN = "mean", width = smooth_window, partial = TRUE, align = "right")
merged <- cbind(merged,CO2_smooth)

rm(CO2_smooth)

################################################################################
#Data Annotation (adds status and cycle info)
################################################################################

#Remove the redundant Time column in the merged dataframe and creates status, 
#sample number,and cycle number columns which are then assigned dummy values
merged$Time <- NULL
merged$status <- 0
merged$sample_num <- 0
merged$cycle <- NA
merged$cycle_bound <- FALSE


#Determines the status of the analysis (stored in the merged dataframe) at each
#time point and assigns sample numbers to jar-specific steps. The status column is 
#mutated to replace the dummy value with a status value assigned based on
#the combination of valves open. For the sake of this process, it is assumed that the default 
#manifold configuration is used and the first position on each manifold is used to 
#flush while the last position on each is an "idle" position in which no jar is sampled.

merged <- merged %>% 
  dplyr::mutate(status = replace(status, is.element(Active_relay1, flush_positions) & is.element(Active_relay2, idle_positions),"flushing lines")) %>%
  dplyr::mutate(status = replace(status, is.element(Active_relay1, flush_positions) & is.na(Active_relay2),"flushing lines")) %>%
  dplyr::mutate(status = replace(status, is.element(Active_relay1,flush_positions) & !is.element(Active_relay2,idle_positions) & !is.na(Active_relay2),"flushing jar")) %>% 
  dplyr::mutate(status = replace(status, !is.element(Active_relay1,flush_positions) & !is.element(Active_relay1,idle_positions),"measuring jar")) %>%
  dplyr::mutate(status = replace(status, is.element(Active_relay1,idle_positions),"idle")) %>%
  dplyr::mutate(status = replace(status, status != lag(status,1),"boundary")) %>%
  dplyr::mutate(status = replace(status, is.na(lag(Date,1)),"boundary")) %>%
  dplyr::mutate(sample_num = case_when(.$status == "measuring jar" ~ .$Active_relay1,
                                .$status == "flushing jar" ~ .$Active_relay2,
                                .$status == "boundary" ~ as.integer(0),
                                .$status == "flushing lines" ~ as.integer(0),
                                .$status == "idle" ~ as.integer(0))) %>%
  dplyr::mutate(cycle_bound = 
           replace(cycle_bound, status == "boundary" &
                     lag(status) == "flushing jar" &
                     lead(status) == "idle", TRUE))


#!!!NOTE: This function assigns cycle numbers based on boundaries between each 
#measurement cycle. When using this function it is assumed that cycles were run 
#to completetion and following the suggested analysis timing template. If this
#is not the case, this function must be modified to accurately label cycle numbers

cycle_counter <- function(x){
                  boundary_vec <- which(x$cycle_bound == TRUE)
                  boundary_vec <- c(1, boundary_vec)
                  #print(boundary_vec)
                  for(i in 2:length(boundary_vec)){
                   # print(i)
                   x[c(boundary_vec[i-1]:boundary_vec[i]),]$cycle <- i -1 
                  }
                  return(x)
}

merged_complete<-cycle_counter(merged)

rm(merged)


################################################################################
#Summarize Data, correct for dilution, and adjust for standards  
################################################################################

#Remove boundary and idle statuses
jars <- merged_complete %>% group_by(sample_num, cycle) %>% 
  filter(cycle != 0) %>% filter(status != "boundary") %>% filter(status != "idle")

rm(merged_complete)

########

#In this version of the dilute selection, the minimum point in 
#the CO2 measurement (the most dilute) has been taken and the points within nslice/2 of it have been selected

most_dilute = jars %>% filter(status == "measuring jar") %>% 
  group_by(sample_num, cycle) %>% summarize(cycle_minimum = min(CO2_smooth))

min_time <- jars %>% filter(status == "measuring jar") %>%  
  group_by(sample_num,cycle) %>%  merge(.,most_dilute) %>%
  filter(CO2_smooth == cycle_minimum) %>% 
  select(sample_num, cycle, minimum_time = EPOCH_TIME)


CO2_dilute = jars%>% 
  # Grab the datapoints right when the jar is being measured
  filter(status == "measuring jar") %>%
  # Create a series of columns identifying the difference between neighboring values,
  # This is not really necessary for the smoothed data but is not computationally
  #intensive so there is little reason to remove it yet.
  select(sample_num,cycle,CO2_smooth, EPOCH_TIME,status) %>%
  mutate(NeighborDist1 = c(NA,diff(CO2_smooth,lag=1)),
         NeighborDist2 = c(NA,NA,diff(CO2_smooth,lag=2)),
         NeighborDist3 = c(NA,NA,NA,diff(CO2_smooth,lag=3)),
         NeighborDist4 = c(diff(CO2_smooth,lead=1),NA),
         NeighborDist5 = c(diff(CO2_smooth,lead=2),NA),
         NeighborDist6 = c(diff(CO2_smooth,lead=3),NA)) %>%
  # Ungroup everything
  group_by(EPOCH_TIME)%>%
  # Create an outlier flag column, that pings when the distance between nearby points is above the outlier cutoff
  mutate(Flag = ifelse(max(abs(NeighborDist1),abs(NeighborDist2),abs(NeighborDist3),abs(NeighborDist4),abs(NeighborDist5),abs(NeighborDist6),na.rm=TRUE)>outliercutoff,1,0))%>%
  # Group them up again
  group_by(sample_num,cycle) %>%  merge(.,min_time) %>%  
  mutate(diff = minimum_time - EPOCH_TIME) %>% filter(abs(diff) <= nslice/2) %>%
  filter(Flag != 1) %>% 
  group_by(sample_num,cycle) %>%
  # Take the mean value of that raw CO2 input
  # End result is table with sample number, cycle number, and mean CO2
  summarize(dilute_CO2=mean(CO2_smooth))


#Do the same thing at the end of the flush step to get the purge values (concentration)
#of CO2 in system after flushing. This is also where the summarized sample's time 
#entry comes from - the last point in the purge step is when accumulation begins.
CO2_purged = jars%>%
  filter(status == "flushing jar") %>%
  select(sample_num,cycle,CO2_smooth, EPOCH_TIME,status) %>%
  mutate(NeighborDist1 = c(NA,diff(CO2_smooth,lag=1)),
         NeighborDist2 = c(NA,NA,diff(CO2_smooth,lag=2)),
         NeighborDist3 = c(NA,NA,NA,diff(CO2_smooth,lag=3)),
         NeighborDist4 = c(diff(CO2_smooth,lead=1),NA),
         NeighborDist5 = c(diff(CO2_smooth,lead=2),NA),
         NeighborDist6 = c(diff(CO2_smooth,lead=3),NA)) %>%
  group_by(EPOCH_TIME)%>%
  mutate(Flag = ifelse(max(abs(NeighborDist1),abs(NeighborDist2),abs(NeighborDist3),
                           abs(NeighborDist4),abs(NeighborDist5),abs(NeighborDist6)
                           ,na.rm=TRUE)>outliercutoff,1,0)) %>%
  group_by(sample_num,cycle)%>%
  slice((n()-nslice+1):n())   %>%
  filter(Flag != 1) %>%
  summarize(purged_CO2=mean(CO2_smooth), time = last(EPOCH_TIME))

#Peak selection takes the last point of the "measuring jar" status and points preceding
#this point by 10 seconds to be representative of the peak value: this tends 
#to represent the most stable part of the "peak" and makes the measurement
#much less likely to be erroneous when the previous sample was very concentrated

meas_end <- jars %>% group_by(sample_num, cycle) %>% 
  filter(status == "measuring jar") %>% summarize(last_point = last(EPOCH_TIME))

CO2_peak = jars%>%
  filter(status == "measuring jar") %>%
  select(sample_num,cycle,CO2_smooth, EPOCH_TIME,status) %>%
  mutate(NeighborDist1 = c(NA,diff(CO2_smooth,lag=1)),
         NeighborDist2 = c(NA,NA,diff(CO2_smooth,lag=2)),
         NeighborDist3 = c(NA,NA,NA,diff(CO2_smooth,lag=3)),
         NeighborDist4 = c(diff(CO2_smooth,lead=1),NA),
         NeighborDist5 = c(diff(CO2_smooth,lead=2),NA),
         NeighborDist6 = c(diff(CO2_smooth,lead=3),NA)) %>%
  group_by(EPOCH_TIME) %>%
  mutate(Flag = ifelse(max(abs(NeighborDist1),abs(NeighborDist2),abs(NeighborDist3),
                           abs(NeighborDist4),abs(NeighborDist5),abs(NeighborDist6),
                           na.rm=TRUE)>outliercutoff,1,0)) %>% 
  group_by(sample_num, cycle) %>% merge(.,meas_end) %>%
  mutate(diff = last_point - EPOCH_TIME) %>% filter(diff <= nslice) %>%
  filter(Flag != 1) %>% 
  group_by(sample_num,cycle) %>%
  summarize(peak_CO2 = mean(CO2_smooth))

####Check your cycle assignments####
test <-jars %>% filter(cycle == 8, sample_num < 25)

ggplot(data = test, aes(x = EPOCH_TIME, y = CO2_smooth, color = status)) +
  geom_point() + facet_wrap(~cycle, scale = "free")

#Merge these 3 sets of values together and remove unnecessary intermediate files  
summary <- merge(CO2_dilute, CO2_peak) %>%
  merge(.,CO2_purged) %>%
  arrange(sample_num,cycle)

summary_long <- summary %>% rename("dilute" = dilute_CO2, "peak" = peak_CO2, "purged" = purged_CO2) %>%
  gather(key = "landmark", value = "CO2", dilute, peak, purged)

#Plot figures to see the CO2 concentration across cycles
#This is an easy way to visualize trends in respiration of samples and to assess
#whether your relay command file is adequate. Ideally the dilute and purge values
#should remain consistent throughout the course of the experiment. Changes in
#these values over the course of the experiment likely mean insufficient time is
#being given for the "flushing jars" step. The script can account for this deviations
#but to maximize data quality the flush time should be increased.

respiration_profiles <- ggplot(data = summary_long) + 
  geom_col(aes(x = cycle, y = CO2, fill = landmark), position = "dodge") + 
  facet_wrap(~sample_num, scale = "free")

respiration_profiles

#Remove unnecessary intermediate files to free up memory

rm(meas_end, min_time, most_dilute)


##########################################
#Calculations
##########################################

#Dilution parameters are determined by injecting known volumnes of pure CO2 into
#purged vessels to determine the fraction of the contribution of the sample jar's gas (fj) to observed gas

#The total CO2 measured at the peak represents a combination of what's in the system
#before you start measuring and what's in the jar. The total volume of the system is jar
#+ gas lines, instrument, manifolds, etc. Thus, the CO2 in an equilibrated system will 
#be whatever  CO2 came from the system (dilute_CO2) multiplied by its fraction (1-fj),
#plus the CO2 that was in the jar multiplied by its fraction (fj).
#That results in the following equation:

corrected <- summary %>% mutate(peak_CO2_cor = (peak_CO2 - dilute_CO2 + (dilute_CO2 * fj))/fj)

#Subtract the background purged CO2 (previous cycle) from the corrected peak CO2 
#to get the CO2 accumulation from PyOM mineralization 
corrected <- corrected %>% 
  group_by(sample_num) %>% 
  arrange(cycle) %>%
  mutate(final_CO2 = peak_CO2_cor - lag(purged_CO2))
  
# Covert CO2 values from ppm to a mass basis for control CO2 data.
# We make some simple assumptions that the 
# sample is at standard temperature and pressure
# We also adjust time into time since incubations started and convenient units
df = corrected %>%
  dplyr::select(sample_num,cycle,time,final_CO2)%>%
  dplyr::group_by(sample_num) %>%
  dplyr::mutate(time_s = time-min(time)) %>%
  dplyr::mutate(time_min = time_s/60)%>%
  dplyr::mutate(time_hr = time_min/60)%>%
  dplyr::arrange(time_hr)%>%
  dplyr::mutate(CO2C_mg = final_CO2/1000000/22.4*jarvol_L*12.01*1000)

#Reading in the jar position information. Which valves are connected to samples that are inoculated or controls or blanks?
treatments=data.frame(read.csv("treatments_Feb.csv"))

#Merging with the CO2 dataframe
df = merge(df,treatments,by="sample_num")

#Cleaning up the dataframe, removing unwanted controls and setting the CO2 mg mineralized in the 1st cycle to 0
df$sample_num = as.factor(df$sample_num)
df$sample_num = as.factor(df$sample_num)
unwanted_levels = c("350_BIO_CON_C_1", "350_BIO_CON_C_2", "550_BIO_CON_C_1", "550_BIO_CON_C_2")
df[df$cycle==1,]$CO2C_mg = 0

#Calculating cumulative CO2-C mineralized over each measurement time CO2 data
df.c = df%>%
  dplyr::filter(!sample_name %in% unwanted_levels) %>%
  dplyr::group_by(sample_num)%>%
  dplyr::arrange(cycle)%>%
  dplyr::mutate(Total_CO2C_mg.cumulative = cumsum(CO2C_mg))

#saveRDS(df.c, "../global_r_objects/cumulative_co2cmg")

#################################################################################
#Plotting
#################################################################################

######## Figure 2 #############################

library(ggplot2)
library(wesanderson)


#Creating confidence interval functions 

lower_ci <- function(mean, se, zvalue = 1.960){
  res <- mean - (zvalue * se)
  return(res)
}
upper_ci <- function(mean, se, zvalue = 1.960){
  res <- mean + (zvalue * se)
  return(res)
}

# reading CO2 data 
df.co2 = readRDS("../global_r_objects/cumulative_co2cmg")
# replicate is a reserved word, so rename the column 
df.co2$repl = df.co2$replicate
df.co2$replicate =NULL
# remove white space
df.co2$treatment = trimws(df.co2$treatment)


# create a data frame with means of blank replicates within each treatment, temperature and cycle 
df.blank.means = df.co2 %>%
  dplyr::filter(blank == TRUE) %>%
  dplyr::filter(treatment != "Blank") %>%
  dplyr::group_by(cycle, treatment, temperature) %>%
  dplyr::summarise(mean_blank = mean(Total_CO2C_mg.cumulative))

# subtract the blank mean from the sample mean within each treatment, temperature and cycle
df.co2.corr = df.co2 %>%
  dplyr::filter(blank == FALSE) %>%
  dplyr::filter(treatment != "Blank") %>% 
  dplyr::left_join(df.blank.means) %>%
  dplyr::mutate(corr_co2c_mg.cumulative = Total_CO2C_mg.cumulative - mean_blank)

df.co2.corr$treatment= as.factor(df.co2.corr$treatment)
df.co2.corr$temperature= as.factor(df.co2.corr$temperature)
#saveRDS(df.co2.corr, "../global_r_objects/corr_cumulative_co2cmg")

#normalize with biochar C

# read elemental data
elemental.data = read.csv("raw_data/elemental.csv")


# calculating means and S.Ds
c.means = elemental.data %>%
  dplyr::select(- c(tot_n, tot_h,  ash)) %>%
  dplyr::group_by(treatment,temperature) %>%
  dplyr::summarize(mean_c = mean(tot_c),
                   sd_c = sd(tot_c))

# read biochar mass and media info
b.mass = read.csv("raw_data/biochar_mass.csv")

c.means = merge(c.means, b.mass, by= c("treatment", "temperature"))
c.means$treatment = trimws(c.means$treatment)
c.means$treatment= as.factor(c.means$treatment)
c.means$temperature= as.factor(c.means$temperature)

# calculate C content in biochar mass used
c.means = c.means %>%
  dplyr::mutate (c_mass_g = (biochar_mass_g*mean_c)/100) %>%
  dplyr::mutate (c_conc = c_mass_g/media_ml) %>%
  dplyr::mutate (c_g_per_jar = c_conc*40) %>%
  dplyr::select(-media_ml, -biochar_mass_g, -c_mass_g, -c_conc)

# merge c.means with df.co2.corr and calculate normalized values

df.co2.norm =  df.co2.corr %>%
  dplyr::left_join(c.means) %>%
  dplyr::mutate(norm_co2c_mg.cumulative = corr_co2c_mg.cumulative/c_g_per_jar)
#saveRDS(df.co2.norm, "../global_r_objects/norm_co2cmg")

#plotting the CO2 respired time-series
# calculate the means and C.Is at each time point
df.co2.norm.mean = df.co2.norm %>%
  dplyr::group_by(cycle,treatment,temperature) %>%
  dplyr::mutate(
    mean_co2c_mg_cumulative = mean(norm_co2c_mg.cumulative),
    sd_co2c_mg.cumulative = sd(norm_co2c_mg.cumulative),
    mean_time_hr = mean(time_hr),
    n_co2c_mg.cumulative = n()) %>%
  dplyr::mutate(se_co2c_mg.cumulative = sd_co2c_mg.cumulative/sqrt(n_co2c_mg.cumulative),
                lower_ci_co2c_mg.cumulative = lower_ci(mean_co2c_mg_cumulative, se_co2c_mg.cumulative),
                upper_ci_co2c_mg.cumulative = upper_ci(mean_co2c_mg_cumulative, se_co2c_mg.cumulative))


p = ggplot(df.co2.norm.mean)
p = p + geom_errorbar(aes(x=mean_time_hr,
                          ymin=lower_ci_co2c_mg.cumulative,
                          ymax=upper_ci_co2c_mg.cumulative,
                          colour = treatment),
                      width=5)
p = p + geom_point(aes(x = mean_time_hr, y = mean_co2c_mg_cumulative, shape = treatment, colour = treatment), size = 3)
p = p + geom_line(aes(x = mean_time_hr, y = mean_co2c_mg_cumulative, linetype = treatment, colour = treatment), size = 1)
p = p + xlab("Time since inoculation (hours)") 
p = p + ylab(expression(atop(Cumulative~C~mineralized, 
                             ~(mg~CO[2]-C/g~initial~biochar~C)))) 
p = p + theme_bw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p = p + scale_shape_manual(values=c(15,19,17,18))
p = p + scale_colour_manual(values = wes_palette("Darjeeling1", n = 4))
p = p + scale_linetype_manual(values = c("dashed", "dotted", "dotdash", "solid"))
p = p + theme(legend.title = element_blank(),legend.text=element_text(size=14), legend.position = "top")
p = p + theme(axis.text=element_text(size=14),axis.title = element_text(size=14),strip.text=element_text(size=16))
p = p + facet_wrap(~temperature)
p

#ggsave("../plots/Rplot_incubation_norm.png", dpi=300, device="png", width=8, height=7)

# total co2 respired over the incubation period
total.co2 = df.co2.norm %>% #change the data frame to df.co2.corr if needed                      
  dplyr::filter(cycle==8) %>%
  dplyr::select(sample_name, corr_co2c_mg.cumulative, norm_co2c_mg.cumulative, treatment, temperature, repl, c_g_per_jar) %>%
  dplyr::select(-sample_num)

#saveRDS(total.co2, "../global_r_objects/total_norm_co2cmg")