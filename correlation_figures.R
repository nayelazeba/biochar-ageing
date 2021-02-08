#Authors: Nayela Zeba, Thea L. Whitman

#Load necessary packages
library(ggplot2)
library(dplyr)
require(tidyr)
require(plyr)
library(wesanderson)


# Reading in the cumulative CO2-c mg data
total.co2 = readRDS("../global_r_objects/total_norm_co2cmg")

################# Tukey tests for CO2-C mg data #######################
total.co2.550= total.co2 %>%
  filter(temperature==550)
total.co2.350= total.co2 %>%
  filter(temperature==350)

aov.out = aov(norm_co2c_mg.cumulative ~ treatment, data=total.co2.550)
options(show.signif.stars=F)
summary(aov.out)
#capture.output(summary(aov.out), file= "../plots/anova_results_norm.txt")
TukeyHSD(aov.out)
#capture.output(TukeyHSD(aov.out), file= "../plots/Tukey_results_interaction.txt")


################################################################################################  Correlation between O/C ratio and C mineralized #####################

#Reading in elemental data
elemental.data = read.csv("raw_data/elemental.csv")
elemental.data$treatment = as.factor(elemental.data$treatment)
elemental.data$temperature = as.factor(elemental.data$temperature)

#Calculating means 
elemental.means = elemental.data %>%
  dplyr::group_by(treatment,temperature) %>%
  dplyr::summarize(mean_c = mean(tot_c,), sd_c = sd(tot_c),
                   mean_n = mean(tot_n), sd_n = sd(tot_n),
                   mean_h = mean(tot_h, na.rm = T), sd_h = sd(tot_h, na.rm = T),
                   mean_ash = mean(ash, na.rm = T), sd_ash = sd(ash, na.rm = T))

#Calculating mean O content as per Enders et al., 2012
elemental.means = elemental.means %>%
  dplyr::mutate(mean_o = 100- (mean_c +mean_n + mean_h +mean_ash))

#Calculating the molar O/C and H/C ratio                  
elemental.molar = elemental.means %>%
  dplyr::mutate(moles_mean_c = mean_c/12.011) %>%
  dplyr::mutate(moles_mean_h = mean_h/1.008) %>%
  dplyr::mutate(moles_mean_n = mean_n/14.007) %>%
  dplyr::mutate(moles_mean_o = mean_o/16) %>%
  dplyr::mutate(molar_oc_ratio = moles_mean_o/moles_mean_c) %>%
  dplyr::mutate(oc_ratio = mean_o/mean_c) %>%
  dplyr::mutate(molar_hc_ratio = moles_mean_h/moles_mean_c) %>%
  dplyr::mutate(hc_ratio = mean_h/mean_c)
  
#Merging with the CO2-C mg data frame
elemental.merged = merge(elemental.molar, total.co2, by.x=c("treatment", "temperature"), by.y=c("treatment", "temperature"))

# Separating into two data frames to plot 350 and 550 data sets separately
elemental.merged.350= elemental.merged %>%
  dplyr::filter(temperature==350) %>%
  dplyr::select(-c_g_per_jar, -sample_num)

elemental.merged.550= elemental.merged %>%
  dplyr::filter(temperature==550) %>%
  dplyr::select(-c_g_per_jar, -sample_num)

#Defining function to extract R2 and p- values for the two data sets
mod_eqn = function(df, temperature){
  m1 = lm(norm_co2c_mg.cumulative ~ molar_oc_ratio, df)
  r2 = summary(m1)$r.squared
  m2 = aov(norm_co2c_mg.cumulative ~ molar_oc_ratio, df)
  p_val = summary(m2)[[1]][["Pr(>F)"]][[1]]
  if (p_val < 0.001) {
    eq <- substitute({italic(R)^{2}}[temperature]~"="~r2~","~{italic(p)}~"<"~0.001, 
                     list(r2 = format(r2, digits = 3), temperature = temperature))
  } else {
    eq <- substitute({italic(R)^{2}}[temperature]~"="~r2~","~{italic(p)}~"="~p_val, 
                     list(r2 = format(r2, digits = 3), p_val = format(p_val, digits = 2), temperature = temperature))  
  }
  as.character(as.expression(eq)) 
}

#Extracting R2 and p-vales into a data frame
eq_350 <- ddply(elemental.merged.350, .(temperature), mod_eqn, 350)
eq_550 <- ddply(elemental.merged.550, .(temperature), mod_eqn, 550)

########## Plotting Fig. 3 ####################
el_reg_plot = ggplot(NULL, mapping = aes())

#Adding 350 data
el_reg_plot = el_reg_plot + geom_point(elemental.merged.350, 
                                       mapping = aes(x=molar_oc_ratio,  y=norm_co2c_mg.cumulative, shape=treatment, alpha= temperature), 
                                       colour= "black", size=4)
#Adding 550 data
el_reg_plot = el_reg_plot + geom_point(elemental.merged.550, 
                                       mapping = aes(x=molar_oc_ratio,  y=norm_co2c_mg.cumulative, shape=treatment, alpha= temperature), 
                                       colour = "black", size=4)
el_reg_plot = el_reg_plot + xlab("O/C ratio") 
el_reg_plot = el_reg_plot + ylab(expression(atop(Cumulative~C~mineralized, 
                                                     ~(mg~CO[2]-C/g~initial~biochar~C)))) 
el_reg_plot = el_reg_plot + theme_bw() + theme(panel.grid.major = element_blank(),
                                               panel.grid.minor = element_blank(), 
                                               axis.line = element_line(colour = "black"))
el_reg_plot = el_reg_plot + theme(legend.title = element_blank(),legend.text=element_text(size=14), 
                                  legend.position = "top") 
el_reg_plot = el_reg_plot + theme(axis.text=element_text(size=14),axis.title = element_text(size=14),
                                  strip.text=element_text(size=16))
el_reg_plot = el_reg_plot + scale_shape_manual(values=c(15,19,17,18))
el_reg_plot = el_reg_plot + scale_alpha_manual(values=c(0.4,1))
el_reg_plot = el_reg_plot + scale_linetype_manual(values=c("dashed","solid"))

#Adding regression line for the 350 and 550 data sets
el_reg_plot = el_reg_plot + geom_smooth(method="lm", 
                                        data=elemental.merged.350, level=FALSE, 
                                        mapping=aes(x=molar_oc_ratio, y=norm_co2c_mg.cumulative, linetype = temperature), 
                                        colour = "black", size=1)
el_reg_plot = el_reg_plot + geom_smooth(method="lm", data=elemental.merged.550, level=FALSE, 
                                        mapping=aes(x=molar_oc_ratio, y=norm_co2c_mg.cumulative, linetype = temperature), 
                                        colour = "black", size=1)

#Adding the statistics for the 350 and 550 data sets
el_reg_plot = el_reg_plot + geom_text(data=eq_350, 
                                      aes(x = Inf, y = Inf, label=V1), 
                                      hjust=1.6, vjust=1, parse = TRUE, 
                                      inherit.aes=FALSE, color="black", size = 5)
el_reg_plot = el_reg_plot + geom_text(data=eq_550 , 
                                      aes(x = Inf, y = Inf, label=V1), 
                                      hjust=1.6, vjust=2.2, parse = TRUE, 
                                      inherit.aes=FALSE, colour="black", size = 5)

el_reg_plot

#ggsave("../plots/cmin_norm_vs_oc.png", dpi=300, device="png", width=7.5, height=7)



######################################################################################
#############  Correlation between FTIR peaks and C mineralized #####################

#Reading in the raw FTIR data
ftir.data= read.csv("raw_data/ftir_peaks.csv")

#Removing the unwanted peaks (3350_OH removed to avoid interfrence from water molecules)
#and calculating relative peak heights 
extr = c("3350_OH", "3055_ar_CH")

ftir.data = ftir.data %>%
  dplyr::filter(!peak_name %in% extr) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(sum_height = sum(corr_height)) %>%
  dplyr::mutate(rel_height = corr_height / sum_height) %>%
  dplyr::select(temperature, treatment, peak_name, rel_height)

#Checking the formatting 
ftir.data$temperature = as.factor(ftir.data$temperature)
ftir.data$treatment = as.factor(ftir.data$treatment)
ftir.data$treatment = trimws(ftir.data$treatment)
ftir.data$peak_name = as.factor(ftir.data$peak_name)


#Merging with the CO2-C mg data frame
ftir.merged = merge(ftir.data, total.co2, by.x=c("treatment", "temperature"), by.y=c("treatment", "temperature"))

#Separating into two data frames to plot 350 and 550 data sets separately
ftir.merged.350= ftir.merged %>%
  dplyr::filter(temperature==350) %>%
  dplyr::select(-c_g_per_jar, -sample_num)

ftir.merged.550= ftir.merged %>%
  dplyr::filter(temperature==550) %>%
  dplyr::select(-c_g_per_jar, -sample_num)

#Renaming the facet labels
ftir.merged.350$facets = factor(ftir.merged.350$peak_name, labels = c(
  "paste(CO, '*')~(1200~cm^{-1})",
  "al~C-H~(1413~cm^{-1})",
  "ar~C==C~(1593~cm^{-1})",
  "C==0~(1701~cm^{-1})",
  "al~C-H~(2932~cm^{-1})",
  "ar~C-H~(810~cm^{-1})"
))

ftir.merged.550$facets = factor(ftir.merged.550$peak_name, labels = c(
  "paste(CO, '*')~(1200~cm^{-1})",
  "al~C-H~(1413~cm^{-1})",
  "ar~C==C~(1593~cm^{-1})",
  "C==0~(1701~cm^{-1})",
  "al~C-H~(2932~cm^{-1})",
  "ar~C-H~(810~cm^{-1})"
))

#Defining function to extract R2 and p - values for the two data sets
mod_eqn = function(df, temperature){
  m1 = lm(norm_co2c_mg.cumulative ~ rel_height, df)
  r2 = summary(m1)$r.squared
  m2 = aov(norm_co2c_mg.cumulative ~ rel_height, df)
  p_val = summary(m2)[[1]][["Pr(>F)"]][[1]]
  if (p_val < 0.001) {
    eq <- substitute({italic(R)^{2}}[temperature]~"="~r2~","~{italic(p)}~"<"~0.001, 
                     list(r2 = format(r2, digits = 3), temperature = temperature))
  } else {
    eq <- substitute({italic(R)^{2}}[temperature]~"="~r2~","~{italic(p)}~"="~p_val, 
                     list(r2 = format(r2, digits = 3), p_val = format(p_val, digits = 2), temperature = temperature))  
  }
  as.character(as.expression(eq)) 
}

#Extracting R2 and p-vales into a data frame
eq_350 <- ddply(ftir.merged.350,.(facets), mod_eqn, 350)
eq_550 <- ddply(ftir.merged.550,.(facets), mod_eqn, 550)

########## Plotting Fig. 4 ####################
ftir_reg_plot = ggplot(NULL, mapping = aes())

#Adding 350 data
ftir_reg_plot = ftir_reg_plot + geom_point(ftir.merged.350, 
                                           mapping = aes(x=rel_height,  y=norm_co2c_mg.cumulative,  shape=treatment, alpha = temperature), 
                                           colour = "black", 
                                           size=4)

#Adding 550 data
ftir_reg_plot = ftir_reg_plot + geom_point(ftir.merged.550, 
                                           mapping = aes(x=rel_height,  y=norm_co2c_mg.cumulative,  shape=treatment,  alpha = temperature), 
                                           colour = "black", 
                                           size=4)

ftir_reg_plot = ftir_reg_plot + facet_wrap(~facets, scale="free", labeller = label_parsed)
ftir_reg_plot = ftir_reg_plot + xlab("Relative peak height") 
ftir_reg_plot = ftir_reg_plot + ylab(expression(atop(Cumulative~C~mineralized, 
                                                     ~(mg~CO[2]-C/g~initial~biochar~C)))) 
ftir_reg_plot = ftir_reg_plot + theme_bw() + theme(panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), 
                                                   axis.line = element_line(colour = "black"))
ftir_reg_plot = ftir_reg_plot + scale_shape_manual(values=c(15,19,17,18))
ftir_reg_plot = ftir_reg_plot + scale_alpha_manual(values= c(0.4, 1))
ftir_reg_plot = ftir_reg_plot + scale_linetype_manual(values=c("dashed","solid"))
ftir_reg_plot = ftir_reg_plot + theme(legend.title = element_blank(),
                                      legend.text=element_text(size=16), 
                                      legend.position = "top") 
ftir_reg_plot = ftir_reg_plot + theme(axis.text=element_text(size=16),
                                      axis.title = element_text(size=16),
                                      strip.text=element_text(size=16))

#Adding regression line for the 350 and 550 data sets
ftir_reg_plot = ftir_reg_plot + geom_smooth(method="lm", 
                                            data=ftir.merged.350, 
                                            level=FALSE, 
                                            mapping=aes(x=rel_height, y=norm_co2c_mg.cumulative, linetype = temperature),
                                            colour = "black", 
                                            size=1)
ftir_reg_plot = ftir_reg_plot + geom_smooth(method="lm", 
                                            data=ftir.merged.550, 
                                            level=FALSE, 
                                            mapping=aes(x=rel_height, y=norm_co2c_mg.cumulative, linetype= temperature), 
                                            colour= "black", 
                                            size=1)

#Adding the statistics for the 350 and 550 data sets
ftir_reg_plot = ftir_reg_plot + geom_text(data=eq_350 , 
                                          aes(x = Inf, y = Inf, label=V1), 
                                          hjust=1.5, vjust=1, 
                                          parse = TRUE, 
                                          inherit.aes=FALSE, color="black", size = 5)
ftir_reg_plot = ftir_reg_plot + geom_text(data=eq_550 , 
                                          aes(x = Inf, y = Inf, label=V1), 
                                          hjust=1.5, vjust=2.2, 
                                          parse = TRUE, inherit.aes=FALSE, 
                                          colour="black", size = 5)
ftir_reg_plot

#ggsave("../plots/cmin_norm_ftir_reg_comb.png", dpi=300, device="png", width=14, height=12)


#####################################################################################################
########################### surface colony count #####################################################
#Reading the colony count data
growth.data = read.csv("raw_data/growth_data_full.csv")
growth.data$repl = growth.data$replicate
growth.data$replicate =NULL
growth.data$treatment = trimws(growth.data$treatment)

#plotting supplemetary Fig. S2
growth_plot_all = ggplot(growth.data,aes(x=treatment, y=pct_area))
growth_plot_all = growth_plot_all + geom_point(aes(colour=treatment, shape=treatment), size=4)
growth_plot_all = growth_plot_all + labs(x="Treatments", y = "Surface growth area (%)")
growth_plot_all = growth_plot_all + theme_bw()+ theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
growth_plot_all = growth_plot_all + scale_shape_manual(values=c(15,19,17,18))
growth_plot_all = growth_plot_all + theme(legend.title = element_blank(), legend.position = "none")+
  theme(axis.text=element_text(size=16),axis.title = element_text(size=18),strip.text=element_text(size=16))
growth_plot_all = growth_plot_all + scale_colour_manual(values = wes_palette("Darjeeling1", n = 4))
growth_plot_all = growth_plot_all + facet_wrap(~temperature)
growth_plot_all

#ggsave("../plots/ageing_growth.png", dpi=300, device="png", width=9, height=8)




