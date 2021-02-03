library(dplyr)
library(tidyr)

#Reading the raw FTIR data
spectra_all = read.csv("raw_data/ftir_raw.csv")

#Cleaning up the dataset
spectra_all$X= NULL
spectra_all = spectra_all %>%
  dplyr::rename("wavenumber" = V1) %>%
  dplyr::rename("height" = V2)
spectra_all$temperature= as.factor(spectra_all$temperature)


######### For Fig. 1b ##################
library(ggplot2)
library(wesanderson)

#Filtering the "350" temperature dataset and adjusting the distnace between the different spectra
spectra_350 = spectra_all %>%
  dplyr::filter(temperature == 350)

spectra_350$wavenumber=as.numeric(paste(spectra_350$wavenumber))
spectra_350$height=as.numeric(paste(spectra_350$height))
spectra_350[spectra_350$treatment=="Physical",]$height=spectra_350[spectra_350$treatment=="Physical",]$height+2
spectra_350[spectra_350$treatment=="Chemical",]$height=spectra_350[spectra_350$treatment=="Chemical",]$height+4
spectra_350[spectra_350$treatment=="Biological",]$height=spectra_350[spectra_350$treatment=="Biological",]$height+6

#Filtering the "550" temperature dataset and adjusting the distnace between the different spectra
spectra_550 = spectra_all %>%
  dplyr::filter(temperature == 550)

spectra_550$wavenumber=as.numeric(paste(spectra_550$wavenumber))
spectra_550$height=as.numeric(paste(spectra_550$height))
spectra_550[spectra_550$treatment=="Physical",]$height=spectra_550[spectra_550$treatment=="Physical",]$height+2
spectra_550[spectra_550$treatment=="Chemical",]$height=spectra_550[spectra_550$treatment=="Chemical",]$height+4
spectra_550[spectra_550$treatment=="Biological",]$height=spectra_550[spectra_550$treatment=="Biological",]$height+6

#Plotting
plot_all = ggplot(NULL, mapping = aes())
plot_all = plot_all + geom_line(data=spectra_550, aes(x=wavenumber,y=height, color=treatment), size=1)
plot_all = plot_all + geom_line(data=spectra_350, aes(x=wavenumber,y=height, color=treatment), size=1)
plot_all = plot_all + scale_x_reverse()
plot_all = plot_all + labs(x = expression(paste("Wavenumber ", (cm^-1))))
plot_all = plot_all + theme_bw() + theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(), 
                                         axis.line = element_line(colour = "black"),
                                         legend.position = "none",
                                         axis.text=element_text(size=18),
                                         axis.title = element_text(size=20), 
                                         strip.text = element_text(size= 20),
                                         axis.title.y =element_blank(),
                                         axis.text.y=element_blank(),
                                         axis.ticks.y=element_blank())
plot_all = plot_all + scale_colour_manual(values = wes_palette("Darjeeling1", n = 4))
plot_all = plot_all + facet_wrap(~temperature)
#annotate treatments
plot_all = plot_all + annotate(geom="text",x = 2400, y=0.9, label = "Unaged", colour = "black", size = 8)
plot_all = plot_all + annotate(geom="text",x = 2400, y=2.6, label = "Physical", colour = "black", size = 8)
plot_all = plot_all + annotate(geom="text",x = 2400, y=5, label = "Chemical",  colour = "black", size = 8)
plot_all = plot_all + annotate(geom="text",x = 2400, y=7, label = "Biological", colour = "black", size = 8)
# 3300 OH
plot_all = plot_all + geom_segment(aes(x = 3350, y=0, xend = 3350, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + geom_segment(aes(x = 3350, y=0, xend = 3350, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + annotate(geom="text",x = 3350, y=8.7, label = "O-H", colour = "black", size = 6)
# 2900 CH
plot_all = plot_all + geom_segment(aes(x = 2932, y=0, xend = 2932, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + geom_segment(aes(x = 2932, y=0, xend = 2932, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + annotate(geom="text",x = 2932, y=8.9, label = "al", colour = "black", size = 6)
plot_all = plot_all + annotate(geom="text",x = 2932, y=8.7, label = "C-H", colour = "black", size = 6)
# 1700 CO
plot_all = plot_all + geom_segment(aes(x = 1700, y=0, xend = 1700, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + geom_segment(aes(x = 1700, y=0, xend = 1700, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + annotate(geom="text",x = 1775, y=8.7, label = "C=O", colour = "black", size= 6)
#1600 CC
plot_all = plot_all + geom_segment(aes(x = 1590, y=0, xend = 1590, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + geom_segment(aes(x = 1600, y=0, xend = 1600, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + annotate(geom="text",x = 1580, y=8.9, label = "ar", colour = "black", size = 6)
plot_all = plot_all + annotate(geom="text",x = 1580, y=8.7, label = "C=C", colour = "black", size = 6)
#1400 CH
plot_all = plot_all + geom_segment(aes(x = 1435, y=0, xend = 1435, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + geom_segment(aes(x = 1415, y=0, xend = 1415, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + annotate(geom="text",x = 1445, y=8.9, label = "al", colour = "black", size = 6)
plot_all = plot_all + annotate(geom="text",x = 1445, y=8.7, label = "C-H", colour = "black", size = 6)
#1200 CO*
plot_all = plot_all + geom_segment(aes(x = 1200, y=0, xend = 1200, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + geom_segment(aes(x = 1200, y=0, xend = 1200, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + annotate(geom="text",x = 1200, y=8.7, label = "CO*", colour = "black", size = 6)
#800 CH
plot_all = plot_all + geom_segment(aes(x = 750, y=0, xend = 750, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + geom_segment(aes(x = 750, y=0, xend = 750, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + geom_segment(aes(x = 900, y=0, xend = 900, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + geom_segment(aes(x = 900, y=0, xend = 900, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.0 )
plot_all = plot_all + annotate(geom="text",x = 820, y=8.9, label = "ar", colour = "black", size = 6)
plot_all = plot_all + annotate(geom="text",x = 820, y=8.7, label = "C-H", colour = "black", size = 6)
plot_all

######## For Fig. 2b ####################
library(dendextend)
library(viridis)

#Creating another column to combine temp. and treatment information
spectra_all$sample_name = paste(spectra_all$temperature, spectra_all$treatment, sep=" ")
spectra_all$temperature= NULL
spectra_all$treatment = NULL

#Rounding up the wavenumbers and removing data past the 3000 cm-1 wavenumber to avoid interference from water
spectra_all$wavenumber = as.numeric(as.character(spectra_all$wavenumber))
spectra_all$wavenumber = round(spectra_all$wavenumber)
spectra_all = spectra_all %>%
  filter(wavenumber <=3000)

#Pivot the dataframe for dendrogram analysis
spectra_all_wide = pivot_wider(spectra_all, names_from = wavenumber, values_from = height, values_fn = list(height = mean))

#pivot creates a tibble, not a dataframe
spectra_all_wide = as.data.frame(spectra_all_wide)
rownames(spectra_all_wide) = spectra_all_wide$sample_name
spectra_all_wide$sample_name = NULL

#Plotting and saving the dendrogram
png(file="../plots/ftir_norm_3000.png",
    width=600, height=350)

par(mar =  c(4, 4, 4, 6) + 0.1)

dend <- as.dendrogram(hclust(dist(spectra_all_wide)))
colors_to_use = as.numeric(as.factor(substr(rownames(spectra_all_wide), 0, 3))) #colour with temperature, add +1,2.... to change colours
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend) <- colors_to_use

dend %>% highlight_branches_col %>% plot(horiz= TRUE, xlab = "distance between biochar samples", cex.lab = 0.9)


dend %>% rect.dendrogram(k=3, horiz = TRUE,
                         border = 8, lty = 5, lwd = 2)

dev.off()
