#Authors: Nayela Zeba, Thea L. Whitman

#Load necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(dendextend)
library(viridis)

#Reading the raw FTIR data
spectra_all = read.csv("raw_data/ftir_raw.csv")

#Cleaning up the dataset
spectra_all$X= NULL
spectra_all = spectra_all %>%
  dplyr::rename("wavenumber" = V1) %>%
  dplyr::rename("height" = V2)
spectra_all$temperature= as.factor(spectra_all$temperature)


######### For Fig. 1a ##################


#Filtering the "350" temperature dataset and adjusting the distance between the different spectra
spectra_350 = spectra_all %>%
  dplyr::filter(temperature == 350)

spectra_350$wavenumber=as.numeric(paste(spectra_350$wavenumber))
spectra_350$height=as.numeric(paste(spectra_350$height))
spectra_350[spectra_350$treatment=="Physical",]$height=spectra_350[spectra_350$treatment=="Physical",]$height+2
spectra_350[spectra_350$treatment=="Chemical",]$height=spectra_350[spectra_350$treatment=="Chemical",]$height+4
spectra_350[spectra_350$treatment=="Biological",]$height=spectra_350[spectra_350$treatment=="Biological",]$height+6

#Filtering the "550" temperature dataset and adjusting the distance between the different spectra
spectra_550 = spectra_all %>%
  dplyr::filter(temperature == 550)

spectra_550$wavenumber=as.numeric(paste(spectra_550$wavenumber))
spectra_550$height=as.numeric(paste(spectra_550$height))
spectra_550[spectra_550$treatment=="Physical",]$height=spectra_550[spectra_550$treatment=="Physical",]$height+2
spectra_550[spectra_550$treatment=="Chemical",]$height=spectra_550[spectra_550$treatment=="Chemical",]$height+4
spectra_550[spectra_550$treatment=="Biological",]$height=spectra_550[spectra_550$treatment=="Biological",]$height+6

#Plotting
plot_all = ggplot(NULL, mapping = aes())
plot_all = plot_all + geom_line(data=spectra_550, aes(x=wavenumber,y=height, color=treatment), size=1.5)
plot_all = plot_all + geom_line(data=spectra_350, aes(x=wavenumber,y=height, color=treatment), size=1.5)
plot_all = plot_all + scale_x_reverse()
plot_all = plot_all + labs(x = expression(paste("Wavenumber ", (cm^-1))))
plot_all = plot_all + theme_bw() + theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(), 
                                         axis.line.y = element_line(colour = "black", size = 1.5),
                                         axis.line.x = element_line(colour = "black", size=1.5),
                                         legend.position = "none",
                                         axis.text=element_text(size=26),
                                         axis.title = element_text(size=26), 
                                         strip.text = element_text(size= 26),
                                         axis.title.y =element_blank(),
                                         axis.text.y=element_blank(),
                                         axis.ticks.y=element_blank(),
                                         panel.border = element_rect(colour = "black", size = 1.5))
plot_all = plot_all + scale_colour_manual(values = wes_palette("Darjeeling1", n = 4))
plot_all = plot_all + facet_wrap(~temperature)
#annotate treatments
plot_all = plot_all + annotate(geom="text",x = 2400, y=0.9, label = "Unaged", colour = "black", size = 8)
plot_all = plot_all + annotate(geom="text",x = 2400, y=2.6, label = "Physical", colour = "black", size = 8)
plot_all = plot_all + annotate(geom="text",x = 2400, y=5, label = "Chemical",  colour = "black", size = 8)
plot_all = plot_all + annotate(geom="text",x = 2400, y=7, label = "Biological", colour = "black", size = 8)
# 3300 OH
plot_all = plot_all + geom_segment(aes(x = 3350, y=0, xend = 3350, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + geom_segment(aes(x = 3350, y=0, xend = 3350, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + annotate(geom="text",x = 3350, y=8.7, label = "O-H", colour = "black", size = 8)
# 2900 CH
plot_all = plot_all + geom_segment(aes(x = 2932, y=0, xend = 2932, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + geom_segment(aes(x = 2932, y=0, xend = 2932, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + annotate(geom="text",x = 2932, y=8.9, label = "al", colour = "black", size = 8)
plot_all = plot_all + annotate(geom="text",x = 2932, y=8.7, label = "C-H", colour = "black", size = 8)
# 1700 CO
plot_all = plot_all + geom_segment(aes(x = 1700, y=0, xend = 1700, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + geom_segment(aes(x = 1700, y=0, xend = 1700, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + annotate(geom="text",x = 1785, y=8.7, label = "C=O", colour = "black", size= 8)
#1600 CC
plot_all = plot_all + geom_segment(aes(x = 1590, y=0, xend = 1590, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + geom_segment(aes(x = 1600, y=0, xend = 1600, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + annotate(geom="text",x = 1580, y=8.9, label = "ar", colour = "black", size = 8)
plot_all = plot_all + annotate(geom="text",x = 1580, y=8.7, label = "C=C", colour = "black", size = 8)
#1400 CH
plot_all = plot_all + geom_segment(aes(x = 1435, y=0, xend = 1435, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + geom_segment(aes(x = 1415, y=0, xend = 1415, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + annotate(geom="text",x = 1395, y=8.9, label = "al", colour = "black", size = 8)
plot_all = plot_all + annotate(geom="text",x = 1395, y=8.7, label = "C-H", colour = "black", size = 8)
#1200 CO*
plot_all = plot_all + geom_segment(aes(x = 1200, y=0, xend = 1200, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + geom_segment(aes(x = 1200, y=0, xend = 1200, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + annotate(geom="text",x = 1180, y=8.7, label = "CO*", colour = "black", size = 8)
#800 CH
plot_all = plot_all + geom_segment(aes(x = 750, y=0, xend = 750, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + geom_segment(aes(x = 750, y=0, xend = 750, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + geom_segment(aes(x = 900, y=0, xend = 900, yend = 8.5), data = spectra_550, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + geom_segment(aes(x = 900, y=0, xend = 900, yend = 8.5), data = spectra_350, linetype="dotted", color = "black", size=1.5 )
plot_all = plot_all + annotate(geom="text",x = 820, y=8.9, label = "ar", colour = "black", size = 8)
plot_all = plot_all + annotate(geom="text",x = 820, y=8.7, label = "C-H", colour = "black", size = 8)
plot_all

ggsave("../plots/Rplot_FTIR.png", dpi = 300, device="png", width=24, height=14)

######## For Fig. 1b ####################


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
    width=300, height=175, res=300, units="mm")

par(mar =  c(4, 4, 4, 6) + 0.1, cex=1.5)

dend <- as.dendrogram(hclust(dist(spectra_all_wide)))
colors_to_use = as.numeric(as.factor(substr(rownames(spectra_all_wide), 0, 3))) #colour with temperature, add +1,2.... to change colours
colors_to_use <- colors_to_use[order.dendrogram(dend)]
labels_colors(dend) <- colors_to_use

#dend %>% highlight_branches_col %>% plot(horiz= TRUE, xlab = "distance between biochar samples", lwd= 4)
dend %>%
  highlight_branches() %>%
  plot(main = "Emphasizing color\n and line-width", horiz= TRUE, xlab = "Distance between biochar samples")


dend %>% rect.dendrogram(k=3, horiz = TRUE,
                         border = 8, lty = 5, lwd = 4)

dev.off()

