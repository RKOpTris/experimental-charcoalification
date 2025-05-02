library(stringr)
library(dplyr)
library(ggplot2)
library(Ternary)

###################
######### MAIN TEXT
###################

### Need to delete this
setwd("/Users/mukappa/Documents/Teratologies/Charred lycopodium")
source("Charcoalified Lycopodium analysis FUNCTIONS.R", local = T)

# Load in data
md <- read.csv("Charcoalified Lycopodium METADATA.csv")
sp <- read.csv("Charcoalified Lycopodium SPECTRA.csv")
pdi <- read.csv("Charcoalified Lycopodium PDI.csv")
lyc <- read.csv("Charcoalified Lycopodium SIZE.csv")
wl <- read.csv("Charcoalified Lycopodium MASS LOSS.csv")
rouxhet <- read.csv("Charcoalified Lycopodium Rouxhet chemical data TABLE.csv")

# Remove ambient air CO2 spikes from 2410-2200 cm-1
sp <- remove.co2.signal(sp)

# Convert data to long format and calculate sample-wise statistics
sp$sample <- md$sample
sp_long <- sp %>% 
  tidyr::gather("wavenumber", "absorbance", -sample) %>%
  mutate(wavenumber = as.numeric(str_replace_all(wavenumber, "[^0-9.-]", "")))
sp_means <- sp_long %>% group_by(sample, wavenumber) %>% dplyr::summarise(mean = mean(absorbance),
                                                                          sd = sd(absorbance),
                                                                          upper = mean + sd,
                                                                          lower = mean - sd)

# Plot preferences
publication_labels <- paste(strip_non_numbers(unique(sp_means$sample)), "\u00b0C")
publication_labels[1] <- "Unpyrolysed"
pointCol <- "#000000BB"
smoothCol <- scico::scico(2, palette = "roma", alpha = 0.3)[2]
overplot_palette <- c(scico::scico(length(unique(sp_means$sample)) + 6, palette = "roma")[rev(c(1:4, 12:17))], "black")


sp_peaks_orig <- sp_means %>% 
  select(-c(sd, upper, lower)) %>% 
  tidyr::spread("wavenumber", "mean")
sp_peaks <- sp_peaks_orig
wavenumbers <- strip_non_numbers(names(sp_peaks_orig[-1]))
sp_peaks_upper <- sp_means %>% 
  select(-c(mean, lower, sd)) %>% 
  tidyr::spread("wavenumber", "upper")
sp_peaks_lower <- sp_means %>% 
  select(-c(mean, sd, upper)) %>% 
  tidyr::spread("wavenumber", "lower")

for(i in 1:nrow(sp_peaks)){
  sp_peaks[i, -1] <- sp_peaks[i, -1] + ((i - 1) * 0.04)
  sp_peaks_upper[i, -1] <- sp_peaks_upper[i, -1] + ((i - 1) * 0.04)
  sp_peaks_lower[i, -1] <- sp_peaks_lower[i, -1] + ((i - 1) * 0.04)
}


######## MAIN TEXT FIGURE 1
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 1.1, 1.1), pty = "m")
plot(wavenumbers, as.numeric(sp_peaks[1, -1]), 
     xlim = c(3600, 900),
     ylim = range(sp_peaks[-1]),
     type = "n", axes = F,
     xlab = expression(paste("Wavenumber (cm"^"-1"*")")),
     ylab = "",
     cex.lab = 2)

arrows(3500, pull(sp_peaks[, ncol(sp_peaks)]), 
       950, pull(sp_peaks[, ncol(sp_peaks)]), col = "lightgrey",
       length = 0)

for(i in 1:nrow(sp_peaks)){
  polygon(c(wavenumbers, rev(wavenumbers)), 
          c(as.numeric(sp_peaks_upper[i, -1]), rev(as.numeric(sp_peaks_lower[i, -1]))), col = smoothCol, border = NA)
  lines(wavenumbers, as.numeric(sp_peaks[i, -1]), col = pointCol)
  spec_peaks <- find_peaks(as.numeric(sp_peaks[i, -1]), type = "peak", nups = 5, nudge = 0.003)
  spec_peaks <- spec_peaks %>% filter(peak_wavenumber < 3050, peak_wavenumber != 1654.6)
  spec_peaks <- spec_peaks[!(spec_peaks$peak_wavenumber < 2800 & spec_peaks$peak_wavenumber > 1900), ]
  text(nudge ~ peak_wavenumber, spec_peaks, labels = round(peak_wavenumber, 0), cex = 0.6)
}
text(3600, pull(sp_peaks[1585]) - 0.0075, labels = publication_labels, adj = 0, col = "black", cex = 1.5)
internal_axis_ticks(1, 1, cex.axis = 1.5)
mtext(2, 2, at = inRange(as.matrix(sp_peaks[-1]), 0.5), text = "Absorbance", cex = 2)
state_x <- 875
text(state_x, sp_peaks[1, 2], "I", cex = 2)
text(state_x, sp_peaks[1, 2], "O", cex = 5, col = pointCol)
text(state_x, sp_peaks[7, 2], "II", cex = 2)
text(state_x, sp_peaks[7, 2], "O", cex = 5, col = pointCol)
text(state_x, sp_peaks[11, 2], "III", cex = 2)
text(state_x, sp_peaks[11, 2], "O", cex = 5, col = pointCol)
arrows(state_x, as.numeric(sp_peaks[7, 2]) + 0.015,
       state_x, as.numeric(sp_peaks[10, 2]), angle = 90, lwd = 2, length = 0.15)
arrows(state_x, as.numeric(sp_peaks[1, 2]) + 0.015, 
       state_x, as.numeric(sp_peaks[2, 2]), angle = 90, lwd = 2, length = 0, lty = 3)
arrows(state_x, as.numeric(sp_peaks[2, 2]), 
       state_x, as.numeric(sp_peaks[6, 2]), angle = 90, lwd = 2, length = 0.15)


######## MAIN TEXT FIGURE 2a
sp_basic <- sp %>% group_by(sample) %>% summarise_all("mean")

sp_pca <- prcomp(sp_basic[-1], center = T, scale = F)
summary(sp_pca)
par(mfrow = c(1, 1), mar = c(4.1, 5.1, 1.1, 1.1), pty = "s")
plot(sp_pca$x[, 1], sp_pca$x[, 2], type = "n", axes = F, xlab = "PC1 score (82.2%)", ylab = "PC2 score (11.7%)", cex.lab = 1.5, cex.axis = 1.1,
     xlim = c(inRange(sp_pca$x[, 1], c(-0.05, 1.05))),
     ylim = c(inRange(sp_pca$x[, 2], c(-0.05, 1.05))))
abline(h = 0, col = "darkgrey", lty = 1)
abline(v = 0, col = "darkgrey", lty = 1)
text(sp_pca$x[, 1], sp_pca$x[, 2], labels = strip_non_numbers(sp_peaks$sample), col = overplot_palette, cex = 1.8)
internal_axis_ticks(1:2, 1:2, cex.axis = 1.3)
plotrix::draw.ellipse(x = -0.075, y = 0, a = 0.035, b = 0.035, angle = 30, border = pointCol, lwd = 2)
plotrix::draw.ellipse(x = 0.080, y = -0.024, a = 0.023, b = 0.065, angle = -80, border = pointCol, lwd = 2)
plotrix::draw.ellipse(x = 0.113, y = 0.085, a = 0.015, b = 0.015, angle = -80, border = pointCol, lwd = 2)
text(-0.075, -0.02, "I", cex = 2)
text(-0.075, -0.02, "O", cex = 5, col = pointCol)
text(x = 0.075, y = -0.018, "II", cex = 2)
text(x = 0.075, y = -0.018, "O", cex = 5, col = pointCol)
text(x = 0.113, y = 0.08, "III", cex = 2)
text(x = 0.113, y = 0.08, "O", cex = 5, col = pointCol)
mtext(2, 4, at = 0.1, text = "a)", las = 1, cex = 2)

######## MAIN TEXT FIGURE 2b
par(mar = c(4.1, 6.1, 1.1, 1.1), pty = "m")
plot(wavenumbers, sp_pca$rotation[, 1], 
     xlim = c(3600, 950), type = "l", axes = F,
     ylim = inRange(sp_pca$rotation[, 1], c(-0.1, 1.1)),
     ylab = "PC1 loading (82.2%)",
     xlab = expression(paste("Wavenumber (cm"^"-1"*")")),
     cex.lab = 1.5)
abline(h = 0, col = "darkgrey", lty = 1)
internal_axis_ticks(c(1, 2), c(1, 2), cex.axis = 1.3)
pc1_peaks <- find_peaks(sp_pca$rotation[, 1], type = "both", nups = 6, ndowns = 6, nudge = 0.01)
pc1_peaks <- pc1_peaks %>% filter(!peak_wavenumber %in% c(2879.2, 2740.4))
text(nudge ~ peak_wavenumber, pc1_peaks, labels = round(peak_wavenumber, 0), cex = 0.9)
mtext(2, 4, at = max(pc1_peaks$nudge), text = "b)", las = 1, cex = 2)

######## MAIN TEXT FIGURE 2c
par(mar = c(4.1, 6.1, 1.1, 1.1), pty = "m")
plot(wavenumbers, sp_pca$rotation[, 2], 
     xlim = c(3600, 950), type = "l", axes = F,
     ylim = inRange(sp_pca$rotation[, 2], c(-0.1, 1.1)),
     ylab = "PC2 loading (11.7%)",
     xlab = expression(paste("Wavenumber (cm"^"-1"*")")),
     cex.lab = 1.5)
abline(h = 0, col = "darkgrey", lty = 1)
internal_axis_ticks(1:2, 1:2, cex.axis = 1.3)
pc2_peaks <- find_peaks(sp_pca$rotation[, 2], type = "both", nups = 6, ndowns = 6, nudge = 0.005)
pc2_peaks <- pc2_peaks %>% filter(peak_wavenumber < 3010, !peak_wavenumber %in% c(2973.7, 2946.7, 2883.1, 2838.7, 2690.2))
text(nudge ~ peak_wavenumber, pc2_peaks, labels = round(peak_wavenumber, 0), cex = 0.9)
mtext(2, 4, at = max(pc2_peaks$nudge), text = "c)", las = 1, cex = 2)


######## MAIN TEXT FIGURE 3a
sp_means2 <- sp_means
sp_means2$temp <- strip_non_numbers(sp_means2$sample)
sp_means2 <- sp_means2 %>% filter(temp > 225 & temp < 375)
sp_means3 <- sp_means2
for(i in 1:length(unique(sp_means3$temp))){
  sp_means3[sp_means3$temp == unique(sp_means3$temp)[i], ]$mean <- sp_means3[sp_means3$temp == unique(sp_means3$temp)[i], ]$mean - sp_means2[sp_means2$temp == 250, ]$mean
}
sp_means3 <- sp_means3 %>% filter(wavenumber < 3600)
sp_means4 <- sp_means
sp_means4$temp <- strip_non_numbers(sp_means4$sample)
sp_means4 <- sp_means4 %>% filter(temp > 325)
sp_means5 <- sp_means4
for(i in 1:length(unique(sp_means5$temp))){
  sp_means5[sp_means5$temp == unique(sp_means5$temp)[i], ]$mean <- sp_means5[sp_means5$temp == unique(sp_means5$temp)[i], ]$mean - sp_means4[sp_means4$temp == 350, ]$mean
}
sp_means5 <- sp_means5 %>% filter(wavenumber < 3600)
sp_means6 <- sp_means
sp_means6$temp <- strip_non_numbers(sp_means6$sample)
sp_means6 <- sp_means6 %>% filter(temp < 275)
sp_means7 <- sp_means6
for(i in 1:length(unique(sp_means7$temp))){
  sp_means7[sp_means7$temp == unique(sp_means7$temp)[i], ]$mean <- sp_means7[sp_means7$temp == unique(sp_means7$temp)[i], ]$mean - sp_means6[sp_means6$temp == 0, ]$mean
}
sp_means7 <- sp_means7 %>% filter(wavenumber < 3600)
sp_means_all_range <- c(sp_means3$mean, sp_means5$mean, sp_means7$mean) %>% range()
sp_means_all_range <- inRange(sp_means_all_range, c(-0.05, 1.2))

absorbance_bands <- data.frame(xmin = c(2840, 1510, 1680, 1407, 1068, 1225, 1300, 3100),  
                               ymin = -1, 
                               xmax = c(3010, 1630, 1775, 1473, 1159, 1294, 1367, 3450),
                               ymax = 1
)

overplot_palette <- c(scico::scico(length(unique(sp_means$sample)) + 6, palette = "roma")[rev(c(1:4, 12:17))], "black")

sp_means7_pn_max <- sp_means7 %>% 
  group_by(wavenumber) %>% 
  summarise(max = pn_max(mean))
sp_means7_peaks <- find_peaks(sp_means7_pn_max$max, sp_means7_pn_max$wavenumber, type = "both", nups = 5, ndowns = 5, nudge = 0.001)  %>% 
  filter(peak_wavenumber < 3008) %>%
  filter(!(peak_wavenumber < 2800 & peak_wavenumber > 1800)) %>%
  filter((peak_type == "trough" & abs < 0) | (peak_type == "peak" & abs > 0))

mfrow(1, 1)
mar(4.1, 6.1, 1.1, 1.1)
plot(mean ~ wavenumber, sp_means7, xlim = c(3600, 950), 
     ylim = sp_means_all_range, 
     type = "n", 
     axes = F, 
     xlab = expression(paste("Wavenumber (cm"^"-1"*")")),
     ylab = "",
     cex.lab = 1.5)
draw_bands()
labels_for_bands(max(sp_means_all_range))
abline(h = 0, col = "darkgrey")
for(i in 2:length(unique(sp_means7$temp))){
  lines(mean ~ wavenumber, sp_means7[sp_means7$temp == unique(sp_means7$temp)[i], ], 
        col = overplot_palette[i],
        lwd = 1)
}
mtext(2, 3.3, at = inRange(sp_means_all_range), text = "Absorbance relative to unpyrolysed\n(arbitrary units)", cex = 1.5)
legend("topleft", lty = 1, lwd = 2, cex = 1, col = overplot_palette[2:6], legend = publication_labels[2:6])
mtext(2, 4, at = 0.046, text = "a)", las = 1, cex = 2)
internal_axis_ticks(1:4, 1:2)

######## MAIN TEXT FIGURE 3b
sp_means3_pn_max <- sp_means3 %>% 
  group_by(wavenumber) %>% 
  summarise(max = pn_max(mean))
sp_means3_peaks <- find_peaks(sp_means3_pn_max$max, sp_means3_pn_max$wavenumber, type = "both", nups = 5, ndowns = 5, nudge = 0.001)  %>% 
  filter(peak_wavenumber < 3008) %>%
  filter(!(peak_wavenumber < 2800 & peak_wavenumber > 1800)) %>%
  filter((peak_type == "trough" & abs < 0) | (peak_type == "peak" & abs > 0)) %>%
  filter(peak_wavenumber != 2989.1)

mar(4.1, 6.1, 1.1, 1.1)
plot(mean ~ wavenumber, sp_means3, xlim = c(3600, 950), type = "n", 
     ylim = sp_means_all_range,
     axes = F, 
     xlab = expression(paste("Wavenumber (cm"^"-1"*")")),
     ylab = "",
     cex.lab = 1.5)
draw_bands()
labels_for_bands(max(sp_means_all_range))
abline(h = 0, col = "darkgrey")
for(i in 2:length(unique(sp_means3$temp))){
  lines(mean ~ wavenumber, sp_means3[sp_means3$temp == unique(sp_means3$temp)[i], ], 
        col = overplot_palette[6:13][i],
        lwd = 1)
}
mtext(2, 3.3, at = inRange(sp_means_all_range), text = "Absorbance relative to 250 \u00b0C\n(arbitrary units)", cex = 1.5)
legend("topleft", lty = 1, lwd = 2, cex = 1, col = overplot_palette[7:10], legend = publication_labels[7:10])
internal_axis_ticks(1:4, 1:2)
mtext(2, 4, at = 0.046, text = "b)", las = 1, cex = 2)

######## MAIN TEXT FIGURE 3c
sp_means5_pn_max <- sp_means5 %>% 
  group_by(wavenumber) %>% 
  summarise(max = pn_max(mean))
sp_means5_peaks <- find_peaks(sp_means5_pn_max$max, sp_means5_pn_max$wavenumber, type = "both", nups = 5, ndowns = 5, nudge = 0.001)  %>% 
  filter(peak_wavenumber < 3008) %>%
  filter(!(peak_wavenumber < 2800 & peak_wavenumber > 1800)) %>%
  filter((peak_type == "trough" & abs < 0) | (peak_type == "peak" & abs > 0)) %>%
  filter(!peak_wavenumber %in% c(1060.7, 1309.4, 985.4, 2838.7, 2883.1, 2944.8, 2971.8))

mar(4.1, 6.1, 1.1, 1.1)
plot(mean ~ wavenumber, sp_means5, xlim = c(3600, 950), type = "n", 
     axes = F, 
     xlab = expression(paste("Wavenumber (cm"^"-1"*")")),
     ylim = sp_means_all_range,
     ylab = "",
     cex.lab = 1.5)
draw_bands()
labels_for_bands(max(sp_means_all_range))
abline(h = 0, col = "darkgrey")
for(i in 2:length(unique(sp_means5$temp))){
  lines(mean ~ wavenumber, sp_means5[sp_means5$temp == unique(sp_means5$temp)[i], ], 
        col = overplot_palette[10:11][i],
        lwd = 1)
}
legend("topleft", lty = 1, lwd = 2, cex = 1, col = overplot_palette[11], legend = publication_labels[11])
internal_axis_ticks(1:4, 1:2)
mtext(2, 3.3, at = inRange(sp_means_all_range), text = "Absorbance relative to 350 \u00b0C\n(arbitrary units)", cex = 1.5)
mtext(2, 4, at = 0.046, text = "c)", las = 1, cex = 2)


######## MAIN TEXT FIGURE 4
pdi_short <- pdi %>% filter(Taxon == "Lycopodium", Temp < 1000, !Temp %in% c(100, 125, 400)) %>% arrange(Temp)
pdi_short$plot_col <- scico::scico(length(unique(pdi_short$Temp)), palette = "lajolla")[as.factor(pdi_short$Temp)]
pdi_short <- pdi_short %>% mutate(actual_col = paste0("#", as.character(as.hexmode(R)), as.character(as.hexmode(G)), str_pad(as.character(as.hexmode(B)), pad = "0", width = 2)))

pdi_pca <- pdi_short %>% dplyr::select(R, G, B) %>% prcomp(scale = T, center = T)
centroids <- pdi_pca$x %>% 
  data.frame() %>% 
  mutate(Temp = pdi_short$Temp) %>% 
  group_by(Temp) %>% 
  summarise_all("mean")
mean_cols <- pdi_short[, c("Temp", "R", "G", "B")] %>% group_by(Temp) %>% dplyr::summarise_all("mean")
mean_cols[2:4] <- lapply(mean_cols[2:4], function(x){as.hexmode(round(x, 0))})
mean_cols <- mean_cols %>% mutate(RGB_col = paste0("#", R, G, B))
centroids <- left_join(centroids, dplyr::select(mean_cols, c(Temp, RGB_col)), by = "Temp")

mar(5.1, 5.1, 1.1, 1.1)
plot(pdi_pca$x[, 1], pdi_pca$x[, 2], type = "n", axes = FALSE,
     xlab = "PC1 score (89.5%)", ylab = "PC2 score (9.8%)", cex.lab = 2, las = 1)
abline(h = 0, v = 0, col = "darkgrey")
points(pdi_pca$x[, 1], pdi_pca$x[, 2], col = pdi_short$actual_col, pch = ".", cex = 3)
internal_axis_ticks(1:2, 1:2, cex.axis = 1.5)
text(centroids$PC1, centroids$PC2, labels = centroids$Temp, col = centroids$RGB_col, cex = 2.5)

points(rep(0.25, 13), rev(seq(0.1, max(pdi_pca$x[, 2]), length.out = 13)), pch = 16, cex = 3, col = centroids$RGB_col)
text(rep(0.5, 13), rev(seq(0.1, max(pdi_pca$x[, 2]), length.out = 13)), c("Non-pyrolysed (0 \u00b0C)", paste0(centroids$Temp[-1], " \u00b0C")), adj = 0)


######## MAIN TEXT FIGURE 5a
pdi %>% filter(Taxon == "Lycopodium") %>% 
  dplyr::select(-c(Sample, Spore_rep, Taxon)) %>% 
  tidyr::gather("Band", "Value", -Temp) %>%
  summarySE(measurevar = "Value", groupvars = c("Temp", "Band")) %>%
  filter(Band == "PDI") %>%
  ggplot(aes(Temp, Value), color = "black") + 
  stat_smooth(method = "loess", col = smoothCol, fill = smoothCol, alpha = 0.2) +
  geom_errorbar(aes(ymin = Value - sd, ymax = Value + sd), width = 0, col = pointCol) +
  geom_point(aes(Temp, Value), size = 2, col = pointCol) +
  scale_x_continuous(breaks = seq(0, 800, by = 100)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_bw() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "top",
        title = element_text(size = 25)) +
  labs(x = "Pyrolysis temperature (\u00b0C)", y = "Palynomorph Darkness Index (%)", title = "a)")

######## MAIN TEXT FIGURE 5b
pdi %>% 
  filter(Taxon == "Lycopodium") %>% 
  dplyr::select(-c(Sample, Spore_rep, Taxon)) %>% 
  ggplot() + 
  geom_density(aes(PDI, fill = as.factor(Temp)), alpha = 0.8) +
  scale_fill_manual(values = centroids$RGB_col, name = "Pyrolysis\ntemperature\n(\u00b0C)") +
  theme_bw() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = "right",
        title = element_text(size = 25)) +
  labs(x = "Palynomorph Darkness Index (%)", y = "Density", title = "b)")

######## MAIN TEXT FIGURE 6a
PDIvals <- pdi %>% group_by(Temp) %>% dplyr::summarise(PDI_sd = sd(PDI),
                                                       PDI = mean(PDI))

lycPDI <- lyc %>% 
  group_by(Temperature) %>% 
  dplyr::summarise(Length_Mean = mean(Length),
                   Length_SD = sd(Length)) %>%
  inner_join(PDIvals, by = c("Temperature" = "Temp")) 

lycPDI %>%
  ggplot(aes(Temperature, Length_Mean)) + 
  geom_point(size = 2, col = pointCol) +
  geom_errorbar(aes(ymin = Length_Mean - Length_SD,
                    ymax = Length_Mean + Length_SD),
                width = 0, col = pointCol) +
  stat_smooth(method = "lm", col = smoothCol, fill = smoothCol, alpha = 0.2) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 40, 5), limits = c(14, 40)) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20),
        title = element_text(size = 25)) +
  labs(x = paste0("Pyrolysis temperature (", "\U00B0", "C)"),
       y = paste0("Spore length (", "\U03BC", "m)"),
       title = "a)") +
  annotate("text", x = 0, y = 14, label = "Spore length = −0.015 × Pyrolysis temperature + 33.15", hjust = 0, size = 5)

robust_model_6a <- robustbase::lmrob(Length_Mean ~ Temperature, data = lycPDI)
summary(robust_model_6a)


### Estimate of N required in order to obtain more accurate mean
meanSD <- lyc %>% 
  group_by(Temperature) %>% 
  dplyr::summarise(sd = sd(Length, na.rm = T)) %>% 
  pull(sd) %>% 
  mean %>% 
  round(1)
samplingbook::sample.size.mean(e = 0.5, S = meanSD, N = Inf, level = 0.95)

######## MAIN TEXT FIGURE 6b
lycPDI %>%
  ggplot(aes(Length_Mean, PDI)) + 
  geom_point(size = 2, col = pointCol) +
  geom_errorbarh(aes(xmin = Length_Mean - Length_SD,
                     xmax = Length_Mean + Length_SD),
                 height = 0, col = pointCol) +
  geom_errorbar(aes(ymin = PDI - PDI_sd,
                    ymax = PDI + PDI_sd),
                width = 0, col = pointCol) +
  stat_smooth(method = "lm", col = smoothCol, fill = smoothCol, alpha = 0.2) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 40, 5), limits = c(14, 40)) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20),
        title = element_text(size = 25)) +
  labs(y = paste0("Palynomorph Darkness Index (%)"),
       x = paste0("Spore length (", "\U03BC", "m)"),
       title = "b)") +
  annotate("text", x = 14, y = -12.5, label = "PDI = -6.51 × Spore length + 216.2", hjust = 0, size = 5)

robust_model_6b <- robustbase::lmrob(PDI ~ Length_Mean, data = lycPDI)
summary(robust_model_6b)

######## MAIN TEXT FIGURE 7
wl_long <- wl %>% tidyr::gather("Study", "Mass_Loss", -Temp) %>% mutate(Mass_Loss = as.numeric(Mass_Loss))
wl_long <- wl_long %>% tidyr::separate(Study, c("Taxon", "Study"))

par(pty = "s")
par(mfrow = c(1, 1), mar = c(5.1, 5.1, 1.1, 6.1))
PDI_data <- pdi %>% group_by(Temp) %>% dplyr::summarise(PDI = mean(PDI))

wl_long_sig <- wl_long %>% na.omit()
x <- wl_long_sig$Temp
y <- wl_long_sig$Mass_Loss
logistic_fit <- nls(y ~ a / (1 + exp(-b * (x - c))),
                    start = list(a = max(y), b = 0.1, c = median(x)))

xPDI <- PDI_data$Temp
yPDI <- PDI_data$PDI

logistic_fit_PDI <- nls(yPDI ~ a / (1 + exp(-b * (xPDI - c))),
                        start = list(a = max(yPDI), b = 0.1, c = median(xPDI)))

plot(Mass_Loss ~ Temp, wl_long_sig, 
     type = "n",
     las = 1,
     cex.lab = 2,
     cex.axis = 1.5,
     xlab = "Pyrolysis temperature (\u00b0C)",
     ylab = "Mass loss (%)",
     axes = FALSE)
curve(coef(logistic_fit)["a"] / (1 + exp(-coef(logistic_fit)["b"] * (x - coef(logistic_fit)["c"]))),
      add = TRUE, lwd = 2, lty = 2)
points(Mass_Loss ~ Temp, wl_long_sig[wl_long_sig$Study == "This", ], pch = 16, cex = 2)
points(Mass_Loss ~ Temp, wl_long_sig[wl_long_sig$Study == "Rouxhet", ], cex = 2, lwd = 3)
legend("bottomright", 
       pch = c(16, 1, NA, NA, 4, NA), 
       col = c("black", "black", "black", NA, "darkred", "darkred"), 
       lty = c(NA, NA, 2, NA, NA, 2), 
       legend = c("Lycopodium, this study",
                  "Lycopodium sporopollenin (Rouxhet et al., 1979)",
                  "Logistic model of both mass loss datasets",
                  "",
                  "Palynomorph Darkness Index (PDI) data, this study",
                  "Logistic model of PDI data"), cex = 1.1, lwd = 2)
internal_axis_ticks(1:2, 1:2, cex.axis = 1.5)
par(new = T)
plot(PDI ~ Temp, PDI_data, axes = FALSE, xlab = "", ylab = "", type = "n")
curve(coef(logistic_fit_PDI)["a"] / (1 + exp(-coef(logistic_fit_PDI)["b"] * (x - coef(logistic_fit_PDI)["c"]))),
      add = TRUE, lwd = 2, lty = 2, col = "darkred")
points(PDI ~ Temp, PDI_data, col = "darkred", cex = 2, pch = 4, lwd = 3)
internal_axis_ticks(4, 4, cex.axis = 1.5)
mtext(4, 3.5, at = ((77-5.59)/2) + 5.59, cex = 2, text = "Palynomorph Darkness Index (%)", col = "darkred")

# Summary statistics of mass loss logistic model
summary(logistic_fit)

# Summary statistics of PDI logistic model
summary(logistic_fit_PDI)

#######################
######### SUPPLEMENTARY
#######################

######## SUPPLEMENTARY FIGURE S1
sp_cor_mat <- sp_cor_pearson <- sp_means %>% 
  select(-c(sd, upper, lower)) %>% 
  tidyr::spread("wavenumber", "mean") %>% 
  mutate(sample = strip_non_numbers(sample)) %>%
  tibble::column_to_rownames("sample") %>%
  prospectr::savitzkyGolay(p = 2, m = 2, w = 7) %>%
  t() %>% 
  cor(method = "pearson")

corrplot::corrplot(sp_cor_mat, method = "square", tl.col = "white", addCoef.col = "white", diag = F)
text(seq(11, 1, -1), seq(1, 11, 1), labels = rev(sort(c(300, 275, 350, 325, 375, 200, 150, 225, 175, 0, 250))), col = rev(overplot_palette), cex = 2.5)

rect(10.5, 1.5, 11.5, 0.5, border = "red", lwd = 5, lty = 2)
rect(0.5, 11.5, 6.5, 5.5, border = "red", lwd = 5, lty = 2)
rect(6.5, 5.5, 10.5, 1.5, border = "red", lwd = 5, lty = 2)

rect(0.5, 11.5, 3.5, 10.5, border = "red", lwd = 10)
rect(3.5, 6.5, 6.5, 5.5, border = "red", lwd = 10)
rect(6.5, 5.5, 9.5, 4.5, border = "red", lwd = 10)
rect(7.5, 2.5, 10.5, 1.5, border = "red", lwd = 10)
rect(9.5, 1.5, 11.5, 0.5, border = "red", lwd = 10)

mtext(3, -1.5, at = 1:6, text = "I", cex = 2)
mtext(3, -2, at = 1:6, text = "O", cex = 4)
mtext(3, -1.5, at = 7:10, text = "II", cex = 2)
mtext(3, -2, at = 7:10, text = "O", cex = 4)
mtext(3, -1.5, at = 11, text = "III", cex = 2)
mtext(3, -2, at = 11, text = "O", cex = 4)

######## SUPPLEMENTARY FIGURE S2
lutzkelyco <- read.csv("Lutzke Lycopodium sporopollenin spectrum.csv")
mylyco <- sp_means %>% filter(sample == "HFT0") %>% ungroup() %>% select(wavenumber, mean) %>% rename("absorbance" = "mean") %>% filter(wavenumber < 3600 & wavenumber > 948)

plot(absorbance ~ wavenumber, lutzkelyco, xlim = rev(range(wavenumber)), type = "n", axes = F, ylab = "Relative absorbance", xlab = expression(paste("Wavenumber (cm"^"-1"*")")), cex.lab = 2)
lines(absorbance ~ wavenumber, lutzkelyco, lwd = 2)
par(new = T)
plot(absorbance ~ wavenumber, mylyco, xlim = rev(range(wavenumber)), type = "n", axes = F, xlab = "", ylab = "")
lines(absorbance ~ wavenumber, mylyco, col = "darkred", lwd = 2)
internal_axis_ticks(1, 1, cex.axis = 1.5)
legend("topleft", lty = 1, , lwd = 2, col = c("black", "darkred"), legend = c("Lycopodium sporopollenin (Lutzke et al., 2019)", "Lycopodium exine (this study)"), cex = 1.5)

plot_deriv <- select(sp_means, -c(sd, upper, lower)) %>% 
  tidyr::spread("wavenumber", "mean") %>% 
  tibble::column_to_rownames("sample") %>% 
  prospectr::savitzkyGolay(w = 7, p = 2, m = 2) %>%
  data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  tidyr::gather("wavenumber", "mean", -sample) %>%
  mutate(wavenumber = strip_non_numbers(wavenumber))

######## SUPPLEMENTARY FIGURE S3
oleic <- read.csv("Lutzke oleic acid spectrum.CSV", header = F)
oleic <- to_absorbance(oleic)
names(oleic) <- c("wavenumber", "absorbance")
oleic <- oleic %>% filter(wavenumber < 1850 & wavenumber > 900)

g_tri <- read.csv("Lutzke glyceryl trioliate spectrum.CSV", header = F)
g_tri <- to_absorbance(g_tri)
names(g_tri) <- c("wavenumber", "absorbance")
g_tri <- g_tri %>% filter(wavenumber < 1850 & wavenumber > 900)

range_x <- c(1800, 950) 

mar(n = 1.1)
mfrow(1, 1)
par(pty = "m")
plot(absorbance ~ wavenumber, oleic, 
     xlim = range_x,
     ylim = inRange(oleic$absorbance, c(-1, 1.1)), type = "n", axes = F, 
     col = scico::scico(2, palette = "roma")[2], lwd = 2, 
     xlab = expression(paste("Wavenumber (cm"^"-1"*")")), cex.lab = 2, ylab = "")
rect(1738, 10, 1746, -10, col = "#00000055", border = NA)
rect(1705, 10, 1713, -10, col = "#00000055", border = NA)
rect(1157, 10, 1163, -10, col = "#00000055", border = NA)
lines(absorbance ~ wavenumber, oleic, 
      col = scico::scico(2, palette = "roma")[2], lwd = 2)

text(absorbance ~ wavenumber, filter(oleic, wavenumber == 1707.658), label = paste0(round(wavenumber, 0), "\nv(C=O)"), pos = 3, col = scico::scico(2, palette = "roma")[2])
par(new = T)
plot(absorbance ~ wavenumber, g_tri, 
     xlim = range_x,
     type = "l", axes = F, col = scico::scico(2, palette = "roma")[1], lwd = 2, 
     ylim = inRange(g_tri$absorbance, c(-1, 1.1)), 
     xlab = "", ylab = "")
text(absorbance ~ wavenumber, filter(g_tri, wavenumber == 1159.9730), label = paste0(round(wavenumber, 0), "\nv(CO)"), pos = 3, col = scico::scico(2, palette = "roma")[1])
text(absorbance ~ wavenumber, filter(g_tri, wavenumber == 1743.334), label = paste0(round(wavenumber, 0), "\nv(C=O)"), pos = 3, col = scico::scico(2, palette = "roma")[1])

internal_axis_ticks(1, 1, cex.axis = 1.5)
legend("topright", lty = c(1, 1, NA, NA), lwd = 2, col = rev(scico::scico(2, palette = "roma")), legend = c("Oleic acid", "Glyceryl trioliate", "", ""), cex = 1.5)
text(1000, 0.22, labels = "Monomer spectra from\nLutzke et al., 2020a")

overplot_palette <- c(scico::scico(length(unique(sp_means$sample)) + 6, palette = "roma")[rev(c(1:4, 12:17))], "black")
par(new = T)
plot(mean ~ wavenumber, sp_means[sp_means$sample == unique(sp_means$sample)[1], ], 
     type = "n",
     xlim = range_x,
     ylim = inRange(sp_means$mean, c(0, 2.05)),
     las = 1,
     axes = F,
     xlab = "",
     ylab = "",
     cex.lab = 2,
     cex.axis = 1.5)

for(i in 1:(length(unique(sp_means$sample)) - 1)){
  lines(mean ~ wavenumber, sp_means[sp_means$sample == unique(sp_means$sample)[i], ],
        col = overplot_palette[i],
        lwd = 2)
}
mtext(2, 1.5, at = inRange(sp_means$mean, 1.025), text = "Relative absorbance", cex = 2)
legend("bottomright", lty = 1, lwd = 2, col = overplot_palette[-length(overplot_palette)], legend = publication_labels[-length(overplot_palette)], cex = 0.8)


######## SUPPLEMENTARY FIGURE S4 PANE 1
plot_all_means(plot_deriv, 3025, 2825, plot_legend = "bottomleft")
plot_wn_label(3025, 2825, 3006, label = "D", temp = "375")
plot_wn_label(3025, 2825, 2975, label = "O", temp = "250", arrow_size = 1)
plot_wn_label(3025, 2825, 2908, label = "D", temp = "275")
plot_wn_label(3025, 2825, 2894, label = "A", temp = "275")
plot_wn_label(3025, 2825, 2843, label = "O", temp = "250", arrow_size = 3)
plot_wn_shift(wn_hi = 3025, wn_lo = 2825, wn_from = 2874, wn_to = 2870, temp_from = 0, temp_to = 375, bar_space = 36, low_text = 375)
other_peaks1 <- c(2957, 2925, 2854)
purrr::walk(other_peaks1, .f = function(x){plot_peak_label(wn_hi = 3025, wn_lo = 2825, wn = x)})

######## SUPPLEMENTARY FIGURE S4 PANE 2
plot_all_means(plot_deriv, 1775, 1575, plot_legend = "none")
plot_wn_label(1775, 1575, 1770, label = "A", temp = "375")
plot_wn_label(1775, 1575, 1727, label = "D", temp = "275")
plot_wn_label(1775, 1575, 1676, label = "O", temp = "100")
plot_wn_shift(wn_hi = 1775, wn_lo = 1575, wn_from = 1745, wn_to = 1739, temp_from = 250, temp_to = 375, low_text = "275")
plot_wn_shift(wn_hi = 1775, wn_lo = 1575, wn_from = 1604, wn_to = 1593, temp_from = 0, temp_to = 375, bar_space = 12, low_text = "275")
other_peaks2 <- c(1710, 1654, 1635)
purrr::walk(other_peaks2, .f = function(x){plot_peak_label(wn_hi = 1775, wn_lo = 1575, wn = x)})

######## SUPPLEMENTARY FIGURE S4 PANE 3
plot_all_means(plot_deriv, 1600, 1375, plot_legend = "none")
plot_wn_label(1600, 1375, 1558, label = "A", temp = "275")
plot_wn_label(1600, 1375, 1541, label = "A", temp = "250", arrow_size = 2)
plot_wn_label(1600, 1375, 1531, label = "D", temp = "300")
plot_wn_label(1600, 1375, 1508, label = "A", temp = "375")
plot_wn_label(1600, 1375, 1402, label = "A", temp = "300")
plot_wn_shift(wn_hi = 1600, wn_lo = 1375, wn_from = 1517, wn_to = 1514, temp_from = 175, temp_to = 350, bar_space = 12, low_text = "275")
other_peaks3 <- c(1490, 1466, 1457, 1436, 1415)
purrr::walk(other_peaks3, .f = function(x){plot_peak_label(wn_hi = 1600, wn_lo = 1375, wn = x)})

######## SUPPLEMENTARY FIGURE S4 PANE 4
plot_all_means(plot_deriv, 1400, 1200, plot_legend = "none")
plot_wn_shift(wn_hi = 1400, wn_lo = 1200, wn_from = 1350, wn_to = 1342, temp_from = 175, temp_to = 375, bar_space = 12, l2r = F, low_text = "275")
plot_wn_shift(wn_hi = 1400, wn_lo = 1200, wn_from = 1308, wn_to = 1300, temp_from = 175, temp_to = 375, bar_space = 12, l2r = T, low_text = "275")
plot_wn_shift(wn_hi = 1400, wn_lo = 1200, wn_from = 1323, wn_to = 1320, temp_from = 175, temp_to = 375, bar_space = 12, l2r = T, low_text = "275")
plot_wn_shift(wn_hi = 1400, wn_lo = 1200, wn_from = 1285, wn_to = 1275, temp_from = 175, temp_to = 375, bar_space = 15, l2r = T, low_text = "225")
plot_wn_shift(wn_hi = 1400, wn_lo = 1200, wn_from = 1240, wn_to = 1234, temp_from = 0, temp_to = 350, bar_space = 36, l2r = F, low_text = "275")
plot_wn_label(1400, 1200, 1260, label = "D", temp = "350")
plot_wn_label(1400, 1200, 1248, label = "A", temp = "375")
other_peaks4 <- c(1379)
purrr::walk(other_peaks4, .f = function(x){plot_peak_label(wn_hi = 1400, wn_lo = 1200, wn = x)})

######## SUPPLEMENTARY FIGURE S4 PANE 5
plot_all_means(plot_deriv, 1200, 950, plot_legend = "none")
plot_wn_label(1200, 950, 1154, label = "A", temp = "375", arrow_size = 2)
plot_wn_label(1200, 950, 968, label = "A", temp = "250")
plot_wn_shift(wn_hi = 1200, wn_lo = 950, wn_from = 1119, wn_to = 1112, temp_from = 0, temp_to = 375, bar_space = 20, l2r = F, low_text = "375")
plot_wn_shift(wn_hi = 1200, wn_lo = 950, wn_from = 1096, wn_to = 1092, temp_from = 0, temp_to = 350, bar_space = 12, l2r = T, low_text = "275")
plot_wn_shift(wn_hi = 1200, wn_lo = 950, wn_from = 1060, wn_to = 1051, temp_from = 0, temp_to = 375, bar_space = 8, l2r = F, low_text = "275")
plot_wn_shift(wn_hi = 1200, wn_lo = 950, wn_from = 996, wn_to = 987, temp_from = 0, temp_to = 350, bar_space = 12, l2r = T, low_text = "275")
other_peaks5 <- c(1167, 1141, 1034)
purrr::walk(other_peaks5, .f = function(x){plot_peak_label(wn_hi = 1200, wn_lo = 950, wn = x)})

######## SUPPLEMENTARY FIGURE S5
HCO <- rouxhet %>% select(Temp, WeightLoss, Perc_H, Perc_C, Perc_O) %>% mutate(Perc_Other = 100 - Perc_H - Perc_C - Perc_O)
HCO %>% 
  tidyr::gather("Element", "Perc", -Temp, -WeightLoss) %>%
  mutate(Element = str_replace_all(Element, "Perc_", "")) %>%
  mutate(Element = factor(Element, levels = c("Other", "H", "O", "C"))) %>%
  na.omit() %>% 
  ggplot(aes(Temp, Perc, fill = Element)) + 
  geom_bar(stat = "identity", position = "stack", width = 20, colour = "darkgrey") +
  scale_fill_manual(values = scico::scico(4, palette = "roma")) +
  theme_bw() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = "top") +
  labs(x = "Pyrolysis temperature (\u00b0C)",
       y = "Proportion of total mass (%)") +
  scale_y_continuous(breaks = seq(0, 100, 10))

######## SUPPLEMENTARY FIGURE S6

######## SUPPLEMENTARY FIGURE S7
TernaryPlot(alab = "Redder \u2192", blab = "\u2190 Greener", clab = "Bluer \u2192",
            lab.col = c("red", "darkgreen", "blue"),
            main = "",
            point = "right", lab.cex = 0.8, grid.minor.lines = 0,
            grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), grid.col = "white", 
            axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
            axis.rotate = FALSE,
            padding = 0.06,
            isometric = T)

cols <- TernaryPointValues(rgb)
cols2 <- cols
cols2[3, ] <- paste0(cols2[3, ], "55")
ColourTernary(cols2, spectrum = NULL)

pdi_ternary <- vector("list", nrow(centroids))
names(pdi_ternary) <- centroids$Temp
rgb_means <- pdi_short[, c("Temp", "R", "G", "B")] %>% group_by(Temp) %>% summarise_all("mean")
for(i in 1:nrow(centroids)){
  pdi_ternary[[i]] <- as.numeric(rgb_means[i, -1])
}
data_points <- pdi_ternary[1:11]

AddToTernary(graphics::points, data_points, pch = 21, cex = 2, 
             bg = vapply(data_points, 
                         function (x) rgb(x[1], x[2], x[3], 128,
                                          maxColorValue = 255),
                         character(1))
)

