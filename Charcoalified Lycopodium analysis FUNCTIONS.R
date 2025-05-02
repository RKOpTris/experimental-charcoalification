##################################
##### MAIN TEXT
##################################

remove.co2.signal <- function(x, co2.hi = 2410, co2.lo = 2200){
  if(is.data.frame(x)){
    co2.hi.index <- nearest(names(x), co2.hi)$index
    co2.lo.index <- nearest(names(x), co2.lo)$index
    y <- x
  } else {
    co2.hi.index <- which.min(abs(co2.hi - x@wavelength))
    co2.lo.index <- which.min(abs(co2.lo - x@wavelength))
    y <- x$spc
  }
  for(r in 1:nrow(y)){
    y[r, co2.lo.index:co2.hi.index] <- seq(from = y[r, co2.lo.index], to = y[r, co2.hi.index], length.out = length(y[r, co2.lo.index:co2.hi.index]))
  }
  y
}

nearest <- function(vector, val){
  stopifnot(is.numeric(vector) | is.character(vector))
  if(is.numeric(vector)){
    ind <- which.min(abs(val - vector))
    data.frame(
      index = ind,
      value = vector[ind],
      stringsAsFactors = F
    )
  } else {
    ind <- which.min(abs(val - strip_non_numbers(vector)))
    data.frame(
      index = ind,
      name = vector[ind],
      stringsAsFactors = F
    )
  }
}

strip_non_numbers <- function(vector){
  gsub("[^0-9.-]", "", vector) %>% as.numeric()
}

find_peaks <- function(x, wns = wavenumbers, type = "both", nudge = 0, nups = 5, ndowns = 1, ...){
  if(type %in% c("both", "peak")){
    spec_diff_peaks <- pracma::findpeaks(x, nups = nups, ...) %>% data.frame()
    names(spec_diff_peaks) <- c("abs", "peak_ind", "peak_start", "peak_end")
    spec_diff_peaks$peak_wavenumber <- wns[spec_diff_peaks$peak_ind]
    spec_diff_peaks$peak_type <- "peak"
    spec_diff_peaks$nudge <- spec_diff_peaks$abs + nudge
    peaks <- spec_diff_peaks
  }
  if(type %in% c("both", "trough")){
    spec_diff_peaks <- pracma::findpeaks(-1 * x, nups = ndowns, ...) %>% data.frame()
    names(spec_diff_peaks) <- c("abs", "peak_ind", "peak_start", "peak_end")
    spec_diff_peaks$peak_wavenumber <- wns[spec_diff_peaks$peak_ind]
    spec_diff_peaks$peak_type <- "trough"
    spec_diff_peaks$nudge <- (spec_diff_peaks$abs + nudge) * -1
    spec_diff_peaks$abs <- spec_diff_peaks$abs * -1
  }
  if(type == "both"){
    bind_rows(peaks, spec_diff_peaks) %>% arrange(peak_wavenumber)
  } else {
    spec_diff_peaks
  }
}

inRange <- function(x, value = 0.5){
  stopifnot(is.numeric(x) | lubridate::is.Date(x))
  if(lubridate::is.Date(x)){
    was.Date <- T
    x <- as.numeric(x)
  } else {
    was.Date <- F
  }
  out <- sapply(value, function(v){min(x) + ((max(x) - min(x)) * v)})
  if(was.Date){
    as.Date(out, origin = "1970-01-01")
  } else {
    out
  }
}

mar <- function(s = 5.1, w = 4.1, n = 4.1, e = 2.1){
  par(mar = c(s, w, n, e))
}

mfrow <- function(r = 1, c = 1){
  par(mfrow = c(r, c))
}

plot_pars <- function(r = 1, c = 1, s = 5.1, w = 4.1, n = 1.1, e = 1.1, pty = "m", reset = F){
  if(reset){
    mfrow()
    mar()
    par()
  } else {
    mfrow(r = r, c = c)
    mar(s = s, w = w, n = n, e = e)
    par(pty = pty)
  }
}

internal_axis_ticks <- function(draw.axes = 1:4, number.axes = 1:2, ...){
  for(n in draw.axes){
    axis(n, tck = 0.01, labels = NA, ...)
  }
  for(n in number.axes){
    axis(n, lwd = 0, line = -0.5, las = 1, ...)
  }
  box(lwd = 1.5)
}

draw_bands <- function(){                               
  rect(absorbance_bands$xmin,
       absorbance_bands$ymin,
       absorbance_bands$xmax,
       absorbance_bands$ymax,
       border = NA,
       col = "#00000022")
}

labels_for_bands <- function(y){
  text(inRange(c(absorbance_bands$xmin[1], absorbance_bands$xmax[1]), 0.5), y, 
       labels = expression(paste(nu, "(CH"[2] * ")")), srt = 90, adj = 1)
  text(inRange(c(absorbance_bands$xmin[2], absorbance_bands$xmax[2]), 0.5), y, 
       labels = expression(paste(nu, "(C=C)")), srt = 90, adj = 1)
  text(inRange(c(absorbance_bands$xmin[3], absorbance_bands$xmax[3]), 0.5), y, 
       labels = expression(paste(nu, "(C=O)")), srt = 90, adj = 1)
  text(inRange(c(absorbance_bands$xmin[4], absorbance_bands$xmax[4]), 0.5), y, 
       labels = expression(paste(delta["scissors"]*"(CH"[2]*")")), srt = 90, adj = 1)
  text(inRange(c(absorbance_bands$xmin[5], absorbance_bands$xmax[5]), 0.5), y, 
       labels = expression(paste(nu, "(CO)")), srt = 90, adj = 1)
  text(inRange(c(absorbance_bands$xmin[6], absorbance_bands$xmax[6]), 0.5), y, 
       labels = expression(paste(nu, "(CO)")), srt = 90, adj = 1)
  text(inRange(c(absorbance_bands$xmin[7], absorbance_bands$xmax[7]), 0.5), y, 
       labels = expression(paste("Methine/", gamma["wagging"]*"(CH"[2]*")")), srt = 90, adj = 1)
  text(inRange(c(absorbance_bands$xmin[8], absorbance_bands$xmax[8]), 0.5), y, 
       labels = expression(paste(nu, "(OH)")), srt = 90, adj = 1)
}

pn_max <- function(x){
  if(mean(x) > 0){
    max(x)
  } else {
    min(x)
  }
}

##### summarySE function from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##################################
##### SUPPLEMENTARY
##################################

plot_wn_label <- function(wn_hi = 4000, wn_lo = 950, tgt_wn, col = "black", label = "A", temp = "100", location = "below", arrow_size = 1){
  range_x <- c(wn_hi, wn_lo)
  plot_deriv_subset <- plot_deriv %>% filter(wavenumber <= range_x[1] & wavenumber >= range_x[2])
  plot_unit <- diff(inRange(plot_deriv_subset$mean, c(-0.1, 1.0))) / 100
  peak_pos <- min(plot_deriv_subset$mean[plot_deriv_subset$wavenumber == nearest(plot_deriv_subset$wavenumber, tgt_wn)$value])
  arrow_hi <- peak_pos - (3 * plot_unit)
  arrow_lo <- arrow_hi - (arrow_size * 10 * plot_unit)
  shape::Arrows(tgt_wn, arrow_lo, tgt_wn, arrow_hi, arr.type = "triangle", lwd = 3)
  text(tgt_wn, arrow_lo - (2 * plot_unit), label = tgt_wn, cex = 0.8)
  text(tgt_wn, arrow_lo - (4 * plot_unit), label = label, cex = 1.5)
  text(tgt_wn, arrow_lo - (6 * plot_unit), label = paste(temp, "\u00b0C"), cex = 0.8)
}

plot_wn_shift <- function(wn_hi = 4000, wn_lo = 950, wn_from, wn_to, temp_from, temp_to, bar_space = 8, l2r = T, low_text = ""){
  range_x <- c(wn_hi, wn_lo)
  plot_deriv_subset <- plot_deriv %>% filter(wavenumber <= range_x[1] & wavenumber >= range_x[2]) %>%
    mutate(sample = strip_non_numbers(sample))
  range_y <- inRange(plot_deriv_subset$mean, c(-0.1, 1.0))
  plot_unit <- diff(range_y) / 100
  wn_from <- nearest(plot_deriv_subset$wavenumber, wn_from)$value
  wn_to <- nearest(plot_deriv_subset$wavenumber, wn_to)$value
  lowest_peak <- min(plot_deriv_subset[plot_deriv_subset$wavenumber < wn_from & plot_deriv_subset$wavenumber > wn_to, ]$mean)
  bar_y <- lowest_peak - (bar_space * plot_unit)
  bar_l <- plot_deriv_subset[plot_deriv_subset$wavenumber == wn_from & plot_deriv_subset$sample == temp_from, ]$mean
  bar_r <- plot_deriv_subset[plot_deriv_subset$wavenumber == wn_to & plot_deriv_subset$sample == temp_to, ]$mean
  if(l2r){
    arrows(wn_from, bar_y, wn_to, bar_y, length = 0.15, lwd = 3, lty = 1)
    bar_l <- plot_deriv_subset[plot_deriv_subset$wavenumber == wn_from & plot_deriv_subset$sample == temp_from, ]$mean
    bar_r <- plot_deriv_subset[plot_deriv_subset$wavenumber == wn_to & plot_deriv_subset$sample == temp_to, ]$mean
  } else {
    arrows(wn_to, bar_y, wn_from, bar_y, length = 0.15, lwd = 3, lty = 1)
    bar_r <- plot_deriv_subset[plot_deriv_subset$wavenumber == wn_to & plot_deriv_subset$sample == temp_from, ]$mean
    bar_l <- plot_deriv_subset[plot_deriv_subset$wavenumber == wn_from & plot_deriv_subset$sample == temp_to, ]$mean
  }
  arrows(wn_from, bar_y, wn_from, bar_l, length = 0, lwd = 2, lty = 3)
  arrows(wn_to, bar_y, wn_to, bar_r, length = 0, lwd = 2, lty = 3)
  text(wn_from, bar_y, pos = 2, labels = round(wn_from, 0), cex = 0.7)
  text(wn_to, bar_y, pos = 4, labels = round(wn_to, 0), cex = 0.7)
  text(inRange(rev(c(wn_from, wn_to)), 0.5), bar_y - plot_unit, pos = 1, labels = paste0(low_text, " \u00b0C"), cex = 0.8)
}

plot_peak_label <- function(wn_hi = 4000, wn_lo = 950, wn, ...){
  range_x <- c(wn_hi, wn_lo)
  plot_deriv_subset <- plot_deriv %>% filter(wavenumber <= range_x[1] & wavenumber >= range_x[2]) %>%
    mutate(sample = strip_non_numbers(sample))
  range_y <- inRange(plot_deriv_subset$mean, c(-0.1, 1.0))
  plot_unit <- diff(range_y) / 100
  lowest_peak <- min(plot_deriv_subset[plot_deriv_subset$wavenumber == nearest(plot_deriv_subset$wavenumber, wn)$value, ]$mean)
  text(wn, lowest_peak - plot_unit, pos = 1, labels = round(wn, 0), cex = 0.8, ...)
}

plot_all_means <- function(tidy_spectra, wn_hi = 4000, wn_lo = 950, plot_legend = "none", overlay = F, invert = F){
  range_x <- c(wn_hi, wn_lo) 
  mar(w = 4.1, n = 1.1, e = 1.1)
  par(pty = "m")
  overplot_palette <- c(scico::scico(length(unique(tidy_spectra$sample)) + 6, palette = "roma")[rev(c(1:4, 12:17))], "black")
  tidy_spectra_subset <- tidy_spectra %>% filter(wavenumber <= range_x[1] & wavenumber >= range_x[2])
  if(invert){
    tidy_spectra_subset$mean <- tidy_spectra_subset$mean * -1
  }
  if(overlay){
    par(new = T)
  }
  plot(mean ~ wavenumber, tidy_spectra_subset[tidy_spectra_subset$sample == unique(tidy_spectra_subset$sample)[1], ], 
       type = "n",
       xlim = range_x,
       ylim = inRange(range(tidy_spectra_subset$mean), c(-0.1, 1.0)),
       las = 1,
       axes = F,
       xlab = expression(paste("Wavenumber (cm"^"-1"*")")),
       ylab = "",
       cex.lab = 2)
  abline(h = 0, col = "lightgrey")
  
  for(i in 1:(length(unique(tidy_spectra_subset$sample)))){
    lines(mean ~ wavenumber, tidy_spectra_subset[tidy_spectra_subset$sample == unique(tidy_spectra_subset$sample)[i], ],
          col = overplot_palette[i],
          lwd = 2)
  }
  internal_axis_ticks(1, 1)
  mtext(2, 1.5, at = inRange(tidy_spectra_subset$mean, 0.5), text = "Relative absorbance", cex = 2)
  
  if(plot_legend != "none"){
    legend(plot_legend, lty = 1, lwd = 3, col = c(overplot_palette[1:6], NA, overplot_palette[7:10], NA, overplot_palette[11]), legend = c(publication_labels[1:6], NA, publication_labels[7:10], NA, publication_labels[11]), cex = 1.2)
  }
}

to_absorbance <- function(df, ...){
  # assuming wavenumbers are 1st column!
  df[-1] <- lapply(df[-1], trans_to_abs, ...)
  df
}

trans_to_abs <- function(vector, infinite_as = NA){
  out <- -log10(vector / 100)
  out[is.infinite(out)] <- infinite_as
  out
}

