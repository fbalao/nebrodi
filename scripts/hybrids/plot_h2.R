plot_h2<-function (data, test.subject = "INDLABEL", POPID.name = "POPID", 
          mean.h.by = NULL, sort.by = "h_posterior_mode", col.group = NULL, 
          group.sep = NULL, fill.source = FALSE, basic.lines = TRUE, 
          source.col = NULL, source.limits = NULL, custom.abline = NULL, 
          labelx="Individual",labely="Hybrid index", ...) 
{
  setkey(data, NULL)
  if (is.null(col.group) == FALSE) {
    rf = colorRampPalette(rev(RColorBrewer::brewer.pal(8, 
                                                       "Dark2")))
    col.Dark2 = rf(nrow(unique(data[, col.group, with = F])))
    col.Dark2 = data.table(col.Dark2)
    setkeyv(data, col.group)
    pop = data.table(unique(data[, col.group, with = F]))
    pop.col = cbind(pop, col.Dark2)
    pop.col2<-pop.col
    pop.col2[pop.col$col.Dark=="#1B9E77", 2] <-"#D7191C"
    pop.col2[pop.col$col.Dark=="#9A58A5", 2] <-"#ABDDA4"
    pop.col<-pop.col2
    setkeyv(pop.col, col.group)
    data = data[pop.col]
    setkey(data, NULL)
    setkey(data, col.Dark2)
    col.summary = data[, head(.SD, 1), by = col.Dark2]
  }
  setkey(data, NULL)
  if (is.na(chmatch("Source", sort.by)) == FALSE) {
    setkey(data, Source)
    data["S0", `:=`(Source2, "A")]
    data["TEST", `:=`(Source2, "B")]
    data["S1", `:=`(Source2, "C")]
    sort.by[chmatch("Source", sort.by)] = "Source2"
  }
 
  
  setkey(data, NULL)
  setkeyv(data, test.subject)
  if (is.null(mean.h.by)) {
    setkeyv(data, sort.by)
    data[, `:=`(rn, row(data)[, 1])]
  }
  else {
    setkeyv(data, mean.h.by)
    data[, `:=`(mean_h, mean(h_posterior_mode)), by = mean.h.by]
    setkeyv(data, sort.by)
    data[, `:=`(rn, row(data)[, 1])]
  }
  setkey(data, NULL)
  plot(data[, h_posterior_mode] ~ data[, rn], type = "n", 
       xlab = labelx, ylab = labely, ...)
  if (is.null(group.sep) == FALSE) {
    setkeyv(data, group.sep)
    abline(v = 0.5, col = "grey")
    abline(v = data[, (max(rn) + 0.5), by = group.sep]$V1, 
           col = "grey")
  }
  setkey(data, NULL)
  setnames(data, POPID.name, "POPID", skip_absent = TRUE)
  if (fill.source == TRUE) {
    setkey(data, Source)
    for (i in 1:(data[Source == "S0", length(unique(POPID))])) {
      rect(data[Source == "S0", (min(rn)), by = POPID]$V1[i] - 
             0.5, par("usr")[3], data[Source == "S0", (max(rn)), 
                                      by = POPID]$V1[i] + 0.5, par("usr")[4], col = "grey95")
    }
    for (i in 1:(data[Source == "S1", length(unique(POPID))])) {
      rect(data[Source == "S1", (min(rn)), by = POPID]$V1[i] - 
             0.5, par("usr")[3], data[Source == "S1", (max(rn)), 
                                      by = POPID]$V1[i] + 0.5, par("usr")[4], col = "grey95")
    }
  }
  setkey(data, NULL)
  if (basic.lines == TRUE) {
    abline(h = 0.5, col = "grey", lty = 2, lwd = 2)
  }
  if (is.null(source.col) == FALSE) {
    setkey(data, Source)
    if (is.null(col.group)) {
      data[, `:=`(col.Dark2, "black")]
    }
    data["S0", `:=`(col.Dark2, source.col[1])]
    data["S1", `:=`(col.Dark2, source.col[2])]
    setkey(data, NULL)
    setkey(data, col.Dark2)
    col.summary = data[, head(.SD, 1), by = col.Dark2]
    setkey(col.summary, Source)
    col.summary["S0", `:=`(POPID, "S0")]
    col.summary["S1", `:=`(POPID, "S1")]
  }
  setkey(data, NULL)
  if (is.null(col.group) & is.null(source.col)) {
    points(data[, h_posterior_mode] ~ data[, rn], ...)
    arrows(data[, rn], data[, h_cred_int_lower], data[, 
                                                      rn], data[, h_cred_int_upper], angle = 90, code = 3, 
           length = 0)
  }
  else {
    points(data[, h_posterior_mode] ~ data[, rn], col="grey35", bg = data[,col.Dark2], ...)
    arrows(data[, rn], data[, h_cred_int_lower], data[, 
                                                      rn], data[, h_cred_int_upper], angle = 90, code = 3, 
           length = 0, col = data[, col.Dark2])
  }
  if (is.null(source.limits) == FALSE) {
    abline(h = data[Source == "S0", max(h_cred_int_upper)], 
           col = source.limits[1], lty = 2, lwd = 2)
    abline(h = data[Source == "S1", min(h_cred_int_lower)], 
           col = source.limits[2], lty = 2, lwd = 2)
  }
  if (is.null(custom.abline) == FALSE) {
    custom.abline
  }
  setnames(col.summary, "POPID", POPID.name)
  if (is.null(mean.h.by)) {
    col.summary = col.summary[, c(test.subject, "Source", 
                                  "col.Dark2", "rn"), with = F]
  }
  else {
    col.summary = col.summary[, c(test.subject, mean.h.by, 
                                  "Source", "mean_h", "col.Dark2", "rn"), with = F]
  }
  return(col.summary)
}
