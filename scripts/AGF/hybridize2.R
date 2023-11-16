hybridize2<-function (x1, x2, n, pop = "hybrid", res.type = c("genind", "df", 
                                                  "STRUCTURE"), file = NULL, quiet = FALSE, sep = "/", hyb.label = "h") 
{
  if (!is.genind(x1)) 
    stop("x1 is not a valid genind object")
  if (!is.genind(x2)) 
    stop("x2 is not a valid genind object")
  if (!all(ploidy(x1) == ploidy(x1)[1])) 
    stop("varying ploidy (in x1) is not supported for this function")
  if (!all(ploidy(x2) == ploidy(x2)[1])) 
    stop("varying ploidy (in x2) is not supported for this function")
  if (ploidy(x1)[1]%%2 != 0) 
    stop("not implemented for odd levels of ploidy")
  if (ploidy(x1)[1] != ploidy(x2)[1]) 
    stop("x1 and x2 have different ploidy")
  checkType(x1)
  checkType(x2)
  n <- as.integer(n)
  ploidy <- ploidy(x1)[1]
  res.type <- match.arg(res.type)
  popNames(x1) <- "pop1"
  popNames(x2) <- "pop2"
  x1x2 <- repool(x1, x2)
  x1x2<-missingno(x1x2, type = "loci", cutoff = 0.0001, quiet = FALSE, freq = FALSE)
  x1 <- x1x2[pop = 1]
  x2 <- x1x2[pop = 2]
  n1 <- nInd(x1)
  n2 <- nInd(x2)
  k <- nLoc(x1)
  y1 <- genind2genpop(x1, pop = factor(rep(1, n1)), quiet = TRUE)
  freq1 <- tab(y1, freq = TRUE)
  freq1 <- split(freq1, y1@loc.fac)
  freq1 <- freq1[locNames(x1)]
  y2 <- genind2genpop(x2, pop = factor(rep(1, n2)), quiet = TRUE)
  freq2 <-tab(y2, freq = TRUE)
  freq2 <- split(freq2, y2@loc.fac)
  freq2 <- freq2[locNames(x2)]
  kX1 <- lapply(freq2, function(v)  try(t(rmultinom(n, ploidy/2, 
                                               v))))
  
  names(kX1) <- locNames(x1)
  vec.paste1 <- NULL
  Vec.all1 <- NULL
  for (i in 1:k) {
    colnames(kX1[[i]]) <- alleles(x1)[[i]]
    vec.paste1 <- c(vec.paste1, alleles(x1)[[i]])
    Vec.all1 <- c(Vec.all1, length(alleles(x1)[[i]]))
  }
  kX2 <- lapply(freq1, function(v) t(rmultinom(n, ploidy/2, 
                                               v)))
  names(kX2) <- locNames(x2)
  vec.paste2 <- NULL
  Vec.all2 <- NULL
  for (i in 1:k) {
    colnames(kX2[[i]]) <- alleles(x2)[[i]]
    vec.paste2 <- c(vec.paste2, alleles(x2)[[i]])
    Vec.all2 <- c(Vec.all2, length(alleles(x2)[[i]]))
  }
  tab1 <- as.matrix(cbind.data.frame(kX1))
  colnames(tab1) <- paste(rep(locNames(x1), Vec.all1), ".", 
                          vec.paste1, sep = "")
  tab2 <- as.matrix(cbind.data.frame(kX2))
  colnames(tab2) <- paste(rep(locNames(x2), Vec.all2), ".", 
                          vec.paste2, sep = "")
  zyg.rownames <- .genlab(hyb.label, n)
  zyg.colnames <- sort(unique(c(colnames(tab1), colnames(tab2))))
  zyg <- matrix(0, nrow = n, ncol = length(zyg.colnames), dimnames = list(zyg.rownames, 
                                                                          zyg.colnames))
  zyg[, colnames(tab1)] <- zyg[, colnames(tab1)] + tab1
  zyg[, colnames(tab2)] <- zyg[, colnames(tab2)] + tab2
  zyg <- zyg
  zyg <- genind(zyg, type = "codom", ploidy = ploidy)
  if (res.type == "STRUCTURE") {
    temp <- genind2df(repool(x1, x2, zyg), usepop = FALSE, 
                      sep = " ")
    res <- unlist(apply(temp, 1, strsplit, " "))
    res <- as.data.frame(matrix(res, nrow = nrow(temp), byrow = TRUE))
    colnames(res) <- rep(colnames(temp), each = ploidy)
    res[is.na(res)] <- "-9"
    pop <- rep(1:3, c(nrow(x1@tab), nrow(x2@tab), n))
    res <- cbind.data.frame(pop, res, stringsAsFactors = FALSE)
    names(res)[1] <- ""
    if (is.null(file)) {
      file <- gsub("[[:space:]]|:", "-", date())
      file <- paste("hybrid", file, sep = "_")
      file <- paste(file, "str", sep = ".")
    }
    write.table(res, file = file, row.names = TRUE, col.names = TRUE, 
                quote = FALSE)
    if (!quiet) 
      cat("\nWrote results to file", file, "\n")
    return(invisible())
  }
  if (res.type == "df") {
    res <- genind2df(zyg, sep = sep)
    return(res)
  }
  if (res.type == "genind") {
    pop <- factor(rep(pop, n))
    res <- zyg
    pop(res) <- pop
    res@call <- match.call()
    return(res)
  }
}


