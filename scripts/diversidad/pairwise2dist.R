closeAllConnections()

x<-rela[,c(2,3,9)]
# get the unique names
x.names <- sort(unique(c(x[[1]], x[[2]])))
 # create a matrix of the right size and put names on it
x.dist <- matrix(0, length(x.names), length(x.names))
dimnames(x.dist) <- list(x.names, x.names)
# create indices by converting names to numbers and create the normal and reversed
# to fill in all the matrix
 x.ind <- rbind(cbind(match(x[[1]], x.names), match(x[[2]], x.names)),
                                     cbind(match(x[[2]], x.names), match(x[[1]], x.names)))
x.dist[x.ind] <- rep(x[[3]], 2)
 x.dist
 