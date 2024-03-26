#!/usr/bin/env Rscript

## DOC ##
' Plot a simple histogram from a vector of values taken from stdin.

usage: quick_histogram.R [-v -m <msg> --breaks=<breaks>] -o <outprefix> -t <title> <file_arg> 

options:
 -v        verbose
 -m <msg>  Message
 -b --breaks <breaks> Number of breaks for the histogram' -> doc

#### LIBS ####
library(docopt)
opts <- docopt(doc)

#### ARGS ####
if (opts$v) print(str(opts)) 
if (!is.null(opts$message)) cat("MESSAGE: ", opts$m)

#### FUNS ####
## File Read ##
OpenRead <- function(arg) {
  if (arg %in% c("-", "/dev/stdin")) {
    file("stdin", open = "r")
  } else if (grepl("^/dev/fd/", arg)) {
    fifo(arg, open = "r")
  } else {
    file(arg, open = "r")
  }
}

#### MAIN ####
# read the data
dat.con <- OpenRead(opts$file_arg)
dat <- read.table(dat.con, sep = " ", header = FALSE)

# do something with dat and opts$param
png(paste0(opts$outprefix, ".png"))
# pdf(paste0(opts$outprefix, ".pdf"))
if(!is.null(opts$breaks)) {
  hist(dat$V1,
       main = opts$title,
       breaks = as.integer(opts$breaks),
       # xlim = c(0, 10),
       xlab = NULL)
} else {
  hist(dat$V1,
       main = opts$title,
       # xlim = c(0, 10),
       xlab = NULL)
}
dev.off()

#### END ####