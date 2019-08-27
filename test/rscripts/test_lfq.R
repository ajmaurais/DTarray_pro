
library(ggplot2)

dat <- read.csv('data/DTarray_long.tsv', sep = '\t', stringsAsFactors = F)
dat <- dat[dat$Sequence != 'PROT_SEQ_NOT_FOUND',]

dat$nsaf_scaled <- scales::rescale(dat$nsaf,c(0,100), from = range(dat$nsaf, na.rm = T, finite = T))

#p <- ggplot(dat, aes(x = ))


#tryptic digest of database

