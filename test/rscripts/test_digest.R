
library(peptideUtils)
library(dplyr)

#read and format actual data
dat_path <- "/Users/Aaron/Documents/School_Work/Mass_spec_data/Wolozin_lab_Biotin-PG_samples/AD_and_FTD/data/peptideList_long.tsv"
dat <- read.csv(dat_path, sep = '\t', stringsAsFactors = F)

dat.filtered <- dat[dat$Unique == 1,]
dat.filtered <- dat.filtered[!(grepl('^reverse_', dat.filtered$Protein_ID)),]
dat.filtered <- dplyr::select(dat.filtered, id = Protein_ID, seq = Sequence, sc = 'Spectral_counts')

dat.filtered <- dat.filtered %>% dplyr::group_by(id, seq) %>%
  dplyr::summarise(sc = sum(sc)) %>%
  dplyr::ungroup()

dat.filtered$mc <- stringr::str_count(dat.filtered$seq, '([RK])(?=[^P])')

#make theoretical data
#ids <- unique(dat.filtered$id)
#ids <- read.csv('data/ids.tsv', sep = '\t', stringsAsFactors = F)
#this should work but it crashes rstudio
#sequences <- peptideUtils::getSequences(ids, 'data/UniProt_Human_Cravattlab_nonredundant2_98id_11-05-2012_reversed.fasta')

sequences <- read.csv('data/sequences.tsv', sep = '\t', stringsAsFactors = F)

#get tryptic peptides
peptides.all <- dplyr::select(utils::stack(peptideUtils::digest(sequences$seq, sequences$id, nMissedCleavages = 1)),
                          id = ind, seq = values)
peptides.all$id <- as.character(peptides.all$id)
#select only unique sequences
peptides.unique <- peptides.all %>% dplyr::group_by(seq) %>%
  dplyr::mutate(seq_count = n()) %>%
  dplyr::ungroup()
peptides.unique <- peptides.unique[peptides.unique$seq_count == 1,]
peptides.unique <- dplyr::select(peptides.unique, -seq_count)
peptides_temp <- peptides.unique[peptides.unique$id %in% dat.filtered$id,]

dat.joined <- dplyr::full_join(peptides_temp, dat.filtered, by = c('id', 'seq'))

dat.bad <- dplyr::anti_join(dat.filtered, peptides_temp, by = c('id', 'seq'))
seq_temp <- getSequences(dat.bad$id, 'data/UniProt_Human_Cravattlab_nonredundant2_98id_11-05-2012_reversed.fasta')
dat.bad$before <- nBefore(dat.bad$seq, seq_temp, 1)
dat.bad$after <- nAfter(dat.bad$seq, seq_temp, 1)
dat.bad$last <- substring(dat.bad$seq, nchar(dat.bad$seq))
rm(seq_temp)
dat.bad <- dplyr::left_join(dat.bad, unique(dplyr::select(dat, fs = Full_sequence, id = Protein_ID, seq = Sequence)), by = c('id', 'seq'))
dat.bad <- tidyr::extract(dat.bad, fs, into = c('before_test', 'after_test'), '^([A-Z-])\\.[A-Z]+\\.([A-Z-])$', remove = F)
dat.bad[dat.bad$before_test == '-',]$before_test <- ""
dat.bad[dat.bad$after_test == '-',]$after_test <- ""
dat.bad$good <- (dat.bad$before == dat.bad$before_test & dat.bad$after == dat.bad$after)
if(!all(dat.bad$good)){
  stop("!all(dat.bad$good))")
}

dat.bad <- dat.bad %>% dplyr::select(id, seq, sc, mc, before, after, last)
dat.bad$in_peptides <- F
for(i in nrow(dat.bad)){
  dat.bad$in_peptides[i] <- dat.bad$seq[i] %in% peptides.all[peptides.all$id == dat.bad$id[i],]
}

#stat filtering by explanations
dat.bad <- dat.bad[dat.bad$mc <= 1,] #gt 1 mc
#Remove non trypsin cleavages
dat.bad <- dat.bad[grepl('[RK]{2}', paste0(dat.bad$before, dat.bad$after)),]
dat.bad <- dat.bad[substr(dat.bad$after, 1, 1) != 'P',]

#seq_temp <- getSequences(dat.bad$id, 'data/UniProt_Human_Cravattlab_nonredundant2_98id_11-05-2012_reversed.fasta')
#dat.bad$two_after <- nAfter(dat.bad$seq, seq_temp, 2)

