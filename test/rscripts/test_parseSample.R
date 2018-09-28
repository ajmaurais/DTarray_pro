
library(Rcpp)

Rcpp::sourceCpp('rscripts/test_parseSample.cpp')

sampleNames <- c('Biotin-PG_Tryp_SF_RA_P2_1',
                 'Biotin-PG_Tryp_SF_Healthy_P3_1',
                 'Biotin-PG_Tryp_SF_Healthy_P1_2',
                 'Biotin-PG_Tryp_SF_Healthy_P2_3',
                 'Biotin-PG_Tryp_SF_Healthy_P4_2',
                 'Biotin-PG_Tryp_SF_RA_P3_3',
                 'Biotin-PG_Tryp_SF_RA_P4_1',
                 'Biotin-PG_Tryp_SF_Healthy_P1_3',
                 'Biotin-PG_Tryp_SF_Healthy_P2_2',
                 'Biotin-PG_Tryp_SF_Healthy_P4_3',
                 'Biotin-PG_Tryp_SF_RA_P3_2',
                 'Biotin-PG_Tryp_SF_RA_P1_1',
                 'Biotin-PG_Tryp_SF_RA_P4_2',
                 'Biotin-PG_Tryp_SF_RA_P2_3',
                 'Biotin-PG_Tryp_SF_Healthy_P2_1',
                 'Biotin-PG_Tryp_SF_RA_P1_2',
                 'Biotin-PG_Tryp_SF_RA_P3_1',
                 'Biotin-PG_Tryp_SF_RA_P4_3',
                 'Biotin-PG_Tryp_SF_RA_P2_2',
                 'Biotin-PG_Tryp_SF_Healthy_P1_1',
                 'Biotin-PG_Tryp_SF_Healthy_P3_2',
                 'Biotin-PG_Tryp_SF_Healthy_P4_1',
                 'Biotin-PG_Tryp_SF_RA_P1_3')
prefix = 'P[0-9]_[0-9]$'

argCombinations <- expand.grid(parseSampleName = c(TRUE, FALSE), outputFormat = c(TRUE, FALSE), re = c(TRUE, FALSE))

dat <- data.frame(sampleName = character(0), prefix = character(0), 
                  parseSampleName = logical(0), outputFormat = logical(0), re = logical(0),
                  old = character(0), new = character(0), stringsAsFactors = FALSE)

for(i in 1:nrow(argCombinations))
{
  for(s in sampleNames){
    new <- data.frame(sampleName = s, prefix = prefix,
                      parseSampleName = argCombinations[1,"parseSampleName"], 
                      outputFormat = argCombinations[1,"outputFormat"], re = argCombinations[1,"re"],
                      old = runTest('o', s, prefix, argCombinations[1,"parseSampleName"],
                                    argCombinations[1,"outputFormat"],
                                    argCombinations[1,"re"]),
                      new = runTest('n', s, prefix, argCombinations[1,"parseSampleName"],
                                    argCombinations[1,"outputFormat"],
                                    argCombinations[1,"re"]),
                      stringsAsFactors = FALSE)
    dat <- rbind(dat, new)
  }
}

dat$eq <- dat$old == dat$new
allGood <- all(dat$eq)
