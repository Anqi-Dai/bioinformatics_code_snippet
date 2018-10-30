# This snippet is about using biomaRt package to retrieve ensembl feature information 
# for detailed documentation, please go to https://www.bioconductor.org/packages//2.7/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
library(biomaRt)
library(tidyverse)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# the ensembl_gene_id of the genes that you are interested in 
genes <- RF14.ctrl$Gene

# retrieve the feature information of the listed attributes
symbol <- getBM(filters = "ensembl_gene_id",
                attributes = c('ensembl_gene_id','hgnc_symbol','external_gene_name',
                               'chromosome_name','start_position','end_position', 'description'),
                
                values = genes, 
                mart = mart)

# could write out for easier future use
write_csv(symbol, '../data/feature_info_for_SDE2.csv')

# when reading this dataframe it has to be something like (define the col_types so that the chromosome name will not be read as integers)
featureInfo <-  read_csv('../data/feature_info_for_SDE2.csv', col_types = 'cccciic')
