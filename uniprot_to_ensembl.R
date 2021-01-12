source("https://raw.githubusercontent.com/hmbotelho/convert_gene_compound_IDs/master/initialize.R")

# uniprot_to_ensembl
#
# Converts Uniprot IDs to Ensembl Gene IDs.
# 
# Hugo Botelho & André Falcão
# v0.2
# 12 January 2021
#
# Input: 
#	* genes: character vector, with Uniprot IDs
# Output: data frame, with the input and corresponding Ensembl gene IDs
uniprot_to_ensembl <- function(up, asvector=TRUE) {
    
    output<-getBM(attributes = c('uniprotswissprot', 'ensembl_gene_id'),
                         filters = 'uniprotswissprot', values = up, mart = myMart,
                         useCache = FALSE)
    output<-output[nchar(output$ensembl_gene_id)>0,]
    output <- output[match(up, output$uniprotswissprot),]
    output$uniprotswissprot <- up
    rownames(output) <- 1:nrow(output)
    
    if(asvector){
        output <- output$ensembl_gene_id
        names(output) <- up
    }
    
    output
}
