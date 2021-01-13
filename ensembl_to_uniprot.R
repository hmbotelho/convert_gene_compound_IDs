# ensembl_to_uniprot
#
# Converts human ensembl gene symbols to Uniprot IDs
# 
# Hugo Botelho & André Falcão
# v0.2
# 12 January 2021
#
# Input: 
#	* ens: character vector, with human ensembl gene IDs
#   * asvector: output as character vector? Othwerwise, the output will be a data.frame
# Output: data frame or character, with the input and corresponding Ensembl gene IDs
ensembl_to_uniprot <- function(ens, asvector=TRUE) {
    
    output<-getBM(attributes = c('ensembl_gene_id', 'uniprotswissprot'),
                  filters = 'ensembl_gene_id', values = ens, mart = myMart,
                  useCache = FALSE)
    output<-output[nchar(output$uniprotswissprot)>0,]
    output <- output[match(ens, output$ensembl_gene_id),]
    output$ensembl_gene_id <- ens
    rownames(output) <- 1:nrow(output)
    
    if(asvector){
        output <- output$uniprotswissprot
        names(output) <- ens
    }
    
    output
}

