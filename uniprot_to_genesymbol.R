source("https://raw.githubusercontent.com/hmbotelho/convert_gene_compound_IDs/master/initialize.R")

# uniprot_to_genesymbol
#
# Converts Uniprot IDs to human gene symbols. Parallel processing available
# 
# Hugo Botelho & André Falcão
# v0.2
# 12 January 2021
#
# Input: 
#	* up: character vector, with Uniprot IDs
#   * asvector: output as character vector? Othwerwise, the output will be a data.frame
# Output: data frame or character, with the input and corresponding human gene symbols
uniprot_to_genesymbol <- function(up, asvector=TRUE){
    
    output <- getBM(attributes = c('uniprotswissprot', 'hgnc_symbol'),
                    filters = 'uniprotswissprot', values = up, mart = myMart,
                    useCache = FALSE)
    output <- output[nchar(output$hgnc_symbol)>0,]
    output <- output[match(up, output$uniprotswissprot),]
    output$uniprotswissprot <- up
    rownames(output) <- 1:nrow(output)
    
    if(asvector){
        output <- output$hgnc_symbol
        names(output) <- up
    }
    
    output
}
