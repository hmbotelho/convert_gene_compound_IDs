# genesymbol_to_uniprot
#
# Converts human gene symbols to Uniprot IDs
# 
# Hugo Botelho & André Falcão
# v0.2
# 12 January 2021
#
# Input: 
#	* genes: character vector, with human gene symbols
#   * asvector: output as character vector? Othwerwise, the output will be a data.frame
# Output: data frame or character, with the input and corresponding Uniprot IDs
genesymbol_to_uniprot <- function(genes, asvector=TRUE) {

    output <- getBM(attributes = c('hgnc_symbol', 'uniprotswissprot'), 
                    filters = 'hgnc_symbol', values = genes, mart = myMart,
                    useCache = FALSE)
    output <- output[nchar(output$uniprotswissprot)>0,]
    output <- output[match(genes, output$hgnc_symbol),]
    output$hgnc_symbol <- genes
    rownames(output) <- 1:nrow(output)
    
    if(asvector){
        output <- output$uniprotswissprot
        names(output) <- genes
    }
    
    output
}
