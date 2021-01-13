# genesymbol_to_geneid
#
# Converts human gene symbols to NCBI gene IDs
# 
# Hugo Botelho & André Falcão
# v0.2
# 12 January 2021
#
# Input: 
#	* genes: character vector, with human gene symbols
#   * asvector: output as character vector? Othwerwise, the output will be a data.frame
# Output: data frame or character, with the input and corresponding Uniprot IDs
genesymbol_to_geneid <- function(genes, asvector=TRUE) {

    output <- getBM(attributes = c('hgnc_symbol', 'entrezgene_id'), 
                    filters = 'hgnc_symbol', values = genes, mart = myMart,
                    useCache = FALSE)
    output <- output[nchar(output$entrezgene_id)>0,]
    output <- output[match(genes, output$hgnc_symbol),]
    output$hgnc_symbol <- genes
    rownames(output) <- 1:nrow(output)
    
    if(asvector){
        output <- output$entrezgene_id
        names(output) <- genes
    }
    
    output
}
