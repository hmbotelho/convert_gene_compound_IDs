# genesymbol_to_ensemblgene
#
# Converts human gene symbols to Ensembl gene IDs
# 
# Hugo Botelho & André Falcão
# v0.2
# 12 January 2021
#
# Input: 
#	* genes: character vector, with human gene symbols
#   * asvector: output as character vector? Othwerwise, the output will be a data.frame
# Output: data frame or character, with the input and corresponding Ensembl gene IDs
genesymbol_to_ensemblgene <- function(genes, asvector=TRUE) {

    output<-getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'), 
                  filters = 'hgnc_symbol', values = genes, mart = myMart,
                  useCache = FALSE)
    output<-output[nchar(output$ensembl_gene_id)>0,]
    output <- output[match(genes, output$hgnc_symbol),]
    output$hgnc_symbol <- genes
    rownames(output) <- 1:nrow(output)
    
    if(asvector){
        output <- output$ensembl_gene_id
        names(output) <- genes
    }
    
    output
}
