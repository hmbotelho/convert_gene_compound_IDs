# uniprot_to_geneid
#
# Converts Uniprot IDs to human gene symbols. Parallel processing available
# 
# Hugo Botelho & André Falcão
# v0.2
# 13 January 2021
#
# Input: 
#	* up: character vector, with Uniprot IDs
#   * asvector: output as character vector? Othwerwise, the output will be a data.frame
# Output: data frame or character, with the input and corresponding entrez gene IDs
uniprot_to_geneid <- function(up, asvector=TRUE){
    
    output <- getBM(attributes = c('uniprotswissprot', 'entrezgene_id'),
                    filters = 'uniprotswissprot', values = up, mart = myMart,
                    useCache = FALSE)
    output <- output[nchar(output$entrezgene_id)>0,]
    output <- output[match(up, output$uniprotswissprot),]
    output$uniprotswissprot <- up
    rownames(output) <- 1:nrow(output)
    
    if(asvector){
        output <- output$entrezgene_id
        names(output) <- up
    }
    
    output
}
