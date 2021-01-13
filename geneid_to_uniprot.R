# geneid_to_uniprot
#
# Converts NCBI gene IDs to Uniprot IDs.
# 
# Hugo Botelho & André Falcão
# v0.1
# 13 January 2021
#
# Input: 
#	* geneids: vector (character or numeric), with NCBI gene IDs
#   * asvector: output as character vector? Othwerwise, the output will be a data.frame
# Output: data frame or character, with the input and corresponding Uniprot IDs
geneid_to_uniprot <- function(geneids, asvector=TRUE){
    
    output<-getBM(attributes = c('entrezgene_id', 'uniprotswissprot'),
                  filters = 'entrezgene_id', values = geneids, mart = myMart,
                  useCache = FALSE)
    
    output<-output[nchar(output$uniprotswissprot)>0,]
    output <- output[match(geneids, output$entrezgene_id),]
    output$entrezgene_id <- geneids
    rownames(output) <- 1:nrow(output)
    
    if(asvector){
        output <- output$uniprotswissprot
        names(output) <- geneids
    }
    
    output
}

