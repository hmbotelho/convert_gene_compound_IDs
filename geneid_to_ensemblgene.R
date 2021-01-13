# geneid_to_ensemblgene
#
# Converts NCBI gene IDs to Ensembl gene IDs.
# 
# Hugo Botelho & André Falcão
# v0.1
# 13 January 2021
#
# Input: 
#	* geneids: vector (character or numeric), with NCBI gene IDs
#   * asvector: output as character vector? Othwerwise, the output will be a data.frame
# Output: data frame or character, with the input and corresponding Ensembl gene IDs
geneid_to_ensemblgene <- function(geneids, asvector=TRUE){
    
    output<-getBM(attributes = c('entrezgene_id', 'ensembl_gene_id'),
                  filters = 'entrezgene_id', values = geneids, mart = myMart,
                  useCache = FALSE)
    
    output<-output[nchar(output$ensembl_gene_id)>0,]
    output <- output[match(geneids, output$entrezgene_id),]
    output$entrezgene_id <- geneids
    rownames(output) <- 1:nrow(output)
    
    if(asvector){
        output <- output$ensembl_gene_id
        names(output) <- geneids
    }
    
    output
}

