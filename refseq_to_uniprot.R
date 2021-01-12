source("https://raw.githubusercontent.com/hmbotelho/convert_gene_compound_IDs/master/initialize.R")

# refseq_to_uniprot
#
# Converts refseq protein IDs to Uniprot IDs.
# 
# Hugo Botelho
# v0.2
# 12 January 2021
#
# Input: 
#	* refseqs: character vector, with RefSeq protein IDs
#   * asvector: output as character vector? Othwerwise, the output will be a data.frame
#	* simplify: logical, should "P13569.3" be converted to "P13569"?
#	* showProgress: logical, print progress to console?
# Output: data frame or character, with the input and corresponding Uniprot IDs
refseq_to_uniprot <- function(refseqs, asvector=TRUE, simplify = TRUE, showProgress = FALSE){
    
    output <- lapply(refseqs, function(REFSEQ){
        
        if(showProgress) print(paste0("Converting '", REFSEQ, "' to Uniprot ID"), quote=FALSE)
        
        # Check that the refseq identifier exists and is unambiguous
        preliminary_info <- entrez_search(db="protein", term = REFSEQ)
        if(preliminary_info$count != 1) return(NULL)
        
        
        # Get the most updated record for this gene
        summary1 <- entrez_summary(db="protein", id=REFSEQ)
        if("replacedby" %in% names(summary1)){
            refseq_updated <- summary1$replacedby
        } else{
            refseq_updated <- REFSEQ
        }
        
        
        # Get all Uniprot IDs from the record
        ncbi_data <- entrez_fetch(db="protein", id=refseq_updated, rettype="raw")
        
        regex_uniprot <- "UniProtKB/Swiss-Prot \\((.*?)\\)"
        uniprot_loc0 <- gregexpr(regex_uniprot, ncbi_data)
        uniprot_loc  <- data.frame(start = uniprot_loc0[[1]],
                                   end   = uniprot_loc0[[1]] + attr(uniprot_loc0[[1]], which="match.length") - 1)
        uniprot_IDs <- apply(uniprot_loc, 1, function(x){
            sub(regex_uniprot, "\\1", substring(ncbi_data, x["start"], x["end"]))
        })
        
        
        # Select the most likely ID (i.e. the most abundant)
        id_table <- table(uniprot_IDs)
        Uniprot <- names(which(id_table == max(id_table))[1])
        
        
        if(simplify){
            # Remove dots (P13569.3  ->  P13569)
            Uniprot <- sub("\\..*", "", Uniprot)
        }
        
        c(REFSEQ, Uniprot)
        
    })
    
    output <- do.call("rbind", output)
    output <- as.data.frame(output, stringsAsFactors = FALSE)
    colnames(output) <- c("refseq", "uniprotswissprot")
    class(output$refseq) <- "character"
    class(output$uniprotswissprot) <- "character"
    output <- output[match(refseqs, output$refseq),]
    output$refseq <- refseqs
    rownames(output) <- 1:nrow(output)
    
    if(asvector){
        output <- output$uniprotswissprot
        names(output) <- refseqs
    }
    
    output
}
