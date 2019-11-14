library(RCurl)
library(openxlsx)


# Works with a single name
ChEMBLtoZINC <- function(chemblID){
    
    zinc <- sapply(chemblID, function(x){

        chemblURL <- paste0("https://www.ebi.ac.uk/chembl/compound/inspect/", x)

        if(!url.exists(chemblURL)) return("")
        print(paste0("Fetching ZINC id for ", x), quote = FALSE)
        
        # download html
        html <- getURL(chemblURL, followlocation = TRUE)
        
        # get zinc id
        regex <- ".*<a href=\"http://zinc15.docking.org/substances/.*?\">(.*?)</a>.*"
        if(!grepl(regex, html)) return("")
        
        gsub(regex, "\\1", html)
        
    })
    
    zinc
}