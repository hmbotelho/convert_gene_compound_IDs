# Gene/Protein & compound identifier converters

Genes and proteins can be referred to by many gene identifiers:
* Official gene symbol (*e.g.* CFTR)
* Uniprot ID (*e.g.* P13569)
* Ensembl gene ID (*e.g.* ENSG00000001626)
* Entrez gene ID (*e.g.* 1080)
* ...

Similarly, chemical compounds may be referred to by many different identifiers:
* Name (*e.g.* Ibuprofen)
* Trade name
* CAS number: (*e.g.* 15687-27-1)
* ChEMBL ID (*e.g.* CHEMBL521)
* ZINC ID (*e.g.* ZINC2647)
* ChEBI ID (*e.g.* CHEBI:5855)
* ...

The variety of identifiers arises due to the variety of databases hosting this information, as well as the different scope (gene/protein, research/commercial compounds) of each identifier.

This repository contains R scripts which convert between several of these identifiers.
The biomaRt scripts were based on previous work by André Falcão (FCUL).

The following converters are currently available:

|                       | To: **Gene symbol** | To: **Uniprot**   | To: **Ensembl**     | To: **Refseq** |  To: **Gene ID** |
|:---------------------:|:-------------------:|:-----------------:|:-------------------:|:--------------:|:----------------:|
| From: **Gene symbol** |          X          | [Yes](#OGS_to_UP) | [Yes](#OGS_to_ENSG) |                |                  |
| From: **Uniprot**     | [Yes](#UP_to_OGS)   |          X        | [Yes](#UP_to_ENS)   |                |                  |
| From: **Ensembl**     |                     | [Yes](#ENS_to_UP) |          X          |                |                  |
| From: **Refseq**      |                     | [Yes](#RS_to_UP)  |                     |        X       |                  |
| From: **Gene ID**     | [Yes](#GID_to_OGS)  |                   |                     |                |          X       |


## Table of Contents
* [1. Sample data](#sample)
* [2. Gene/protein converters](#GPconverters)
    * [2.1. Gene symbol to Uniprot](#OGS_to_UP)
    * [2.2. Gene symbol to Ensembl](#OGS_to_ENSG)
    * [2.3. Uniprot to Gene symbol](#UP_to_OGS)
    * [2.4. Uniprot to Ensembl](#UP_to_ENS)
    * [2.5. Ensembl to Uniprot](#ENS_to_UP)
    * [2.6. NCBI Gene ID to Gene symbol](#GID_to_OGS)
    * [2.7. RefSeq to Uniprot](#RS_to_UP)
* [3. Compound converters](#Cconverters)



## <a name="sample">1. Sample data</a>

Let us use a few genes to show the converters in action.
```
## Gene Identifiers ####

# Gene symbols
somegenesymbols <- c("CFTR", "CRYZ", "TGFBR3", "PI4KB", "TAZ", "ZNF695", "BRCA1")

# Uniprot IDs
someuniprot <- c("P13569", "Q08257", "Q03167", "Q9UBF8", "Q16635", "Q8IW36", "P38398")

# Ensembl gene IDs
someensembl <- c("ENSG00000001626", "ENSG00000116791", "ENSG00000069702", "ENSG00000143393", "ENSG00000102125", "ENSG00000197472", "ENSG00000012048")

# NCBI gene IDs
somegeneids <- c(1080, 1429, 7049, 5298, 6901, 57116, 672)

# RefSeq protein IDs
somerefseq <- c("NP_000483.2", "NP_001123514.1", "NP_003234.2", "NP_001185702.1", "NP_000107.1", "NP_009225.1")


## Compound Identifiers ####

# Chembl IDs
somechembl <- c("CHEMBL521", "CHEMBL25", "CHEMBL165", "CHEMBL857", "CHEMBL2010601", "CHEMBL925", "CHEMBL112570")

# ZINC IDs
somezinc <- c("ZINC2647", "ZINC53", "ZINC6787", "ZINC35024346", "ZINC52509463", "ZINC266964", "ZINC3875383")
```



## <a name="GPconverters">2. Gene/protein converters</a>

The following sections describe and demonstrate the scripts converting between gene/protein identifiers.

```
source("genesymbol_to_ensemblgene.R")
```

In the following sections we will see these scripts in action.



## <a name="OGS_to_UP">2.1. Gene symbol to Uniprot</a>

**Description:** Converts human gene symbols to Uniprot IDs.  

**Input:**  
    * `genes`: character vector, with human gene symbols.  
    * `showProgress`: logical, print progress to console?  
	* `parallelize`: logical, use multi-processors?  
**Output:**  
    * data frame or character vector  
**Dependencies:** `biomaRt`  

```
> source("genesymbol_to_uniprot.R")
> genesymbol_to_uniprot_df(somegenesymbols)
   hgnc_symbol uniprotswissprot
1          TAZ           Q16635
4         CFTR           P13569
5        BRCA1           P38398
8       TGFBR3           Q03167
10      ZNF695           Q8IW36
11        CRYZ           Q08257
14       PI4KB           Q9UBF8
>
> genesymbol_to_uniprot_chr(somegenesymbols, showProgress = FALSE, parallelize = TRUE)
[1] "P13569" "Q08257" "Q03167" "Q9UBF8" "Q16635" "Q8IW36" "P38398"
```



## <a name="OGS_to_ENSG">2.2. Gene symbol to Ensembl</a>

**Description:** Converts human gene symbols to Ensembl gene IDs.  

**Input:**  
    * `genes`: character vector, with human gene symbols  
    * `showProgress`: logical, print progress to console?  
	* `parallelize`: logical, use multi-processors?  
**Output:**  
    * data frame or character vector  
**Dependencies:** `biomaRt`, `parallel`  

```
> source("genesymbol_to_ensemblgene.R")
> genesymbol_to_ensemblgene_df(somegenesymbols)
  hgnc_symbol ensembl_gene_id
1       BRCA1 ENSG00000012048
2        CFTR ENSG00000001626
3        CRYZ ENSG00000116791
4       PI4KB ENSG00000143393
5         TAZ ENSG00000102125
6      TGFBR3 ENSG00000069702
7      ZNF695 ENSG00000197472
>
> genesymbol_to_ensemblgene_chr(somegenesymbols, showProgress = FALSE, parallelize = TRUE)
[1] "ENSG00000001626" "ENSG00000116791" "ENSG00000069702" "ENSG00000143393" "ENSG00000102125" "ENSG00000197472"
[7] "ENSG00000012048"
```



## <a name="UP_to_OGS">2.3. Uniprot to Gene symbol</a>

**Description:** Converts Uniprot IDs to human gene symbols.  

**Input:**  
    * `up`: character vector, with Uniprot IDs  
	* `showProgress`: logical, print progress to console?  
	* `parallelize`: logical, use multi-processors?  
**Output:**  
    * character vector  
**Dependencies:** `biomaRt`, `parallel`  

```
> source("uniprot_to_genesymbol.R")
> uniprot_to_genesymbol_chr(someuniprot, showProgress = F, parallelize = T)
[1] "CFTR"   "CRYZ"   "TGFBR3" "PI4KB"  "TAZ"    "ZNF695" "BRCA1"
```



## <a name="UP_to_ENS">2.4. Uniprot to Ensembl</a>

**Description:** Converts Uniprot IDs to Ensembl gene IDs.  

**Input:**  
    * `up`: character vector, with Uniprot IDs  
	* `showProgress`: logical, print progress to console?  
	* `parallelize`: logical, use multi-processors?  
**Output:**  
    * character vector  
**Dependencies:** `biomaRt`, `parallel`  

```
> source("uniprot_to_ensembl.R")
> uniprot_to_ensembl_chr(someuniprot, showProgress = F, parallelize = T)
[1] "ENSG00000001626" "ENSG00000116791" "ENSG00000069702" "ENSG00000143393" "ENSG00000102125" "ENSG00000197472"
[7] "ENSG00000012048"
```



## <a name="ENS_to_UP">2.5. Ensembl to Uniprot</a>

**Description:** Converts human ensembl gene symbols to Uniprot IDs.  

**Input:**  
    * `ens`: character vector, with human ensembl gene IDs  
	* `showProgress`: logical, print progress to console?  
**Output:**  
    * character vector  
**Dependencies:** `rentrez`  

```
> source("ensembl_to_uniprot.R")
> ensembl_to_uniprot_chr(someensembl, showProgress = F, parallelize = T)
[1] "P13569" "Q08257" "Q03167" "Q9UBF8" "Q16635" "Q8IW36" "P38398"
```



## <a name="GID_to_OGS">2.6. NCBI Gene ID to Gene symbol</a>

**Description:** Converts NCBI gene IDs to Gene symbols.  

**Input:**  
    * `geneids`: vector (character or numeric), with NCBI gene IDs  
	* `showProgress`: logical, print progress to console?  
**Output:**  
    * character vector  
**Dependencies:** `rentrez`  

```
> source("geneid_to_genesymbol.R")
> geneid_to_genesymbol_chr(somegeneids, showProgress = F)
[1] "CFTR"   "CRYZ"   "TGFBR3" "PI4KB"  "TAZ"    "ZNF695" "BRCA1" 
```



## <a name="RS_to_UP">2.7. RefSeq to Uniprot</a>

**Description:** NCBI RefSeq IDs to Uniprot IDs.  

**Input:**  
    * `refseqs`: character vector, with RefSeq protein IDs  
	* `simplify`: logical, should "P13569.3" be converted to "P13569"?  
	* `showProgress`: logical, print progress to console?  
**Output:**  
    * character vector  
**Dependencies:** `rentrez`  

```
> source("refseq_to_uniprot.R")
> refseq_to_uniprot_chr(somerefseq, simplify = T, showProgress = F)
[1] "P13569" "Q08257" "Q03167" "Q9UBF8" "Q16635" "P38398"
```



## <a name="Cconverters">3. Compound converters</a>

The following sections describe and demonstrate the scripts converting between compound identifiers.

Coming soon...
