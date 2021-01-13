# Gene/Protein & compound identifier converters in R

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

|                       | To: **Gene symbol** | To: **Uniprot**   | To: **Ensembl**     | To: **Refseq** |  To: **Gene ID**    |
|:---------------------:|:-------------------:|:-----------------:|:-------------------:|:--------------:|:-------------------:|
| From: **Gene symbol** |          X          | [Yes](#OGS_to_UP) | [Yes](#OGS_to_ENSG) |                | [Yes](#OGS_to_GID)  |
| From: **Uniprot**     | [Yes](#UP_to_OGS)   |          X        | [Yes](#UP_to_ENS)   |                | [Yes](#UP_to_GID)   |
| From: **Ensembl**     |                     | [Yes](#ENS_to_UP) |          X          |                |                     |
| From: **Refseq**      |                     | [Yes](#RS_to_UP)  |                     |        X       |                     |
| From: **Gene ID**     | [Yes](#GID_to_OGS)  | [Yes](#GID_to_UP) | [Yes](#GID_to_ENSG) |                |           X         |


## Table of Contents
* [1. Installation](#installation)
* [2. Sample data](#sample)
* [3. Gene/protein converters](#GPconverters)
    * [3.1. Gene symbol to Uniprot](#OGS_to_UP)
    * [3.2. Gene symbol to Ensembl](#OGS_to_ENSG)
    
    * [3.4. Uniprot to Gene symbol](#UP_to_OGS)
    * [3.5. Uniprot to Ensembl](#UP_to_ENS)
    * [3.6. Uniprot to NCBI Gene ID](#UP_to_GID)
    * [3.7. Ensembl to Uniprot](#ENS_to_UP)
    * [3.8. NCBI Gene ID to Gene symbol](#GID_to_OGS)
    * [3.9. NCBI Gene ID to Uniprot](#GID_to_UP)
    * [3.10. NCBI Gene ID to Ensembl](#GID_to_ENSG)
    * [3.11. RefSeq to Uniprot](#RS_to_UP)
* [4. Compound converters](#Cconverters)



## <a name="installation">1. Installation</a>

```
source("https://raw.githubusercontent.com/hmbotelho/convert_gene_compound_IDs/master/initialize.R")
```



## <a name="sample">2. Sample data</a>

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



## <a name="GPconverters">3. Gene/protein converters</a>

The following sections describe and demonstrate the scripts converting between gene/protein identifiers.



## <a name="OGS_to_UP">3.1. Gene symbol to Uniprot</a>

**Description:** Converts human gene symbols to Uniprot IDs.  

**Input:**  
    * `genes`: character vector, with human gene symbols.  
    * `asvector`: logical, output as character? Otherwise, the output will be a data frame.  
**Output:**  
    * data.frame or character vector  
**Dependencies:** `biomaRt`  

```
genesymbol_to_uniprot(somegenesymbols, asvector=FALSE)
```
```
  hgnc_symbol uniprotswissprot
1        CFTR           P13569
2        CRYZ           Q08257
3      TGFBR3           Q03167
4       PI4KB           Q9UBF8
5         TAZ           Q16635
6      ZNF695           Q8IW36
7       BRCA1           P38398
```

```
genesymbol_to_uniprot(somegenesymbols)
```
```
    CFTR     CRYZ   TGFBR3    PI4KB      TAZ   ZNF695    BRCA1 
"P13569" "Q08257" "Q03167" "Q9UBF8" "Q16635" "Q8IW36" "P38398" 
```



## <a name="OGS_to_ENSG">3.2. Gene symbol to Ensembl</a>

**Description:** Converts human gene symbols to Ensembl gene IDs.  

**Input:**  
    * `genes`: character vector, with human gene symbols  
    * `asvector`: logical, output as character? Otherwise, the output will be a data frame.  
**Output:**  
    * data.frame or character vector  
**Dependencies:** `biomaRt`  

```
genesymbol_to_ensemblgene(somegenesymbols, asvector=FALSE)
```
```
  hgnc_symbol ensembl_gene_id
1        CFTR ENSG00000001626
2        CRYZ ENSG00000116791
3      TGFBR3 ENSG00000069702
4       PI4KB ENSG00000143393
5         TAZ ENSG00000102125
6      ZNF695 ENSG00000197472
7       BRCA1 ENSG00000012048
```

```
genesymbol_to_ensemblgene(somegenesymbols)
```
```
             CFTR              CRYZ            TGFBR3             PI4KB               TAZ            ZNF695             BRCA1 
"ENSG00000001626" "ENSG00000116791" "ENSG00000069702" "ENSG00000143393" "ENSG00000102125" "ENSG00000197472" "ENSG00000012048" 
```



## <a name="OGS_to_GID">3.3. Gene symbol to NCBI Gene ID</a>

**Description:** Converts human gene symbols to NCBI gene IDs.  

**Input:**  
    * `genes`: character vector, with human gene symbols  
    * `asvector`: logical, output as character? Otherwise, the output will be a data frame.  
**Output:**  
    * data.frame or character vector  
**Dependencies:** `biomaRt`  

```
genesymbol_to_geneid(somegenesymbols, asvector=FALSE)
```
```
  hgnc_symbol entrezgene_id
1        CFTR          1080
2        CRYZ          1429
3      TGFBR3          7049
4       PI4KB          5298
5         TAZ          6901
6      ZNF695         57116
7       BRCA1           672
```

```
genesymbol_to_geneid(somegenesymbols)
```
```
  CFTR   CRYZ TGFBR3  PI4KB    TAZ ZNF695  BRCA1 
  1080   1429   7049   5298   6901  57116    672 
```



## <a name="UP_to_OGS">3.4. Uniprot to Gene symbol</a>

**Description:** Converts Uniprot IDs to human gene symbols.  

**Input:**  
    * `up`: character vector, with Uniprot IDs  
	* `asvector`: logical, output as character? Otherwise, the output will be a data frame.  
**Output:**  
    * data.frame or character vector  
**Dependencies:** `biomaRt`

```
uniprot_to_genesymbol(someuniprot, asvector=FALSE)
```
```
  uniprotswissprot hgnc_symbol
1           P13569        CFTR
2           Q08257        CRYZ
3           Q03167      TGFBR3
4           Q9UBF8       PI4KB
5           Q16635         TAZ
6           Q8IW36      ZNF695
7           P38398       BRCA1
```

```
uniprot_to_genesymbol(someuniprot)
```
```
  P13569   Q08257   Q03167   Q9UBF8   Q16635   Q8IW36   P38398 
  "CFTR"   "CRYZ" "TGFBR3"  "PI4KB"    "TAZ" "ZNF695"  "BRCA1" 
```



## <a name="UP_to_ENS">3.5. Uniprot to Ensembl</a>

**Description:** Converts Uniprot IDs to Ensembl gene IDs.  

**Input:**  
    * `up`: character vector, with Uniprot IDs  
	* `asvector`: logical, output as character? Otherwise, the output will be a data frame.  
**Output:**  
    * data.frame or character vector  
**Dependencies:** `biomaRt`

```
uniprot_to_ensembl(someuniprot, asvector=FALSE)
```
```
  uniprotswissprot ensembl_gene_id
1           P13569 ENSG00000001626
2           Q08257 ENSG00000116791
3           Q03167 ENSG00000069702
4           Q9UBF8 ENSG00000143393
5           Q16635 ENSG00000102125
6           Q8IW36 ENSG00000197472
7           P38398 ENSG00000012048
```

```
uniprot_to_ensembl(someuniprot)
```
```
           P13569            Q08257            Q03167            Q9UBF8            Q16635            Q8IW36            P38398 
"ENSG00000001626" "ENSG00000116791" "ENSG00000069702" "ENSG00000143393" "ENSG00000102125" "ENSG00000197472" "ENSG00000012048" 
```



## <a name="UP_to_GID">3.6. Uniprot to Gene ID</a>

**Description:** Converts Uniprot IDs to NCBI Gene IDs.  

**Input:**  
    * `up`: character vector, with Uniprot IDs  
	* `asvector`: logical, output as character? Otherwise, the output will be a data frame.  
**Output:**  
    * data.frame or character vector  
**Dependencies:** `biomaRt`

```
uniprot_to_geneid(someuniprot, asvector=FALSE)
```
```
  uniprotswissprot entrezgene_id
1           P13569          1080
2           Q08257          1429
3           Q03167          7049
4           Q9UBF8          5298
5           Q16635          6901
6           Q8IW36         57116
7           P38398           672
```

```
uniprot_to_geneid(someuniprot)
```
```
P13569 Q08257 Q03167 Q9UBF8 Q16635 Q8IW36 P38398 
  1080   1429   7049   5298   6901  57116    672
```



## <a name="ENS_to_UP">3.7. Ensembl to Uniprot</a>

**Description:** Converts human ensembl gene symbols to Uniprot IDs.  

**Input:**  
    * `ens`: character vector, with human ensembl gene IDs  
	* `asvector`: logical, output as character? Otherwise, the output will be a data frame.  
**Output:**  
    * data.frame or character vector  
**Dependencies:** `rentrez`  

```
ensembl_to_uniprot(someensembl, asvector=FALSE)
```
```
  ensembl_gene_id uniprotswissprot
1 ENSG00000001626           P13569
2 ENSG00000116791           Q08257
3 ENSG00000069702           Q03167
4 ENSG00000143393           Q9UBF8
5 ENSG00000102125           Q16635
6 ENSG00000197472           Q8IW36
7 ENSG00000012048           P38398
```

```
ensembl_to_uniprot(someensembl)
```
```
ENSG00000001626 ENSG00000116791 ENSG00000069702 ENSG00000143393 ENSG00000102125 ENSG00000197472 ENSG00000012048 
       "P13569"        "Q08257"        "Q03167"        "Q9UBF8"        "Q16635"        "Q8IW36"        "P38398" 
```



## <a name="GID_to_OGS">3.8. NCBI Gene ID to Gene symbol</a>

**Description:** Converts NCBI gene IDs to Gene symbols.  

**Input:**  
    * `geneids`: vector (character or numeric), with NCBI gene IDs  
	* `asvector`: logical, output as character? Otherwise, the output will be a data frame.  
**Output:**  
    * data.frame or character vector  
**Dependencies:** `biomaRt`, `rentrez`  

```
geneid_to_genesymbol(somegeneids, asvector=FALSE)
```
```
  entrezgene_id hgnc_symbol
1          1080        CFTR
2          1429        CRYZ
3          7049      TGFBR3
4          5298       PI4KB
5          6901         TAZ
6         57116      ZNF695
7           672       BRCA1
```

```
geneid_to_genesymbol(somegeneids)
```
```
    1080     1429     7049     5298     6901    57116      672 
  "CFTR"   "CRYZ" "TGFBR3"  "PI4KB"    "TAZ" "ZNF695"  "BRCA1" 
```



## <a name="GID_to_UP">3.9. NCBI Gene ID to Uniprot</a>

**Description:** Converts NCBI gene IDs to Uniprot IDs.  

**Input:**  
    * `geneids`: vector (character or numeric), with NCBI gene IDs  
	* `asvector`: logical, output as character? Otherwise, the output will be a data frame.  
**Output:**  
    * data.frame or character vector  
**Dependencies:** `biomaRt`, `rentrez`  

```
geneid_to_uniprot(somegeneids, asvector=FALSE)
  
```
```
  entrezgene_id uniprotswissprot
1          1080           P13569
2          1429           Q08257
3          7049           Q03167
4          5298           Q9UBF8
5          6901           Q16635
6         57116           Q8IW36
7           672           P38398
```

```
geneid_to_uniprot(somegeneids)
```
```
    1080     1429     7049     5298     6901    57116      672 
"P13569" "Q08257" "Q03167" "Q9UBF8" "Q16635" "Q8IW36" "P38398" 
```



## <a name="GID_to_ENSG">3.10. NCBI Gene ID to Ensembl gene IDs</a>

**Description:** Converts NCBI gene IDs to Ensembl gene IDs.  

**Input:**  
    * `geneids`: vector (character or numeric), with NCBI gene IDs  
	* `asvector`: logical, output as character? Otherwise, the output will be a data frame.  
**Output:**  
    * data.frame or character vector  
**Dependencies:** `biomaRt`, `rentrez`  

```
geneid_to_ensemblgene(somegeneids, asvector=FALSE)
```
```
  entrezgene_id ensembl_gene_id
1          1080 ENSG00000001626
2          1429 ENSG00000116791
3          7049 ENSG00000069702
4          5298 ENSG00000143393
5          6901 ENSG00000102125
6         57116 ENSG00000197472
7           672 ENSG00000012048
```

```
geneid_to_ensemblgene(somegeneids)
```
```
             1080              1429              7049              5298              6901             57116               672 
"ENSG00000001626" "ENSG00000116791" "ENSG00000069702" "ENSG00000143393" "ENSG00000102125" "ENSG00000197472" "ENSG00000012048"
```



## <a name="RS_to_UP">3.11. RefSeq to Uniprot</a>

**Description:** NCBI RefSeq IDs to Uniprot IDs.  

**Input:**  
    * `refseqs`: character vector, with RefSeq protein IDs  
	* `asvector`: logical, output as character? Otherwise, the output will be a data frame.  
**Output:**  
    * data.frame or character vector  
**Dependencies:** `biomaRt`, `rentrez`  

```
refseq_to_uniprot(somerefseq, asvector=FALSE)
```
```
          refseq uniprotswissprot
1    NP_000483.2           P13569
2 NP_001123514.1           Q08257
3    NP_003234.2           Q03167
4 NP_001185702.1           Q9UBF8
5    NP_000107.1           Q16635
6    NP_009225.1           P38398
```

```
refseq_to_uniprot(somerefseq)
```
```
   NP_000483.2 NP_001123514.1    NP_003234.2 NP_001185702.1    NP_000107.1    NP_009225.1 
      "P13569"       "Q08257"       "Q03167"       "Q9UBF8"       "Q16635"       "P38398" 
```



## <a name="Cconverters">4. Compound converters</a>

The following sections describe and demonstrate the scripts converting between compound identifiers.

Coming soon...
