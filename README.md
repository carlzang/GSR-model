The files contain an R code for computing the GSR indexes, testing datasets using 85 ovarian clear cell carcinoma and 136 normal ovarian tissue gene expression profiles downloaded from GEO datasets in .SOFT format. The sample information was list in the caselist.csv. 

To compute GSR indexs, download the following files to you R working dictionary :

1. go.gmt :  'c5.all.v5.0.symbols.gmt' downloaeded from MSigDB and rename to 'go.gmt' for GO gene set defintions

2. pathway.gmt : ' c2.cp.reactome.v5.0.symbols.gmt' ,downloaeded from MSigDB and rename to 'pathway.gmt' for Reactome pathway gene set defintions

3. names_common : one-column dataframe, gene symbols from the genes of all microarray gene datasets in common

3. conData : gene expression matrix, the row is the normal control samples, the column is the gene symbols of names_common

4. disData : gene expression matrix, the row is the disease samples, the column is the gene symbols of names_common

The SOFT files from different microarray platforms should be intersected to find the genes in common (names_common), only the common genes and the corresponding expression profiles are utilized. 
The outfiles are
'Rvalue_con_GO_OvcaClearCell' = GSR indexes of cotntrol group computing through GO term gene sets ;
'Rvalue_dis_GO_OvcaClearCell' = GSR indexes of clear cell ca group computing through GO terms gene sets ; 
'Rvalue_con_CAN_OvcaClearCell' = GSR indexes of cotntrol group computing through canonical (Reactome) pathway gene sets ;
'Rvalue_dis_CAN_OvcaClearCell' = GSR indexes of clear cell ca group computing through canonical (Reactome) pathway gene sets .
