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




Computating GSR indexes by modified Differential Rank Conservation (DIRAC) : 

The algorithm of computing GSR indexes was modified from DIRAC. The DIRAC measures the perturbation of a gene set by converting gene expression levels to gene expression rankings, and quantifying the regularity of gene expression ranking in the gene set. This ranking is a binary vector containing all possible combinations of one-to-one gene relationships established by pairwise comparison of relative gene expression levels. The algorithm of DIRAC starts with establishing a gene set ranking template of each gene set for either disease or unaffected control group. The next step of DIRAC computes the ranking matching score, which measures the degree of each sample's gene expression ranking of each gene set matching the corresponding gene set ranking template. The ranking matching score of DIRAC measures the perturbation of gene expression ranking in comparison with the most common gene expression ranking in a gene set for either disease or unaffected control group. 
Instead of measuring the perturbation of gene expression ranking, the GSR index measures the change of gene expression ranking between two phenotypes in a gene set. For this purpose, the GSR indexes of both EOC and the normal control are computed by comparing the sample's gene expression ranking with a standard template derived from the most common gene expression ranking in a gene set among the entire normal ovarian tissue control samples. Then the subsequent analyses were carried out based on this same standard with the EOC and normal control GSR indexes. We define a baseline gene set ranking template as a template of the most common gene ranking among the unaffected controls in a gene set; it is used as a standard template for a gene set from the unaffected population. A gene set contains m gene G = {G1, …, Gm}, and the corresponding gene expression profile E = (E1, …, Em), Ei denotes the expression level of gene Gi. Each sample is labeled by a phenotype P ∈ {D, C}, where D or C denotes case (EOC) and unaffected control group, respectively. The baseline gene set raking template for each gene set is established by pairwise comparison between the expression levels of two genes for all possible combinations of gene pair. The baseline gene rank template B for a given gene set G is the binary vector ('A' or 'B') where each component is either 'A' if the probabilities Pr(Ei > Ej |P = C)>0.5 or 'B' if Pr(Ei > Ej |P = C)≤0.5. For the expression profile of a given sample en, the GSR index for a given gene set is the fraction of the m*(m–1)/2 pairs for which the observed gene expression ranking within en matches the baseline gene ranking template B.
