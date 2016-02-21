
### GSR model ####
### need files : in R working dictionary
#1. go.gmt :  'c5.all.v5.0.symbols.gmt' downloaeded from MSigDB and rename to 'go.gmt' for GO gene set defintions
#2. pathway.gmt : ' c2.cp.reactome.v5.0.symbols.gmt' ,downloaeded from MSigDB and rename to 'pathway.gmt' for Reactome pathway gene set defintions
#3. names_common : one-column dataframe, gene symbols from the genes of all microarray gene datasets in common
#3. conData : gene expression matrix, the row is the normal control samples, the column is the gene symbols of names_common
#4. disData : gene expression matrix, the row is the disease samples, the column is the gene symbols of names_common
# conData and disData are derived from the downloaded gene expression profile in SOFT format from GEO database
# the SOFT files from different microarray platforms should be intersected to find the genes in common (names_common), only the common genes 
# and the corresponding expression profiles are utilized. 
# outfile Rvalue_con = GSR indexes of cotntrol group; Rvalue_dis = GSR indexes of disease group
#########################################################

library(combinat)

disease='OvcaClearCell'    # assigning the disease name, ovarian clear cell carcinoma for the dataset 
annotation= 'GO'          # usign GO gene set for annotation 
file= paste('D:/working/DB/geneSet/','go.gmt', sep='')  # location of go.gmt file   
gs.number= nrow(read.csv(file, header=F))    # gene set number 

#########################################################
##  preparing data 

con=read.table('conData', sep=' ')    
dis=read.table('disData', sep=' ')     
gene= read.table('names_common') 

condata= data.frame(gene, con)
condata= condata[order(condata$V1),]
disdata= data.frame(gene, dis)
disdata= disdata[order(disdata$V1),]

geneset= list()
gs.name.t= vector()

for ( i in 1:gs.number){
  geneset[[i]]= scan(file, sep='\t',  what='character', quote=NULL, nlines=1, skip= i-1)
  gs.name= geneset[[i]][1]
  geneset[[i]]= sort(geneset[[i]][3:length(geneset[[i]])])
  gs.name.t= c(gs.name.t, gs.name) 
}
names(geneset)= gs.name.t

save(geneset, file=paste('load', annotation, 'geneset', sep='_'))

########################################################
##  computing baseline gene set template 

load(paste('load', annotation, 'geneset', sep='_'))

GS= geneset
con_template= list()
gs.nullid= vector()

for ( i in 1:length(GS)) {
  
  GS.name= names(GS)[i]
  GS.gene= GS[[i]]
  GS.condata.id= which(is.element(condata$V1, GS.gene))
  GS.gene= GS.gene[is.element(GS.gene, condata$V1)]
  GS.data= condata[GS.condata.id,] 
  rownames(GS.data)= GS.data[,1]
  GS.data= t(GS.data)
  
  if (ncol(GS.data)< 3) {
    gs.nullid= c(gs.nullid, i)
    next }
  
  rownames(GS.data)=NULL
  GS.data=GS.data[2:nrow(GS.data),]
  
  combination= combn(length(GS.gene),2)
  combination.n= ncol(combination)
  
  if (is.null(combination.n)) {
    gs.nullid= c(gs.nullid, i)
    next }
  
  grank.ings= function (m){
    pair= GS.data[,combination[,m]]
    pair.compare= ifelse(pair[,1]>=pair[,2], 'A','B')
    A= length(pair.compare[pair.compare=='A'])
    B= length(pair.compare[pair.compare=='B'])
    pair.max= ifelse(A>=B, 'A','B')
    return(pair.max)
  }
  grank.gs= sapply(1:combination.n, grank.ings)
  cat('i= ', i , '\n')
  con_template[[i]]= grank.gs
}  
geneset.tot= names(GS)
if (length(con_template) != length(geneset.tot))   { geneset.tot= geneset.tot[-length(geneset.tot)] }
names(con_template)=geneset.tot
con_template.ori= con_template
if (length(gs.nullid)>0) {con_template= con_template[-gs.nullid]}
save(con_template, file='load_con_template_GO')
GS.ori= GS
if (length(gs.nullid)>0) {GS= GS[-gs.nullid]}
save(GS, file=paste('load_GS', annotation, sep='_'))

###########################################################
## compute control Gene Set Regularity index (GSR indexes, or Rvalue ) 
load('load_con_template_GO')
load(paste('load_GS', annotation, sep='_'))

casen= ncol(condata)-1
match.score.allcase.con= data.frame()

for ( k in 1: casen) {
  condata.k= condata[,k+1]
  
  matchscore.tot.onecase= vector()
  gsrank_onecase= function(i){       # i = ith geneset 
    
    gs.name= names(GS)[i]
    gs.gene= GS[[i]]
    gs.condata.id= which(is.element(condata$V1, gs.gene))
    gs.gene= gs.gene[is.element(gs.gene, condata$V1)]
    condata.k= condata.k[gs.condata.id]
    
    combination= combn(length(gs.gene), 2)
    combination.n= ncol(combination)
    
    rank.ings.onecase= vector()
    rank.ings= function(m) { 
      condata.gs.rank = ifelse(condata.k[combination[1,m]] >= condata.k[combination[2,m]], 'A','B')  
      rank.ings.onecase= c(rank.ings.onecase, condata.gs.rank)
      return(rank.ings.onecase)
    }
    rank.case= sapply(1:combination.n, rank.ings)
    
    matchone= ifelse(rank.case==con_template[[i]], 'y', 'n' )
    matchscore= length(matchone[matchone=='y'])/length(matchone)
    
    matchscore.tot.onecase= c(matchscore.tot.onecase, matchscore)
    return(matchscore.tot.onecase)
    
  }
  match.score= sapply(1:length(GS), gsrank_onecase)  
  cat('case number=  ', k, '\n')
  match.score.allcase.con= rbind(match.score.allcase.con, match.score)
  
}  
colnames(match.score.allcase.con)=names(GS)  
Rvalue.con_GO= match.score.allcase.con
write.csv(match.score.allcase.con, paste('Rvalue_con', annotation, disease,sep='_'))
write.csv(names(GS),paste('Geneset', annotation, disease,sep='_') )

#####################################################################
## compute disease Gene Set Regularity (GSR) index  

casen= ncol(disdata)-1
match.score.allcase.dis= data.frame()

for ( k in 1: casen) {
  disdata.k= disdata[,k+1]
  
  matchscore.tot.onecase= vector()
  gsrank_onecase= function(i){       # i = ith geneset 
    
    gs.name= names(GS)[i]
    gs.gene= GS[[i]]
    gs.disdata.id= which(is.element(disdata$V1, gs.gene))
    gs.gene= gs.gene[is.element(gs.gene, disdata$V1)]
    disdata.k= disdata.k[gs.disdata.id]
    
    combination= combn(length(gs.gene), 2)
    combination.n= ncol(combination)
    
    rank.ings.onecase= vector()
    rank.ings= function(m) { 
      disdata.gs.rank = ifelse(disdata.k[combination[1,m]] >= disdata.k[combination[2,m]], 'A','B')  
      rank.ings.onecase= c(rank.ings.onecase, disdata.gs.rank)
      return(rank.ings.onecase)
    }
    rank.case= sapply(1:combination.n, rank.ings)
    
    matchone= ifelse(rank.case==con_template[[i]], 'y', 'n' )
    matchscore= length(matchone[matchone=='y'])/length(matchone)
    
    matchscore.tot.onecase= c(matchscore.tot.onecase, matchscore)
    return(matchscore.tot.onecase)
    
  }
  match.score= sapply(1:length(GS), gsrank_onecase)  
  cat('case number=  ', k, '\n')
  match.score.allcase.dis= rbind(match.score.allcase.dis, match.score)
  
}  
colnames(match.score.allcase.dis)=names(GS)  
Rvalue.dis_GO= match.score.allcase.dis
write.csv(match.score.allcase.dis, paste('Rvalue_dis', annotation, disease,sep='_'))

###################################################################################################################
##################################################################################################################

annotation= 'CAN'     # using Reactome canonical pathway gene set 
file= file= paste('D:/working/DB/geneSet/','pathway.gmt', sep='')  # location of Reactome pathway gene set defintion 
gs.number= nrow(read.csv(file, header=F))

#################################################
## preparing data 

con=read.table('conData', sep=' ')
dis=read.table('disData', sep=' ')
gene= read.table('names_common')

condata= data.frame(gene, con)
condata= condata[order(condata$V1),]
disdata= data.frame(gene, dis)
disdata= disdata[order(disdata$V1),]

geneset= list()
gs.name.t= vector()

for ( i in 1:gs.number){
  geneset[[i]]= scan(file, sep='\t',  what='character', quote=NULL, nlines=1, skip= i-1)
  gs.name= geneset[[i]][1]
  geneset[[i]]= sort(geneset[[i]][3:length(geneset[[i]])])
  gs.name.t= c(gs.name.t, gs.name) 
}
names(geneset)= gs.name.t

########################################################
##  computing baseline gene set template 

GS= geneset
con_template= list()
gs.nullid= vector()

for ( i in 1:length(GS)) {
  
  GS.name= names(GS)[i]
  GS.gene= GS[[i]]
  GS.condata.id= which(is.element(condata$V1, GS.gene))
  GS.gene= GS.gene[is.element(GS.gene, condata$V1)]
  GS.data= condata[GS.condata.id,] 
  rownames(GS.data)= GS.data[,1]
  GS.data= t(GS.data)
  
  if (ncol(GS.data)< 3) {
    gs.nullid= c(gs.nullid, i)
    next }

  rownames(GS.data)=NULL
  GS.data=GS.data[2:nrow(GS.data),]
  
  combination= combn(length(GS.gene),2)
  combination.n= ncol(combination)
  
  if (is.null(combination.n)) {
    gs.nullid= c(gs.nullid, i)
    next }
  
  grank.ings= function (m){
    pair= GS.data[,combination[,m]]
    pair.compare= ifelse(pair[,1]>=pair[,2], 'A','B')
    A= length(pair.compare[pair.compare=='A'])
    B= length(pair.compare[pair.compare=='B'])
    pair.max= ifelse(A>=B, 'A','B')
    return(pair.max)
  }
  grank.gs= sapply(1:combination.n, grank.ings)
  cat('i= ', i , '\n')
  con_template[[i]]= grank.gs
}  
geneset.tot= names(GS)
if (length(con_template) != length(geneset.tot))   { geneset.tot= geneset.tot[-length(geneset.tot)] }
names(con_template)=geneset.tot
con_template.ori= con_template
if (length(gs.nullid)>0) {con_template= con_template[-gs.nullid]}
save(con_template, file='load_con_template_GO')
GS.ori= GS
if (length(gs.nullid)>0) {GS= GS[-gs.nullid]}
save(GS, file=paste('load_GS', annotation, sep='_'))

##############################################
## compute control GSR index = Rvalue_con  
load('load_con_template_CAN')
load(paste('load_GS', annotation, sep='_'))

casen= ncol(condata)-1
match.score.allcase.con= data.frame()

for ( k in 1: casen) {
  condata.k= condata[,k+1]
  
  matchscore.tot.onecase= vector()
  gsrank_onecase= function(i){       # i = ith geneset 
    
    gs.name= names(GS)[i]
    gs.gene= GS[[i]]
    gs.condata.id= which(is.element(condata$V1, gs.gene))
    gs.gene= gs.gene[is.element(gs.gene, condata$V1)]
    condata.k= condata.k[gs.condata.id]
    
    combination= combn(length(gs.gene), 2)
    combination.n= ncol(combination)
    
    rank.ings.onecase= vector()
    rank.ings= function(m) { 
      condata.gs.rank = ifelse(condata.k[combination[1,m]] >= condata.k[combination[2,m]], 'A','B')  
      rank.ings.onecase= c(rank.ings.onecase, condata.gs.rank)
      return(rank.ings.onecase)
    }
    rank.case= sapply(1:combination.n, rank.ings)
    
    matchone= ifelse(rank.case==con_template[[i]], 'y', 'n' )
    matchscore= length(matchone[matchone=='y'])/length(matchone)
    
    matchscore.tot.onecase= c(matchscore.tot.onecase, matchscore)
    return(matchscore.tot.onecase)
    
  }
  match.score= sapply(1:length(GS), gsrank_onecase)  
  cat('case number=  ', k, '\n')
  match.score.allcase.con= rbind(match.score.allcase.con, match.score)
  
}  
colnames(match.score.allcase.con)=names(GS)  
Rvalue.con_CAN= match.score.allcase.con
write.csv(match.score.allcase.con, paste('Rvalue_con', annotation, disease,sep='_'))
write.csv(names(GS),paste('Geneset', annotation, disease,sep='_') )

#####################################################################
## compute disease GSR index (Rvalue_dis) 

casen= ncol(disdata)-1
match.score.allcase.dis= data.frame()

for ( k in 1: casen) {
  disdata.k= disdata[,k+1]
  
  matchscore.tot.onecase= vector()
  gsrank_onecase= function(i){       # i = ith geneset 
    
    gs.name= names(GS)[i]
    gs.gene= GS[[i]]
    gs.disdata.id= which(is.element(disdata$V1, gs.gene))
    gs.gene= gs.gene[is.element(gs.gene, disdata$V1)]
    disdata.k= disdata.k[gs.disdata.id]
    
    combination= combn(length(gs.gene), 2)
    combination.n= ncol(combination)
    
    rank.ings.onecase= vector()
    rank.ings= function(m) { 
      disdata.gs.rank = ifelse(disdata.k[combination[1,m]] >= disdata.k[combination[2,m]], 'A','B')  
      rank.ings.onecase= c(rank.ings.onecase, disdata.gs.rank)
      return(rank.ings.onecase)
    }
    rank.case= sapply(1:combination.n, rank.ings)
    
    matchone= ifelse(rank.case==con_template[[i]], 'y', 'n' )
    matchscore= length(matchone[matchone=='y'])/length(matchone)
    
    matchscore.tot.onecase= c(matchscore.tot.onecase, matchscore)
    return(matchscore.tot.onecase)
    
  }
  match.score= sapply(1:length(GS), gsrank_onecase)  
  cat('case number=  ', k, '\n')
  match.score.allcase.dis= rbind(match.score.allcase.dis, match.score)
  
}  
colnames(match.score.allcase.dis)=names(GS)  
Rvalue.dis_CAN= match.score.allcase.dis
write.csv(match.score.allcase.dis, paste('Rvalue_dis', annotation, disease,sep='_'))

###############################################################################################################3
