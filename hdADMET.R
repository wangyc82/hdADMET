# this procedure is to perform hdADMET
# the training data came from GDSCandCTRP combination data
# prepare the input data
## Acd (adjacent matrix representing the interactions between new compounds and known clinical anticancer drugs)
## Adp (adjacent matrix representing the interactions between drug and their ADMET properties)
## App (adjacent matrix representing ADMET inter-relationships)

# take example data for instance
# example data includes the MACCs fingerprints for new compounds (MACCs.fp.mat) and anticancer drugs (MACCs.fp.mat.CA)
# adjacent matrix representing drug-ADMET interaction,adj_ADMET_CA1
# adjacent matrix representing ADMET inter-relationships, ADMET_cor that were obtained by the ADMET profiles of drugs


# preparing Acd
# sim.MACCs.mat.new<-matrix(0,nrow(MACCs.fp.mat),nrow(MACCs.fp.mat.CA))
# rownames(sim.MACCs.mat.new)<-rownames(MACCs.fp.mat)
# colnames(sim.MACCs.mat.new)<-rownames(MACCs.fp.mat.CA)
# for (i in 1:nrow(MACCs.fp.mat)) {for (j in 1:nrow(MACCs.fp.mat.CA)) {sim.MACCs.mat.new[i,j]<-fpSim(MACCs.fp.mat[i,], MACCs.fp.mat.CA[j,], method="Tanimoto")}}
# Acd.full<-matrix(0,nrow(sim.MACCs.mat.new),ncol(sim.MACCs.mat.new))
# rownames(Acd.full)<-rownames(sim.MACCs.mat.new)
# colnames(Acd.full)<-colnames(sim.MACCs.mat.new)
# Acd.full[which(sim.MACCs.mat.new>=0.3)]<-1, threshold 0.3 determined by the distribution of as.vector(sim.MACCs.mat.new) and make sure length(which(rowSums(Acd.full)==0))==0
# check the scale of Acd.full, to make sure Acd do not contain zero column adn row
#length(which(rowSums(Acd.full)==0))
#[1] 0
#length(which(colSums(Acd.full)==0))
#[1] 69
#Acd<-Acd.full[,which(colSums(Acd.full)!=0)]

# prepare Adp
# Adp.full<-adj_ADMET_CA1[colnames(Acd),]
#length(which(rowSums(Adp.full)==0))
#[1] 0
# length(which(colSums(Adp.full)==0))
#[1] 0
# Adp<-Adp.full

# prepare App
# App.full<-matrix(0,ncol(adj_ADMET_CA1),ncol(adj_ADMET_CA1))
# rownames(App.full)<-colnames(adj_ADMET_CA1)
# colnames(App.full)<-colnames(adj_ADMET_CA1)
# App.full[which(abs(ADMET_cor)>0.4)]<-1
#length(which(colSums(App.full)==0))
#[1] 0
# length(which(rowSums(App.full)==0))
# [1] 0
# App<-App.full

# make sure compound and ADMET name has no any delimeter

# set len1, len2, len3, len4, len5 accoding the dimension of Acd, Adp, App
# for this example case, we set len1=10000;len2=1e5;len3=1e5;len4=1000;len5=5*1e+5;
#library(ChemmineOB)
#library(ChemmineR)
hdADMET<-function(Acd, Adp, App,len1,len2,len3,len4,len5) {
  Cd<-lapply(1:nrow(Acd), function(x) colnames(Acd)[which(Acd[x,]==1)])
  names(Cd)<-rownames(Acd)
  Dc<-lapply(1:ncol(Acd), function(x) rownames(Acd)[which(Acd[,x]==1)])
  names(Dc)<-colnames(Acd)
  Dp<-lapply(1:nrow(Adp), function(x) colnames(Adp)[which(Adp[x,]==1)])
  names(Dp)<-rownames(Adp)
  P<-lapply(1:ncol(Adp), function(x) rownames(Adp)[which(Adp[,x]==1)])
  names(P)<-colnames(Adp)
  P1<-lapply(1:nrow(App), function(x) colnames(App)[which(App[x,]==1)])
  names(P1)<-rownames(App)
  P2<-lapply(1:ncol(App), function(x) rownames(App)[which(App[,x]==1)])
  names(P2)<-colnames(App)
  
  #then generate the random walk based on each meta-path
  #for C-D-C meth-path
  #set len1 is a pre-defined parameter that to determine how many random walks will be generated based on this meta-path
  I=rownames(Acd)
  L1<-character()
  for (i in 1:len1) {
    p<-sample(I,1)
    m<-Cd[[p]]
    n<-sample(m,1)
    s<-Dc[[n]]
    L1[i]<-paste(p,n,sample(s,1),sep = " ")
    rm(p,m,n,s)
  }
  
  #for D-P-D meth-path
  #set len2 is a pre-defined parameter that to determine how many random walks will be generated based on this meta-path
  I=rownames(Adp)
  L2<-character()
  for (i in 1:len2) {
    p<-sample(I,1)
    m<-Dp[[p]]
    n<-sample(m,1)
    s<-P[[n]]
    L2[i]<-paste(p,n,sample(s,1),sep = " ")
    rm(p,m,n,s)
  }
  
  #for C-D-P-D-C meth-path
  #set len3 is a pre-defined parameter that to determine how many random walks will be generated based on this meta-path
  I=rownames(Acd)
  L3<-character()
  for (i in 1:len3) {
    p<-sample(I,1)
    m<-Cd[[p]]
    n<-sample(m,1)
    s<-Dp[[n]]
    o<-sample(s,1)
    e<-P[[o]]
    c<-sample(e,1)
    d<-Dc[[c]]
    L3[i]<-paste(p,n,o,c,sample(d,1),sep = " ")
    rm(p,m,n,s)
  }
  
  #if adding PPI network, for P-P-P meth-path
  #set len4 is a pre-defined parameter that to determine how many random walks will be generated based on this meta-path
  I=rownames(App)
  L4<-character()
  for (i in 1:len4) {
    p<-sample(I,1)
    m<-P1[[p]]
    n<-sample(m,1)
    s<-P2[[n]]
    L4[i]<-paste(p,n,sample(s,1),sep = " ")
    rm(p,m,n,s)
  }
  
  #if adding PPI network, for C-D-P-P-D-C meth-path
  #set len5 is a pre-defined parameter that to determine how many random walks will be generated based on this meta-path
  I=rownames(Acd)
  L5<-character()
  for (i in 1:len5) {
    p<-sample(I,1)
    m<-Cd[[p]]
    n<-sample(m,1)
    s<-Dp[[n]]
    o<-sample(s,1)
    w<-P1[[o]]
    e<-sample(w,1)
    f<-P2[[e]]
    a<-sample(f,1)
    b<-P[[a]]
    c<-sample(b,1)
    d<-Dc[[c]]
    L5[i]<-paste(p,n,o,e,a,c,sample(d,1),sep = " ")
  }
  library("word2vec")
  model <- word2vec(x = c(L1,L2,L3,L4,L5), type = "skip-gram", dim = 50, iter = 20)
  embedding <- as.matrix(model)
  labC<-match(rownames(Acd),rownames(embedding))
  labD<-match(rownames(Adp),rownames(embedding))
  labP<-match(colnames(Adp),rownames(embedding))
  Xcomp<-embedding[labC,]
  Xdrug<-embedding[labD,]
  Xpro<-embedding[labP,]
  Xcomppound<-rbind(Xdrug,Xcomp)
  simXcompound<-exp(-0.1*as.matrix(dist(Xcomppound)))
  simXcomp<-simXcompound[rownames(Xcomp),rownames(Xdrug)]
  simXdrug<-simXcompound[rownames(Xdrug),rownames(Xdrug)]
  simXpro<-exp(-0.1*as.matrix(dist(Xpro)))
  A<-which(Adp==1,arr.ind = TRUE)
  B<-which(Adp==0,arr.ind = TRUE)
  Xpos<-cbind(simXdrug[A[,1],],simXpro[A[,2],])
  Xneg<-cbind(simXdrug[B[,1],],simXpro[B[,2],])
  Xrn<-rbind(Xpos,Xneg)
  Yrn<-rep(c(1,0),c(nrow(A),nrow(B)))
  Ast<-matrix(1,nrow(Acd),ncol(App))
  Cst<-which(Ast==1,arr.ind = TRUE)
  Xst<-cbind(simXcomp[Cst[,1],],simXpro[Cst[,2],])
  Yst<-rep(1,nrow(Xst))
  library(h2o)
  h2o.init()
  data_trn=cbind(Xrn,Yrn)
  data_tst=cbind(Xst,Yst)
  colnames(data_tst)<-colnames(data_trn)
  train.hex <- as.h2o(data_trn)
  test.hex <- as.h2o(data_tst)
  model.DNN=h2o.deeplearning(x = 1:ncol(Xrn), y = ncol(Xrn)+1, training_frame = train.hex, hidden=c(200,200), epochs=10, activation="Tanh",seed = 10,reproducible = TRUE)
  prob_test_DNN<-as.data.frame(h2o.predict(model.DNN, test.hex))[,1]
  prdnewdrug<-matrix(prob_test_DNN,ncol=ncol(Ast),byrow=FALSE)
  rownames(prdnewdrug)<-rownames(Acd)
  colnames(prdnewdrug)<-colnames(App)
  h2o.shutdown()
  return(prdnewdrug)
}
