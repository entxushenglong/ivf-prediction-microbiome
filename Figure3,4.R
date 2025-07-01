library(tidyverse)
library(reshape2)
library(glmnet)
library(pROC)
library(dplyr)
library(caret)
library(shapviz)
library(xgboost)
library(ggplot2)
library(fbroc)
library(boot)
library(ggpubr)
library(plot3D)
library(factoextra)
library(FactoMineR)
library(vegan)
library(ggalt)
library(scatterplot3d) 
library(grid)
library(gg.gap)
library(cvAUC)
library(forestplot)


data_new<-read.csv("otu_table_filtered.csv",fill = T)
rownames(data_new)<-as.character(data_new$Genus)
data_new<-data_new[,-1]
sample_class<-c()
for (i in 1:nrow(data_new)) {
  t<-nchar(rownames(data_new)[i])
  sample_class[i]<-substr(rownames(data_new)[i],(t-1),t)
}
clin<-read.csv("clinical_add.csv")


#Figure3A/Figure4A
yd_data<-data_new[which(sample_class=="YD"),]
yd_check<-colSums(yd_data)
yd_data<-yd_data[,which(yd_check != 0)] 

yd<-read.csv("yd.csv",header = F)
rownames(yd)<-yd$V1
yd<-yd[rownames(yd_data),]

for (i in 1:nrow(yd)) {
  yd[i,"sample"]<-substr(yd$V1[i],1,(nchar(yd$V1[i])-2))
}
yd_clin<-clin[which(clin$SampleID %in% yd$sample),]
for (i in 1:nrow(yd_clin)) {
  rownames(yd_clin)[i]<-paste0(yd_clin$SampleID[i],"YD",collapse = "") 
}
cha3_ying<-c()
for (i in 1:ncol(yd_data)) {
  if(length(which(yd_data[,i]>0))>2){
    cha3_ying<-c(cha3_ying,i)
  }
}

x_yd<-yd_data[rownames(yd_clin),cha3_ying]
m3_yd <- glmnet(as.matrix(x_yd),as.matrix(yd_clin$result), family="binomial",standardize=TRUE)
best_lambda=min(m3_yd[["lambda"]])
yd_train_pre<-as.numeric(predict(m3_yd, s = best_lambda, newx = as.matrix(x_yd)))
auc(roc(yd_clin$result,yd_train_pre))
plot(roc(yd_clin$result,yd_train_pre),print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), 
     grid.col=c("green","red"), max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T,
     main="")

GJ_data<-data_new[which(sample_class=="GJ"),]
GJ_check<-colSums(GJ_data)
GJ_data<-GJ_data[,which(GJ_check != 0)]
GJ<-read.csv("gj.csv",header = F)
rownames(GJ)<-GJ$V1
GJ<-GJ[rownames(GJ_data),]

cha3_gong<-c()
for (i in 1:ncol(GJ_data)) {
  if(length(which(GJ_data[,i]>0))>2){
    cha3_gong<-c(cha3_gong,i)
  }
}
for (i in 1:nrow(GJ)) {
  GJ[i,"sample"]<-substr(GJ$V1[i],1,(nchar(GJ$V1[i])-2))
}
GJ_clin<-clin[which(clin$SampleID %in% GJ$sample),]
for (i in 1:nrow(GJ_clin)) {
  rownames(GJ_clin)[i]<-paste0(GJ_clin$SampleID[i],"GJ",collapse = "") 
}
x_GJ<-GJ_data[rownames(GJ_clin),cha3_gong]
m3_GJ <- glmnet(as.matrix(x_GJ),as.matrix(GJ_clin$result), family="binomial",standardize=TRUE)
best_lambda=min(m3_GJ [["lambda"]])
GJ_train_pre<-as.numeric(predict(m3_GJ, s = best_lambda, newx = as.matrix(x_GJ)))
auc(roc(GJ_clin$result,GJ_train_pre))
plot(roc(GJ_clin$result,GJ_train_pre),print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), 
     grid.col=c("green","red"), max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T,
     main="")

GQ_data<-data_new[which(sample_class=="ZG"),]
GQ_check<-colSums(GQ_data)
GQ_data<-GQ_data[,which(GQ_check != 0)]
GQ<-read.csv("GQ.csv",header = F)
GQ_set<-intersect(GQ$V1,rownames(GQ_data))
GQ_data<-GQ_data[GQ_set,]
GQ<-GQ[which(GQ$V1 %in% GQ_set),]
rownames(GQ)<-GQ$V1
GQ<-GQ[GQ_set,]
cha3_qiang<-c()
for (i in 1:nrow(GQ_data)) {
  if(length(which(GQ_data[,i]>0))>2){
    cha3_qiang<-c(cha3_qiang,i)
  }
}
for (i in 1:nrow(GQ)) {
  GQ[i,"sample"]<-substr(GQ$V1[i],1,(nchar(GQ$V1[i])-3))
}
GQ_clin<-clin[which(clin$SampleID %in% GQ$sample),]
for (i in 1:nrow(GQ_clin)) {
  rownames(GQ_clin)[i]<-paste0(GQ_clin$SampleID[i],"YZG",collapse = "") 
}
x_GQ<-GQ_data[rownames(GQ_clin),cha3_qiang]
m3_GQ <- glmnet(as.matrix(x_GQ),as.matrix(GQ_clin$result), family="binomial",standardize=TRUE)
best_lambda <-min(m3_GJ [["lambda"]])
GQ_train_pre<-as.numeric(predict(m3_GQ, s = best_lambda, newx = as.matrix(x_GQ)))
auc(roc(GQ_clin$result,GQ_train_pre))
plot(roc(GQ_clin$result,GQ_train_pre),print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), 
     grid.col=c("green","red"), max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T,
     main="")

#Figure3B
CV=function(n,Z=10,seed=888){
  z=rep(1:Z,ceiling(n/Z))[1:n]
  set.seed(seed)
  z=sample(z,n)
  mm=list()
  #mm[[i]]为第i个下标集
  for (i in 1:Z) mm[[i]]=(1:n)[z==i];return(mm)
}
n=nrow(x);Z=10;mm=CV(n,Z);D=1
yd_fold<-data.frame(n=1:10,train=0,test=0)
for(i in 1:Z){   #循环十次
  m=mm[[i]];
  x1=as.matrix(t(yd_data[cha3_ying,-m]))
  y1=as.matrix(yd[1,-m])
  m1<-glmnet(x1, y1, family="binomial",standardize=TRUE)
  best_lambda <- min(m1[["lambda"]])
  m1_train_pre<-as.numeric(predict(m1, s = best_lambda, newx = x1))
  yd_fold[i,"train"]<-auc(roc(as.matrix(yd[1,-m]),m1_train_pre))
  m1_test_pre<-as.numeric(predict(m1, s = best_lambda, newx = as.matrix(x[m,])))
  yd_fold[i,"test"]<-auc(roc(as.matrix(yd[1,m]),m1_test_pre))
}



yd_roc_5folds<-list()
yd_auc_5folds<-c()
ydclin_roc_5folds<-list()
ydclin_auc_5folds<-c()
x_yd_clinadd<-cbind(yd_train_pre,yd_clin[,2:12])
ydadd_roc_5folds<-list()
ydadd_auc_5folds<-c()
i=1
while (i < 1:100) {
  for(i in 1:Z){
    m=mm[[i]];
    #micro
    x1=as.matrix(x_yd[-m,])
    y1=as.matrix(yd_clin$result[-m])
    m1<-glmnet(x1, y1, family="binomial",standardize=TRUE)
    best_lambda <- min(m1[["lambda"]])
    yd_5folds_pre<-as.numeric(predict(m1, s = best_lambda, newx = as.matrix(x_yd[m,])))
    yd_roc_5folds[[i]]<-roc(yd_clin$result,yd_5folds_pre)[["predictor"]]
    yd_auc_5folds[i]<-roc(yd_clin$result,yd_5folds_pre)[["auc"]]
    #clin
    x2=as.matrix(yd_clin[-m,3:12])
    m2<-glmnet(x2, y1, family="binomial",standardize=TRUE)
    best_lambda <- min(m2[["lambda"]])
    ydclin_5folds_pre<-as.numeric(predict(m2, s = best_lambda, newx = as.matrix(yd_clin[m,3:12])))
    ydclin_roc_5folds[[i]]<-roc(yd_clin$result,ydclin_5folds_pre)[["predictor"]]
    ydclin_auc_5folds[i]<-roc(yd_clin$result,ydclin_5folds_pre)[["auc"]]
    #add
    x3=as.matrix(x_yd_clinadd[-m,-2])
    m3<-glmnet(x3, y1, family="binomial",standardize=TRUE)
    best_lambda <- min(m3[["lambda"]])
    ydadd_5folds_pre<-as.numeric(predict(m3, s = best_lambda, newx = as.matrix(x_yd_clinadd[m,-2])))
    ydadd_roc_5folds[[i]]<-roc(yd_clin$result,ydadd_5folds_pre)[["predictor"]]
    ydadd_auc_5folds[i]<-roc(yd_clin$result,ydadd_5folds_pre)[["auc"]]   
  }
}
fold_yd<-min(yd_auc_5folds[which(yd_auc_5folds>=median(yd_auc_5folds))])
rocyd <- plot.roc(yd_clin$result,yd_roc_5folds[[which(yd_auc_5folds==fold_yd)[1]]],percent=F,ci=TRUE,print.auc=F)
ciyd <- ci.se(rocyd,specificities=seq(0, 1, 0.2))
plot(ci(rocyd, of="thresholds", thresholds="best"))
plot(ciyd, type="shape", col=rgb(248,118,109,60, maxColorValue=255))
fold_ydclin<-min(ydclin_auc_5folds[which(ydclin_auc_5folds>=median(ydclin_auc_5folds))])
fold_ydadd<-min(ydadd_auc_5folds[which(ydadd_auc_5folds>=median(ydadd_auc_5folds))])

clin_yd<-ydclin_roc_5folds[[which(ydclin_auc_5folds==fold_ydclin)[1]]]
micro_yd<-yd_roc_5folds[[which(yd_auc_5folds==fold_yd)[1]]]
miclin_yd<-ydadd_roc_5folds[[which(ydadd_auc_5folds==fold_ydadd)[1]]]
plot(roc(yd_clin$result, clin_yd))
plot(roc(yd_clin$result, micro_yd), add=TRUE, col="#F8766D",lty=2)
plot(roc(yd_clin$result, miclin_yd), add=TRUE, col="#F8766D",lwd=3)
auc(roc(yd_clin$result, miclin_yd))





GJ_roc_5folds<-list()
GJ_auc_5folds<-c()
GJclin_roc_5folds<-list()
GJclin_auc_5folds<-c()
x_GJ_clinadd<-cbind(GJ_train_pre,GJ_clin[,2:12])
GJadd_roc_5folds<-list()
GJadd_auc_5folds<-c()
i=1
while (i < 1:100) {
  for(i in 1:Z){
    m=mm[[i]];
    #micro
    x1=as.matrix(x_GJ[-m,])
    y1=as.matrix(GJ_clin$result[-m])
    m1<-glmnet(x1, y1, family="binomial",standardize=TRUE)
    best_lambda <- min(m1[["lambda"]])
    GJ_5folds_pre<-as.numeric(predict(m1, s = best_lambda, newx = as.matrix(x_GJ[m,])))
    GJ_roc_5folds[[i]]<-roc(GJ_clin$result,GJ_5folds_pre)[["predictor"]]
    GJ_auc_5folds[i]<-roc(GJ_clin$result,GJ_5folds_pre)[["auc"]]
    #clin
    x2=as.matrix(GJ_clin[-m,3:12])
    m2<-glmnet(x2, y1, family="binomial",standardize=TRUE)
    best_lambda <- min(m2[["lambda"]])
    GJclin_5folds_pre<-as.numeric(predict(m2, s = best_lambda, newx = as.matrix(GJ_clin[m,3:12])))
    GJclin_roc_5folds[[i]]<-roc(GJ_clin$result,GJclin_5folds_pre)[["predictor"]]
    GJclin_auc_5folds[i]<-roc(GJ_clin$result,GJclin_5folds_pre)[["auc"]]
    #add
    x3=as.matrix(x_GJ_clinadd[-m,-2])
    m3<-glmnet(x3, y1, family="binomial",standardize=TRUE)
    best_lambda <- min(m3[["lambda"]])
    GJadd_5folds_pre<-as.numeric(predict(m3, s = best_lambda, newx = as.matrix(x_GJ_clinadd[m,-2])))
    GJadd_roc_5folds[[i]]<-roc(GJ_clin$result,GJadd_5folds_pre)[["predictor"]]
    GJadd_auc_5folds[i]<-roc(GJ_clin$result,GJadd_5folds_pre)[["auc"]]   
  }
}
fold_GJ<-min(GJ_auc_5folds[which(GJ_auc_5folds>=median(GJ_auc_5folds))])
rocGJ <- plot.roc(GJ_clin$result,GJ_roc_5folds[[which(GJ_auc_5folds==fold_GJ)[1]]],percent=F,ci=TRUE,print.auc=F)
ciGJ <- ci.se(rocGJ,specificities=seq(0, 1, 0.2))
plot(ci(rocGJ, of="thresholds", thresholds="best"))
plot(ciGJ, type="shape", col=rgb(0, 186, 56,60, maxColorValue=255))
fold_GJclin<-min(GJclin_auc_5folds[which(GJclin_auc_5folds>=median(GJclin_auc_5folds))])
fold_GJadd<-min(GJadd_auc_5folds[which(GJadd_auc_5folds>=median(GJadd_auc_5folds))])

clin_GJ<-GJclin_roc_5folds[[which(GJclin_auc_5folds==fold_GJclin)[1]]]
micro_GJ<-GJ_roc_5folds[[which(GJ_auc_5folds==fold_GJ)]]
miclin_GJ<-GJadd_roc_5folds[[which(GJadd_auc_5folds==fold_GJadd)[1]]]
plot(roc(GJ_clin$result, clin_GJ))
plot(roc(GJ_clin$result, micro_GJ), add=TRUE, col="#00BA38",lty=2)
plot(roc(GJ_clin$result, miclin_GJ), add=TRUE, col="#00BA38",lwd=3)



GQ_roc_5folds<-list()
GQ_auc_5folds<-c()
GQclin_roc_5folds<-list()
GQclin_auc_5folds<-c()
x_GQ_clinadd<-cbind(GQ_train_pre,GQ_clin[,2:12])
GQadd_roc_5folds<-list()
GQadd_auc_5folds<-c()
i=1
while (i < 1:100) {
  for(i in 1:Z){
    m=mm[[i]];
    #micro
    x1=as.matrix(x_GQ[-m,])
    y1=as.matrix(GQ_clin$result[-m])
    m1<-glmnet(x1, y1, family="binomial",standardize=TRUE)
    best_lambda <- min(m1[["lambda"]])
    GQ_5folds_pre<-as.numeric(predict(m1, s = best_lambda, newx = as.matrix(x_GQ[m,])))
    GQ_roc_5folds[[i]]<-roc(GQ_clin$result,GQ_5folds_pre)[["predictor"]]
    GQ_auc_5folds[i]<-roc(GQ_clin$result,GQ_5folds_pre)[["auc"]]
    #clin
    x2=as.matrix(GQ_clin[-m,3:12])
    m2<-glmnet(x2, y1, family="binomial",standardize=TRUE)
    best_lambda <- min(m2[["lambda"]])
    GQclin_5folds_pre<-as.numeric(predict(m2, s = best_lambda, newx = as.matrix(GQ_clin[m,3:12])))
    GQclin_roc_5folds[[i]]<-roc(GQ_clin$result,GQclin_5folds_pre)[["predictor"]]
    GQclin_auc_5folds[i]<-roc(GQ_clin$result,GQclin_5folds_pre)[["auc"]]
    #add
    x3=as.matrix(x_GQ_clinadd[-m,-2])
    m3<-glmnet(x3, y1, family="binomial",standardize=TRUE)
    best_lambda <- min(m3[["lambda"]])
    GQadd_5folds_pre<-as.numeric(predict(m3, s = best_lambda, newx = as.matrix(x_GQ_clinadd[m,-2])))
    GQadd_roc_5folds[[i]]<-roc(GQ_clin$result,GQadd_5folds_pre)[["predictor"]]
    GQadd_auc_5folds[i]<-roc(GQ_clin$result,GQadd_5folds_pre)[["auc"]]   
  }
}
fold_GQ<-min(GQ_auc_5folds[which(GQ_auc_5folds>=median(GQ_auc_5folds))])
rocGQ <- plot.roc(GQ_clin$result,GQ_roc_5folds[[which(GQ_auc_5folds==fold_GQ)[1]]],percent=F,ci=TRUE,print.auc=F)
ciGQ <- ci.se(rocGQ,specificities=seq(0, 1, 0.2)) # over a select set of specificities
plot(ci(rocGQ, of="thresholds", thresholds="best")) # add one threshold
plot(ciGQ, type="shape", col=rgb(97,156,255,60, maxColorValue=255))
fold_GQclin<-min(GQclin_auc_5folds[which(GQclin_auc_5folds>=median(GQclin_auc_5folds))])
fold_GQadd<-min(GQadd_auc_5folds[which(GQadd_auc_5folds>=median(GQadd_auc_5folds))])

clin_GQ<-GQclin_roc_5folds[[which(GQclin_auc_5folds==fold_GQclin)[1]]]
micro_GQ<-GQ_roc_5folds[[which(GQ_auc_5folds==fold_GQ)[1]]]
miclin_GQ<-GQadd_roc_5folds[[which(GQadd_auc_5folds==fold_GQadd)[1]]]
plot(roc(GQ_clin$result, clin_GQ))
plot(roc(GQ_clin$result, micro_GQ), add=TRUE, col="#619CFF",lty=2)
plot(roc(GQ_clin$result, miclin_GQ), add=TRUE, col="#619CFF",lwd=3)
auc(roc(GQ_clin$result, miclin_GQ))

#Figure3C-D
bst_yd <- xgboost(data = as.matrix(x_yd), label = yd_clin$result,
                  max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
                  objective = "binary:logistic")
yd_XG<- predict(bst_yd, as.matrix(x_yd))
roc_XG <- roc(yd_clin$result, yd_XG)
plot(roc_XG, main = "ROC Curve", print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, legacy.axes = TRUE,col="blue")
X_shap <- data.matrix(x_yd[sample(nrow(x_yd), 90), colnames(x_yd)])
shp <- shapviz(bst_yd, X_pred = X_shap)
sv_force(shp,row_id = 2)  # force plot
sv_importance(shp,kind = "beeswarm")
sv_importance(shp)

bst_GJ <- xgboost(data = as.matrix(x_GJ), label = GJ_clin$result,
                  max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
                  objective = "binary:logistic")
GJ_XG<- predict(bst_GJ, as.matrix(x_GJ))
roc_XG <- roc(GJ_clin$result, GJ_XG)
plot(roc_XG, main = "ROC Curve", print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, legacy.axes = TRUE,col="blue")

#SHAP
X_shap <- data.matrix(x_GJ[sample(nrow(x_GJ), 90), colnames(x_GJ)])
shp <- shapviz(bst_GJ, X_pred = X_shap)
sv_force(shp,row_id = 2)  # force plot
sv_importance(shp,kind = "beeswarm")
sv_importance(shp)

bst_GQ <- xgboost(data = as.matrix(x_GQ), label = GQ_clin$result,
                  max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
                  objective = "binary:logistic")
GQ_XG<- predict(bst_GQ, as.matrix(x_GQ))
roc_XG <- roc(GQ_clin$result, GQ_XG)
plot(roc_XG, main = "ROC Curve", print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, legacy.axes = TRUE,col="blue")

#SHAP
X_shap <- data.matrix(x_GQ[sample(nrow(x_GQ), 90,replace = T), colnames(x_GQ)])
shp <- shapviz(bst_GQ, X_pred = X_shap)3
sv_force(shp,row_id = 2)  # force plot
sv_importance(shp,kind = "beeswarm")
sv_importance(shp)

#Figure5B/C
bst_clin <- xgboost(data = as.matrix(clin[,3:12]), label = clin$result,
                    max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
                    objective = "binary:logistic")
clin_XG<- predict(bst_clin, as.matrix(clin[,3:12]))
roc_XG <- roc(clin$result, clin_XG)
plot(roc_XG, main = "ROC Curve", print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, legacy.axes = TRUE,col="blue")

#SHAP
X_shap <- data.matrix(clin[sample(nrow(clin), 90,replace = T), 3:12])
shp <- shapviz(bst_clin, X_pred = X_shap)
sv_force(shp,row_id = 2)  # force plot
sv_importance(shp,kind = "beeswarm")
sv_importance(shp)
