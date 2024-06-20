#Ce code permet de reproduire les résultats présentés durant la conférence de l'axe RSME du vendredi 17 mai 2024
#Les données sont disponibles: https://github.com/borenstein-lab/microbiome-metabolome-curated-data/tree/main/data/processed_data/FRANZOSA_IBD_2019

library(data.table) #Pour importer des gros jeux de données
library(mixOmics) 
library(multiview)
library(glmnet)
library(compositions) #Pour certaines fonctions de prétraitement
library(ggplot2)
library(FactoMineR)

species = fread("species.counts.tsv", header=TRUE, sep="\t")
metabolites = fread("mtb.tsv", header=TRUE, sep="\t")
metadata = read.table("metadata.tsv", header=TRUE, sep="\t")

species = data.frame(species, row.names = species$Sample)
species$Sample = NULL

clr.species = clr(species)

metabolites = data.frame(metabolites, row.names = metabolites$Sample)
metabolites$Sample = NULL


var_by_species = apply(species, 2, var)
top_500_most_var_species = order(var_by_species,decreasing = T)[1:500]
species_filter = species[,top_500_most_var_species] #On retire les espèces avec des faibles variances pour l'illustration

var_by_metabolites = apply(metabolites, 2, var)
top_500_most_var_metabolites = order(var_by_metabolites,decreasing = T)[1:500]
metabolites_filter = metabolites[,top_500_most_var_metabolites] #On retire les espèces avec des faibles variances pour l'illustration

phenotype = metadata[,c("Study.Group", "Sample"), drop=F]
phenotype$Status = ifelse(phenotype$Study.Group %in%c("UC","CD"), 1, 0)
y = phenotype$Status

table(y)#Distribution du phénotype

#1ère étape: Traitement des données
#Données de métagénomiques: compositionnelles --> transformation clr
#Données de métabolomique: log-normalisation

species_filter_clr = as.matrix(clr(species_filter))
metabolites_filter_log = as.matrix(log(metabolites_filter+1))

#Apprentissage non-supervisé
#PCA
pca_species = PCA(clr.species, scale.unit = F)
plot(pca_species, label="none", choix="ind")
plot(pca_species, label="none", choix="var")

pca_metabolites = PCA(metabolites, scale.unit = T)
plot(pca_metabolites, label="none", choix="ind")
plot(pca_species, label="none", choix="var")

#Hierarchical clustering 
d = dist(clr.species)
hclust_species = hclust(d)
plot(hclust_species, labels = F)

d1 = dist(metabolites)
hclust_metabolites = hclust(d1)
plot(hclust_metabolites, labels = F)

#Modèle mono-omique + modèle multivues
#On sépare le fichier de données en données d'entrainement et de test
#Ici 70%/30%
index_train = sample(1:nrow(species_filter_clr), nrow(species_filter_clr)*0.7)
index_test = (1:nrow(species_filter_clr))[-index_train]

species_filter_clr_train = species_filter_clr[index_train,]
metabolites_filter_log_train = metabolites_filter_log[index_train,]

species_filter_clr_test = species_filter_clr[index_test,]
metabolites_filter_log_test = metabolites_filter_log[index_test,]

#Prédiction sur le microbiome-metabolome: early fusion
cv_fit_species_metabolites = cv.glmnet(cbind(metabolites_filter_log_train,species_filter_clr_train),y[index_train],family = "binomial" , trace.it = 1)
plot(cv_fit_species_metabolites)

predict(cv_fit_species, newx=species_filter_clr_test, type="response", s="lambda.min")
predict(cv_fit_metabolites, newx=metabolites_filter_log_test, type="response", s="lambda.min")
predict(cv_fit_species_metabolites, newx=cbind(metabolites_filter_log_test, species_filter_clr_test), type="response", s="lambda.min")

#Courbes ROC et AUC
library(pROC)

roc(y[index_test], predict(cv_fit_metabolites, newx=metabolites_filter_log_test, type="response", s="lambda.min")[,1], percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    # print.thres = c(0.30,0.35, 0.40, 0.45,0.48, 0.50,0.55, 0.60),#
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE, ci.type="bars", print.thres.cex = 0.7, main="Roc Curve for Lasso with Metabolome as predictors")

roc(y[index_test], predict(cv_fit_species, newx=species_filter_clr_test, type="response", s="lambda.min")[,1], percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    # print.thres = c(0.30,0.35, 0.40, 0.45,0.48, 0.50,0.55, 0.60),#
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE, ci.type="bars", print.thres.cex = 0.7, main="Roc Curve for Lasso with Microbiome as predictors")

roc(y[index_test], predict(cv_fit_species_metabolites, newx=cbind(metabolites_filter_log_test,species_filter_clr_test), type="response", s="lambda.min")[,1], percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    # print.thres = c(0.30,0.35, 0.40, 0.45,0.48, 0.50,0.55, 0.60),#
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE, ci.type="bars", print.thres.cex = 0.7, main="Roc Curve for Lasso: Early-fusion")

roc(y[index_test], predict(model_late_fusion_train, newx=df_late_fusion_test, type="response"), percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    # print.thres = c(0.30,0.35, 0.40, 0.45,0.48, 0.50,0.55, 0.60),#
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE, ci.type="bars", print.thres.cex = 0.7, main="Roc Curve for Lasso: Late-fusion")

#Multiview --> Elastic-net: prédiction + sélection de features
cv_fit_train = cv.multiview(list(species_filter_clr_train,metabolites_filter_log_train),y[index_train], rho=0.9,family = binomial() , trace.it = TRUE)

plot(cv_fit_train)
coef(cv.fit, s="lambda.min")

predict(cv_fit_train, newx=list(species_filter_clr_test, metabolites_filter_log_test),type="response", s="lambda.min")

roc(y[index_test], predict(cv_fit_train, newx=list(species_filter_clr_test, metabolites_filter_log_test),type="response", s="lambda.min")[,1], percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    # print.thres = c(0.30,0.35, 0.40, 0.45,0.48, 0.50,0.55, 0.60),#
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE, ci.type="bars", print.thres.cex = 0.7, main="Roc Curve for Lasso with Microbiome and Metabolome as predictors")

library(randomForest)

rf_train = randomForest(cbind(metabolites_filter_log_train,species_filter_clr_train),as.factor(y[index_train]))

roc(y[index_test], predict(cv_fit_species_metabolites, newx=cbind(metabolites_filter_log_test,species_filter_clr_test), type="response", s="lambda.min")[,1], percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    # print.thres = c(0.30,0.35, 0.40, 0.45,0.48, 0.50,0.55, 0.60),#
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE, ci.type="bars", print.thres.cex = 0.7, main="Roc Curve for Elastic-Net: Early-fusion")

roc(y[index_test], predict(cv_fit_train, newx=list(species_filter_clr_test, metabolites_filter_log_test),type="response", s="lambda.min")[,1], percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    # print.thres = c(0.30,0.35, 0.40, 0.45,0.48, 0.50,0.55, 0.60),#
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE, ci.type="bars", print.thres.cex = 0.7, main="Roc Curve for Elastic-Net: Multi-View")

roc(y[index_test], predict(rf_train, cbind(metabolites_filter_log_test,species_filter_clr_test), type = "prob")[,2], percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    # print.thres = c(0.30,0.35, 0.40, 0.45,0.48, 0.50,0.55, 0.60),#
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE, ci.type="bars", print.thres.cex = 0.7, main="Roc Curve for Random-Forest: Early-Fusion")


#sPLS-DA
X_train = list("metabo"= metabolites_filter_log_train, "micro"=species_filter_clr_train)
X_test = list("metabo"= metabolites_filter_log_test, "micro"=species_filter_clr_test)

list.keepX = list(metabo = c(50, 25), micro=c(50,25)) 
result_sparse_diablo_train <-  block.splsda(X=X_train, y[index_train], keepX = list.keepX) 

predict_spls_DA = predict(result_sparse_diablo_train, X_test, type="response")
p = predict_spls_DA$AveragedPredict[,2,1]
p[p>1] = 1

roc(y[index_test], p, percent=F,   boot.n=1000, ci.alpha=0.9, stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    # print.thres = c(0.30,0.35, 0.40, 0.45,0.48, 0.50,0.55, 0.60),#
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE, ci.type="bars", print.thres.cex = 0.7, main="Roc Curve for block sPLS-DA")

#Calibration
library(rms)
val.prob(predict(rf_train, cbind(metabolites_filter_log_test,species_filter_clr_test), type = "prob")[,2], y[index_test])
val.prob(predict(cv_fit_train, newx=list(species_filter_clr_test, metabolites_filter_log_test),type="response", s="lambda.min")[,1], y[index_test])
val.prob(predict(cv_fit_species_metabolites, newx=cbind(metabolites_filter_log_test,species_filter_clr_test),type="response", s="lambda.min")[,1], y[index_test])
val.prob(p, y[index_test])

#Cette partie n'est pas totalement fonctionnelle et ne sert que de structure pour une late fusion
#Prédiction sur le microbiome: late fusion
cv_fit_species = cv.glmnet(species_filter_clr_train,y[index_train],family = "binomial" , trace.it = 1)
plot(cv_fit_species)

#Prédiction sur le metabolome: late fusion
cv_fit_metabolites = cv.glmnet(metabolites_filter_log_train,y[index_train],family = "binomial" , trace.it = 1)
plot(cv_fit_metabolites)


c0 = colnames(metabolites_filter_log[,which(matrix(coef(cv_fit_metabolites, s="lambda.min"))!=0)])
c1 = colnames(species_filter_clr[,which(matrix(coef(cv_fit_species, s="lambda.min"))!=0)])

df_late_fusion_test = cbind(metabolites_filter_log_test, species_filter_clr_test)[,c(c0,c1)]
df_late_fusion_train = cbind(metabolites_filter_log_train, species_filter_clr_train)[,c(c0,c1)]

model_late_fusion_train = glm(y[index_train]~df_late_fusion_train, family=binomial())
