#Code nécessaire pour répliquer les résultats de la conférence du 10 mai 2024
#Les données sont disponibles dans le répertoire 

set.seed(1234) #Pour s'assurer de la reproductibilité des résultats 

library(ggplot2)
library(ggpubr)

#Section valeur-p et inférence
Zscore = 1.4

distribution.nulle = data.frame("Zscore"=rnorm(10000, 0,1))

ggplot(distribution.nulle, aes(x=Zscore))+geom_histogram(bins = 100, color="darkblue", fill="lightblue")+geom_vline(xintercept = c(-1.4,1.4),linetype="dashed", linewidth=1.1)+theme_bw()

valeurp_Zscore = 2*pnorm(abs(Zscore), lower.tail = F)

plot(abs(distribution.nulle), 2*pnorm(abs(distribution.nulle), lower.tail = F))


#Section analyse différentielle
#Les données sont disponibles sur: https://www.kaggle.com/datasets/usharengaraju/indian-women-in-defense?select=airway_scaledcounts.csv
library(DESeq2)

count_genes = read.csv("airway_scaledcounts.csv", header=T)
rownames(count_genes) = count_genes$ensgene
count_genes$ensgene = NULL

design = data.frame("id" = colnames(count_genes), "condition" = rep(c("control", "treated"), 4), "celltype"= c("N61311","N61311",
                                                                                                               "N052611", "N052611", "N080611","N080611", "N061011","N061011"))
df.distrib = data.frame("counts"=count_genes$SRR1039521)
df.distrib = df.distrib[df.distrib$counts<=2000,,drop=F] #Filtre pour faciliter la visualisation
ggplot(df.distrib, aes(x=counts))+geom_histogram(bins=200,color="darkblue", fill="lightblue")+theme_bw()

dds = DESeqDataSetFromMatrix(countData = count_genes, colData = design, design = ~condition)
dds = DESeq(dds)

#Log-Log Plot de moyenne-variance 

df.counts = data.frame("baseMean"=log(mcols(dds)$baseMean+1), "baseVar"=log(mcols(dds)$baseVar+1))
ggplot(df.counts, aes(x=baseMean,y=baseVar))+geom_point(color="darkblue")+geom_abline(slope = 1, intercept = 0, color="red", size=1.1)+geom_smooth(method="lm")+theme_bw()+ylab("log(baseVar)")+xlab("log(baseMean)")

#Relation entre Moyennes des comptes et dispersion
#Valeurs par défaut
plotDispEsts(dds)

#D'autres choix sont cependant disponibles
#Choix d'une régression non-paramétrique locale
dds = estimateDispersionsFit(dds, fitType = "local")
plotDispEsts(dds)

#Choix d'une régression Gamma Poisson
dds = estimateDispersionsFit(dds, fitType = "glmGamPoi")
plotDispEsts(dds)

#Recalcul des MAP
dds = estimateDispersionsFit(dds, fitType = "local")
dds = estimateDispersionsMAP(dds)
plotDispEsts(dds)


#Résultats pénalisés 
#Intuitivement la pénalisation permet de différencier les variations dûent aux faibles comptes (larges variances)
#de ceux résultant d'information biologique pertinente
#Si un LFC est important et que la quantité d'information statistique est importante aussi (nombre de réplicats + counts modérés): pas de pénalisation

#C'est à dire que pour les gènes avec des valeurs de counts faibles on tend 
#à avoir des fold-change plus élevés

#Rappel: Le shrinkage permet de normaliser les gènes et donc de les rendre comparables.

r.def = results(dds)

df.pvalues = data.frame("pvalues"= r.def$pvalue)

gaston::qqplot.pvalues(df.pvalues$pvalues)#Nette inflation 

ggplot(df.pvalues, aes(x=pvalues))+geom_histogram(color="darkblue", fill="lightblue")+theme_bw()
hist(r.def$pvalue, main="Histogramme des valeurs-p", xlab="valeurs-p",ylab="Fréquences")#Nette inflation 

df.for.volcano = na.omit(data.frame("LFC"=r.def$log2FoldChange, "-log10(adjusted-p)"=-log10(r.def$padj)))
df.for.volcano$signi = df.for.volcano$X.log10.adjusted.p.>1

df.for.MA = na.omit(data.frame("LFC"=r.def$log2FoldChange, "Log2-Expression"=log2(r.def$baseMean), "pvalues"=r.def$padj))
df.for.MA$signi = df.for.MA$pvalues<=0.1

ggarrange(ggplot(df.for.MA, aes(y=LFC, x=Log2.Expression, color=signi))+geom_point()+theme_bw()+xlab("Log2 Expression"),
          ggplot(df.for.volcano, aes(x=LFC, y=X.log10.adjusted.p., color=signi))+geom_point()+theme_bw()+ylab("-log10(adjusted p-value)"))


r.def.pen = lfcShrink(dds, type="normal", coef=2)

df.for.volcano = na.omit(data.frame("LFC"=r.def.pen$log2FoldChange, "-log10(adjusted-p)"=-log10(r.def.pen$padj)))
df.for.volcano$signi = df.for.volcano$X.log10.adjusted.p.>1

df.for.MA = na.omit(data.frame("LFC"=r.def.pen$log2FoldChange, "Log2-Expression"=log2(r.def.pen$baseMean), "pvalues"=r.def.pen$padj))
df.for.MA$signi = df.for.MA$pvalues<=0.1

ggarrange(ggplot(df.for.MA, aes(y=LFC, x=Log2.Expression, color=signi))+geom_point()+theme_bw()+xlab("Log2 Expression"),
          ggplot(df.for.volcano, aes(x=LFC, y=X.log10.adjusted.p., color=signi))+geom_point()+theme_bw()+ylab("-log10(adjusted p-value)"))


