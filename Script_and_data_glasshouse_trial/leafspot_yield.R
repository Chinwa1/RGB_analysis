###get the working directory
getwd()
###set the working directory
###load the data, can use the import function
rust=read.csv("GRV.csv",header = TRUE)
###attach the data
is.na(rust)
attach(rust)
###check the structure of the data
str(sample)
###check the first 6 rows
head(sample)
tail(sample)
###check the class of the data
class(sample)

###########open libraries
library(dplyr)
library(agricolae)
library(lme4)
library(lmerTest)
###make the columns factors
rust$Genotypes=as.factor(rust$Genotypes)
sample$Genotype=as.factor(sample$Genotype)
sample$Rep=as.factor(sample$Rep) 
sample$Block=as.factor(sample$Block)
######check the structure of the data
str(sample)
###check for the normality of data
shapiro.test(rust$Pod_yield)
###plot a normality line
library(tidyselect)
library(ggplot2)
qqnorm(sample$Severity)
qqline(sample$Severity)
qqplot(sample$Severity,rust$HTP_Severity)
# Using to explore the distribution of a variable
ggplot(sample, aes(sample = Severity)) +
  stat_qq() +
  stat_qq_line()
ggplot(sample, aes(sample = Severity, colour = factor(Genotype))) +
  stat_qq() +
  stat_qq_line()+
  theme_bw()+
  theme(element_blank())
###compute the model lmer
library(lme4)
library(lmerTest)
library(agricolae)

rustnova=lm(Pod_yield~Genotypes+Replication,rust)
rustnova=lm(maturity$b.~genotype,maturity)
attributes(rustnova)
rustnova$rank
coef(rustnova)
summary(rustnova)
anova(rustnova)

Lightness=anova(rustnova)
write.csv(Severity,"ANOVA FOr rust.csv")
###save the output as word###
sink(" Rust ANOVA.doc") #create/open word file named "output" in working directory
print(Incubation_period)
print(u.)
getwd()
##############################USE lmer with random factors
rustnovalmer=(lmer(Severity~Stage*genotype+(1|Rep),rust))
###compute the ANOVA
rustnovalmer1=anova(rustnovalmer)
## Anova-like table of random-effect terms using likelihood ratio tests:
rustnovalmerandom=ranova(rustnovalmer)
###Extract fixed effects
rustnovafixef=fixef(rustnovalmer,add.dropped=TRUE)

####extract the means
rustnovals=ls_means(lmer(Severity~genotype*Stage+(1|Rep),rust))

## Least-Square means and pairwise differences:
(lsm <- ls_means(rustnovalmer))
SVRTlsm=ls_means(rustnovalmer, which = "genotype", pairwise = TRUE)
envlsm=ls_means(rustnovalmer, which = "Stage", pairwise = TRUE)

## ls_means also have plot and as.data.frame methods:
plot(lsm, which=c("genotype", "Stage"))
plsm=as.data.frame(lsm)
# Extract coefficient table:
coef(rustnova)
coef(summary(rustnova))
ls_means(rustnovalmer)
ranef(rustnovalmer)
?interaction.plot
###compute 2 way interaction plots
attach(sample)
Genotype=sample$Genotype
Stage=sample$Stage
interaction.plot(x.factor = Stage,
                 trace.factor =Genotype,
                 response = Severity,
                 fun = mean,
                 xlab = "Genotype",
                 ylab = "Severity",
                 trace.label = "Growth stage",
                 col = factor(Genotype),
                 lty = 05,
                 lwd=1,
                 main="Interaction plot for disease severity",
                 legend=TRUE,fixed = FALSE,
                 pch = 1,leg.bty = FALSE,
                 xtick = TRUE,axes = TRUE,type = "b")
####types l,p,b,o and c
# Type 3 anova table:
(anov <- anova(rustnovalmer, type="3"))
# Display tests/hypotheses for type 1, 2, and 3 ANOVA tables:
# (and illustrate effects of 'fractions' and 'names' arguments)
T1=show_tests(anova(rustnovalmer, type="1"))
T2=show_tests(anova(rustnovalmer, type="2"), fractions=TRUE, names=TRUE)
T3=show_tests(anova(rustnovalmer, fractions=TRUE))

###mean separation using boxplots
boxplot(Severity ~ Genotype*Stage,
        col=c("white","lightgray"),sample)
boxplot(HTP_Severity ~ Genotype*Stage,
        fill=genotype,rust)
library(ggplot2)
library(dplyr)
sample %>% ggplot(aes(x= Genotype, y=Severity)) +
  geom_boxplot()
sample %>% ggplot(aes(x= Genotype, y=Severity)) +
  geom_boxplot()+
  facet_grid(Genotype~Stage)+ geom_boxplot(varwidth = TRUE)
rust %>% ggplot(aes(x= genotype, y=HTP_Severity)) +
  geom_boxplot() +
  facet_grid(genotype~.)
rust %>% ggplot(aes(x= genotype, y=CSI )) +
  geom_jitter() +
  facet_grid(.~genotype)

#############Make it colourful##########################
###use aes_fill when using bopxplot and aes_colour when using geom_jitter
sample %>% ggplot(aes(x= Genotype, y=Severity)) +
  geom_boxplot(aes(colour=Genotype))+ggtitle("Severity")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(vjust = 0.5))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,size = 10))+
  theme(axis.text.y = element_text(angle = 0,size = 10))+
  theme(axis.title.x = element_text(size = 10))+
  theme(axis.title.y = element_text(size = 10))+
  facet_grid(.~Stage)
##################################make it unique
rust %>% ggplot(aes(x= genotype, y=HTP_Severity)) +
  geom_boxplot(aes(colour=genotype))+ggtitle("HTP Severity score")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(vjust = 0.5))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,size = 10))+
  theme(axis.text.y = element_text(angle = 0,size = 10))+
  theme(axis.title.x = element_text(size = 10))+
  theme(axis.title.y = element_text(size = 10))+
  facet_grid(.~Stage)+geom_violin()+order(decreasing=TRUE)
#######principal component analysis
treat=read.csv("microbes.csv",header = T)
attach(means)
means$Genotype=as.factor(means$Genotype)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggpubr)
df=means[,-1]
row.names(df)=as.factor(`Genotype`)
View(df)
res.PCA <- PCA(df, graph = TRUE)
print(res.PCA)
eig.val <- get_eigenvalue(res.PCA)
eig.val
head(eig.val)
eig.val=as.data.frame(eig.val)
eigenvals=write.csv(eig.val,"Eigen values for the treatments.csv")
write_xls(eig.val,"Eigen values for the parameters.xlsx")
fviz_screeplot(res.PCA)
var <- get_pca_var(res.PCA)
var$coord
var$cor
var$contrib
var$coord=as.data.frame(var$coord)
var$cor=as.data.frame(var$cor)
var$contrib=as.data.frame(var$contrib)

write_xlsx(var$coord,"Coordinates for the parameters.xlsx")

write_xlsx(var$cor,"correlations of the parameters .xlsx")

write_xlsx(var$contrib,"contribution variances for the parameters.xlsx")
write.table(var$coord,"Coord of variables.xls", row.names=TRUE, sep="\t")
write.table(var$cor,"Correlation_variables_dimensions.xls", row.names=TRUE, sep="\t")
write.table(var$contrib,"Contribution_variables.xls", row.names=TRUE, sep="\t")

# Plot of variables
fviz_pca_var(res.PCA, repel = TRUE)
# Contribution to the first dimension
a1=fviz_contrib(res.PCA, "var", axes = 1)
a1
# Contribution to the second dimension
b=fviz_contrib(res.PCA, "var", axes = 2)
b
c=fviz_contrib(res.PCA, "var", axes = 2)
c
ggarrange(a1,b,c,ncol=1,nrow = 2)

# Compute hierarchical clustering on principal components
res.hcpc <- HCPC(res.PCA, graph = FALSE)
res.hcpc
varA=fviz_pca_var(res.PCA,axes = c(1, 2),col.var = "cos2",
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                  repel = TRUE )+
  theme_minimal()
varA
varB=fviz_pca_var(res.PCA,axes = c(1, 2),col.var = "cos2",
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                  repel = TRUE )+
  theme_minimal()

varC=fviz_pca_var(res.PCA,axes = c(1, 2),col.var = "cos2",
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                  repel = TRUE )+
  theme_minimal()

ggarrange(varA,varB,ncol=2,legend = NULL,common.legend = TRUE)
c=fviz_dend(res.hcpc,
            cex = 0.8,                   # Label size
            # Color palette see ?ggpubr::ggpar
            rect = FALSE, rect_fill = FALSE, # Add rectangle around groups
            horiz=TRUE, # Rectangle color
            labels_track_height = 0.5       # Augment the room for labels
)
c
d=fviz_cluster(res.hcpc,axes = c(1, 2),
               repel = TRUE,           # Avoid label overlapping
               show.clust.cent = FALSE,geom=c("point","text"), # Show cluster centers
               # Color palette see ?ggpubr::ggpar
               ggtheme = theme_minimal(),ellipse.type = "convex",
               main = "Factor map"
)
d


res.hcpc
clustdata=res.hcpc$data.clust
write.table(  clustdata,"Dataset with the cluster of the individuals.xls", row.names=TRUE, sep="\t")
quanti=res.hcpc$desc.var$quanti
quanti$`1`
quanti$`2`
quanti$`3`
write.table(  quanti$`1`,"Description of the cluster 1 by the var.xls", row.names=TRUE, sep="\t")
write.table(  quanti$`2`,"Description of the cluster 2 by the var.xls", row.names=TRUE, sep="\t")
write.table(  quanti$`3`,"Description of the cluster 3 by the var.xls", row.names=TRUE, sep="\t")
res.hcpc$desc.axes$quanti
write.table(  res.hcpc$desc.axes$quanti,"Description of the clusters by the axes.xls", row.names=TRUE, sep="\t")

#######Correlation
library(corrplot)
df=Flowering[,4:15]
dp=peggging[,4:15]
dpod=pod_filling[,4:15]
dm=maturity[,4:15]
c=cor(dpod)
col=colorRampPalette(c("blue","white","red"))(20)
heatmap(x=c,col=col,symm = TRUE,order(decreasing = TRUE))
####heatmap split
heatmap(mat,name="map",row_split=2,column_split=3)
heatmap(mat,name="mat",row_split=splitdf,raw_gap=unit(5,"mm"))
#################Heatmaps#########
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
##############################ordinal logistics
data$NSP=as.ordered(data$NSP)
data$Tendance=as.factor(data$Tendance)
xtabs(~NSP+Tendance,data)
###Partitioning dataset
ind=sample(2,nrow(data),replace = TRUE,prob = c(0.8,0.2))
train=data[ind==1,]
test=data[ind==2,]
#####ordinal logistics regression model
library(MASS)
model=polr(NSP~LB+AC+FM,train,Hess = TRUE)
model=polr(NSP~.,train,Hess = TRUE)###if want to use all the variables
model=polr(NSP~.,-Max,train,Hess = TRUE)###if want to remove the variable or + if want to add

summary(model)
###calculating P values
(ctable=coef(summary(model)))
p=pnorm(abs(ctable[,"t value"]),lower.tail = FALSE)*2
(ctable=cbind(ctable,"p value"=p))
####prediction
pred=predict(model,train[1:5,],type='prob')
pred
pred=predict(model,train)###when using the whole data
pred
####confusion matrix and error training
(tab=table(pred,train$NSP))
1-sum(diag(tab))/sum(tab)
####confusion matrix and error for test data
pred1=predict(model,test)
(tab1=table(pred1,test$NSP))
1-sum(diag(tab1))/sum(tab1)
##############
#######path analysis

citation(package = "Raster")
