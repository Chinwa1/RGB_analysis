###get the working directory
getwd()
###set the working directory
setwd("E:/Field trials/Groundnut/data")
###load the data, can use the import function
FA=read.csv("C:/Users/hchin/OneDrive/Atlantic Beach/Documents/PhD/DPhil/DPhil/RGB for indirect selection/Field trial data final/RGB_data_final/RGB_indexes_final.csv",header = TRUE)
###attach the data
attach(FA)
###check the structure of the data
str(FA)
###check the first 6 rows
head(symptomatic)
tail(symptomatic)
###check the class of the data
class(symptomatic)
getwd()
###make the columns factors
FA$AER=as.factor(FA$AER)
FA$Growth_stage=as.factor(FA$Growth_stage)
FA$Replication=as.factor(FA$Replication)
FA$Genotype=as.factor(FA$Genotype)
####Change integers to numeric###
FA$plants_with_frass=as.numeric(FA$plants_with_frass)
FA$Damage_rate=as.numeric(FA$Damage_rate)
###check for the normality of data
shapiro.test(yld$Grain_weight)
###plot a normality line
library(tidyselect)
library(ggplot2)
qqnorm(severity$symptomatic)
qqline(severity$symptomatic)
qqplot(GRV$Genotype,GRV$Severity)
# Using to explore the distribution of a variable
ggplot(severity, aes(sample = severity$symptomatic)) +
  stat_qq() +
  stat_qq_line()
ggplot(GRV, aes(sample = Severity, colour = factor(Genotype))) +
  stat_qq() +
  stat_qq_line()+
  theme_bw()+
  theme(element_blank())
###compute the model lmer
library(lme4)
library(lmerTest)
library(agricolae)
chin=lm(Severity~Genotype*AER+Growth_stage+Replication,FA)
anova(chin)
attributes(chin)
coef(chin)
summary(chin)
anova(chin)
View(severity)
chik=anova(chin)
write.csv(chik,"ANOVA FOR GY,csv")
####generate means using BIB
out<-BIB.test(FA$Agroecological_region, FA$Pattern, FA$GY, test="lsd")
mean=out$means

###export the output to excel
write.csv(mean,"means of GY.csv")
###Extract fixed effects
chim=fixef(chin,add.dropped=TRUE)
###export the fixed effects
write.csv(chim,"fixed effects.csv")

###check if the data is in data frame
class(chim)
####extract the means
chik=ls_means(lmer(h~Pattern*Agroecological_region+(1|Rep),FA))

write.csv(chik,"Means for yield lmer.csv")
###compute the ANOVA
chi=anova(chin)
## Anova-like table of random-effect terms using likelihood ratio tests:
chim=ranova(chin)
###Extract the anova table
write.csv(chim,"ANOVA for random effects with likelihood ratio tests GY.csv")

## Least-Square means and pairwise differences:
(lsm <- ls_means(chin))

GYlsm=ls_means(chin, which = "Pattern", pairwise = TRUE)

envlsm=ls_means(chin, which = "Agroecolocical_region", pairwise = TRUE)

###extract the table
write.csv(sitelsm,"LSM means for the Treatments based on CFU.csv")
write.csv(depthlsm,"LSM means for the depth based on CFU.csv")
write.csv(envlsm,"LSM means for the site based on CFU.csv")


## ls_means also have plot and as.data.frame methods:
## Not run: 
plot(lsm, which=c("Treatment", "Site","Depth"))
plsm=as.data.frame(lsm)

###extract the lsm table 
write.csv(plsm,"least square means of the CFU.csv")
# Extract coefficient table:
coef(chin)
coef(summary(chin))
ls_means(chin)
ranef(chim)
###compute 2 way interaction plots
interaction.plot(x.factor = GRV$Stage,
                 trace.factor =GRV$Genotype,
                 response = Severity,
                 fun = mean,
                 xlab = "Growth stage",
                 ylab = "Severity",
                 trace.label = "Genotype",
                 col = factor(Genotype),
                 lty = 4,
                 lwd=2.5,xtick = TRUE,
                 main="Disease severity")
with(GRV,{interaction.plot(Genotype,Stage,Severity)
  interaction.plot(Stage,Genotype,Severity,cex.axis=0.8)
  Stage=factor(Stage,levels = sort.list(tapply(Severity,Stage,mean)))
  interaction.plot(Stage,Genotype,Severity,col = factor(Genotype),lty = 3)
  
})

# Type 3 anova table:
(an <- anova(chin, type="3"))

# Display tests/hypotheses for type 1, 2, and 3 ANOVA tables:
# (and illustrate effects of 'fractions' and 'names' arguments)
T1=show_tests(anova(chin, type="1"))
T2=show_tests(anova(chin, type="2"), fractions=TRUE, names=TRUE)
T3=show_tests(an, fractions=TRUE)

c
#######principal component analysis
treat=read.csv("microbes.csv",header = T)
attach(treat)
treat$Genotype=as.factor(treat$Genotype)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggpubr)
df=treat[,-1]
row.names(df)=as.factor(`Genotype`)
View(df)
res.PCA <- PCA(df, graph = TRUE)
print(res.PCA)
eig.val <- get_eigenvalue(res.PCA)
eig.val
head(eig.val)
eig.val=as.data.frame(eig.val)
write_xlsx(eig.val,"Eigen values for the parameters.xlsx")
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

###Example 2: correlation between pH, variable 2 and other elements from soil.
library(agricolae)
floweringp=read.csv("floweringp.csv",row.names="Traits")
flowering=read.csv("E:/Glasshouse analysis/Leafspot/flowering.csv",row.names = "Trait")
pegging=read.csv("E:/Glasshouse analysis/Leafspot/pegging.csv",row.names = "Trait")
podfilling=read.csv("E:/Glasshouse analysis/Leafspot/podfilling.csv",row.names = "Trait")
maturity=read.csv("E:/Glasshouse analysis/Leafspot/maturity.csv",row.names = "Trait")
View(maturity)
attach(maturity)
floweringp=as.matrix(floweringp)
pegging=as.matrix(pegging)
podfilling=as.matrix(podfilling)
maturity=as.matrix(maturity)
write.csv(analysis1,file = "E:/Glasshouse analysis/Leafspot/maturity.csv")
analysis1<-correlation(flower[,1:13],method="pearson", alternative="two.sided")
View(analysis1)
write.csv(analysis1,file = "E:/Glasshouse analysis/Leafspot/pegging.csv")
library(ComplexHeatmap)
str(correl)
pegging=as.matrix(pegging)
h1=Heatmap(floweringp, name = "flowering",row_split = 2,column_split = 2)
h2=Heatmap(pegging, name = "pegging",row_split = 2,column_split = 2)
h3=Heatmap(podfilling, name = "podfilling",row_split = 2,column_split = 2)
h4=Heatmap(maturity, name = "maturity",row_split = 2,column_split = 2)
library(circlize)
col_fun = colorRamp2(c(0, 1, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))
Heatmap(floweringp, name = "Flowering", col = col_fun)

final=h1 + h2+h3 + h4
draw(final)
pdf(final,"Final heatmap.pdf")
######################################Heat maps#######
#install.packages('pheatmap')
# Create sample data ===================================================
set.seed(43)
data <- matrix(rnorm(500), 50, 10)
colnames(data) <- paste0("Sample_", 1:10)
rownames(data) <- paste0("Gene_", 1:50)

head(data)

# Annotations ===================================================

# create a data frame for column annotation
ann_df <- data.frame(Group = rep(c("Disease", "Control"), c(5, 5)),
                     Lymphocyte_count = rnorm(10))
row.names(ann_df) <- colnames(data)
head(ann_df)

gene_functions_df <- data.frame(gene_functions = rep(c('Oxidative_phosphorylation', 
                                                       'Cell_cycle',
                                                       'Immune_regulation',
                                                       'Signal_transduction',
                                                       'Transcription'), rep(10, 5)))
row.names(gene_functions_df) <- rownames(data)

ann_colors <- list(
  gene_functions = c("Oxidative_phosphorylation" = "#F46D43",
                     "Cell_cycle" = "#708238",
                     "Immune_regulation" = "#9E0142",
                     "Signal_transduction" = "beige", 
                     "Transcription" = "violet"), 
  Group = c("Disease" = "darkgreen",
            "Control" = "blueviolet"),
  Lymphocyte_count = brewer.pal(5, 'PuBu')
)



# Base heatmap ===================================================
library(pheatmap)
heat_plot <- pheatmap(data, 
                      col = brewer.pal(10, 'RdYlGn'), # choose a colour scale for your data
                      cluster_rows = T, cluster_cols = T, # set to FALSE if you want to remove the dendograms
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D',
                      annotation_row = gene_functions_df, # row (gene) annotations
                      annotation_col = ann_df, # column (sample) annotations
                      annotation_colors = ann_colors, # colours for your annotations
                      annotation_names_row = F, 
                      annotation_names_col = F,
                      fontsize_row = 10,          # row label font size
                      fontsize_col = 7,          # column label font size 
                      angle_col = 45, # sample names at an angle
                      legend_breaks = c(-2, 0, 2), # legend customisation
                      legend_labels = c("Low", "Medium", "High"), # legend customisation
                      show_colnames = T, show_rownames = F, # displaying column and row names
                      main = "Super heatmap with annotations") # a title for our heatmap
##############pheatmap test
library(RColorBrewer)
heat_plot <- pheatmap(floweringp, 
                      col = brewer.pal(10, 'RdYlGn'), # choose a colour scale for your data
                      cluster_rows = T, cluster_cols = T, # set to FALSE if you want to remove the dendograms
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D',
                      fontsize_row = 10,          # row label font size
                      fontsize_col = 10,          # column label font size 
                      angle_col = 45, # sample names at an angle
                      legend_breaks = c(-2, 0, 2), # legend customisation
                      legend_labels = c("Low", "Medium", "High"), # legend customisation
                      show_colnames = T, show_rownames = T, # displaying column and row names
                      main = "flowering heatmap ") # a title for our heatmap


# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(pheatmap) # for our heatmap
library(RColorBrewer) # for a colourful plot

# Set input path
path <- "~/Biostatsquid/Scripts/Heatmap/"
setwd(path)

# Create sample data ===================================================
set.seed(43)
data <- matrix(rnorm(500), 50, 10)
colnames(data) <- paste0("Sample_", 1:10)
rownames(data) <- paste0("Gene_", 1:50)

head(data)
# Annotations

# create a data frame for column annotation
ann_df <- data.frame(Group = rep(c("Disease", "Control"), c(5, 5)),
                     Lymphocyte_count = rnorm(10))
row.names(ann_df) <- colnames(data)
head(ann_df)

gene_functions_df <- data.frame(gene_functions = rep(c('Oxidative_phosphorylation', 
                                                       'Cell_cycle',
                                                       'Immune_regulation',
                                                       'Signal_transduction',
                                                       'Transcription'), rep(10, 5)))
row.names(gene_functions_df) <- rownames(data)
####Step 1: Create a clustered heatmap
## Get heatmap ===================================================
pheatmap(data)
# Clustering ===================================================
pheatmap(data, 
         cluster_rows = T, cluster_cols = T, # set to FALSE if you want to remove the dendograms
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D')
# Customisation ===================================================

## Adding a title
pheatmap(data, 
         main = "Super cool heatmap")

## Showing rows and columns
pheatmap(data,
         main = "Super cool heatmap",
         show_colnames = T, show_rownames = T,
         number_color = "black", 
         fontsize_number = 8)

## Showing values
pheatmap(data,
         fontsize_col = 10,
         fontsize_row = 10,
         display_numbers = TRUE,
         number_color = "black", 
         fontsize_number = 6,#
         border_color = "black") # default is grey60
## Cell colours
pheatmap(data,
         border_color = "black", # default is grey60
         number_color = "black", 
         fontsize_number = 8,
         col = brewer.pal(10, 'RdYlGn')) # https://r-graph-gallery.com/38-rcolorbrewers-palettes.html
## Legend customisation
# add a custom legend
pheatmap(data, 
         legend_breaks = c(-2, 0, 2),
         legend_labels = c("Low", "Medium", "High"))

# or remove the legend
pheatmap(data, 
         legend = FALSE)
# Split heatmap clusters
pheatmap(data, 
         cutree_rows = 2, cutree_cols = 4)
pheatmap(data, 
         border_color = FALSE,      # no border to cell
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size 
         angle_col = 45,             # angle for column labels
         na_col = "black",           # color of the cell with NA values
         legend = FALSE#,            # to draw legend or not (TRUE/FALSE)
)
ann_colors <- list(
  gene_functions = c("Oxidative_phosphorylation" = "#F46D43",
                     "Cell_cycle" = "#708238",
                     "Immune_regulation" = "#9E0142",
                     "Signal_transduction" = "beige", 
                     "Transcription" = "violet"), 
  Group = c("Disease" = "darkgreen",
            "Control" = "blueviolet"),
  Lymphocyte_count = brewer.pal(5, 'PuBu')
)


pheatmap(data, 
         col = brewer.pal(10, 'RdYlGn'),
         annotation_row = gene_functions_df, 
         annotation_col = ann_df, 
         annotation_colors = ann_colors,
         main = "Super heatmap with annotations") 
heat_plot <- pheatmap(data, 
                      col = brewer.pal(10, 'RdYlGn'), # choose a colour scale for your data
                      cluster_rows = T, cluster_cols = T, # set to FALSE if you want to remove the dendograms
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D',
                      annotation_row = gene_functions_df, # row (gene) annotations
                      annotation_col = ann_df, # column (sample) annotations
                      annotation_colors = ann_colors, # colours for your annotations
                      annotation_names_row = F, 
                      annotation_names_col = F,
                      fontsize_row = 10,          # row label font size
                      fontsize_col = 7,          # column label font size 
                      angle_col = 45, # sample names at an angle
                      legend_breaks = c(-2, 0, 2), # legend customisation
                      legend_labels = c("Low", "Medium", "High"), # legend customisation
                      show_colnames = T, show_rownames = F, # displaying column and row names
                      main = "Super heatmap with annotations") # a title for our heatmap
# Save it -----------------------------------------------------------
pdf(heat_plot, height = 10, width = 8)
heat_plot
dev.off()
#####################
####Heatmap correlation
attach(Correlflower)
library(corrplot)
c=cor(Correlflower[,2:10])
col=colorRampPalette(c("blue","white","red"))
heatmap(x=c,color = col)
######################
attach(correlpeg)
library(corrplot)
d=cor(correlpeg[,2:10])
cols=colorRampPalette(c("blue","white","red"))
heatmap(x=c,color=col,symm = T)
###################
attach(correlpod)
library(corrplot)
e=cor(correlpod[,2:10])
cols=colorRampPalette(c("blue","white","red"))
heatmap(x=e,color=col,symm = T)

########################
attach(correlmat)
library(corrplot)
f=cor(correlmat[,2:10])
cols=colorRampPalette(c("blue","white","red"))
flowering=as.matrix(flowering)
heatmap(flowering,color=col,symm = T)
