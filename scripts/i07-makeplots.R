
#
library(ggrepel)
library(ggplot2)
library(readxl)
library(patchwork)
library(RColorBrewer)

#get and prep data
supdata = read_xlsx('../data/41588_2023_1446_MOESM3_ESM.xlsx',sheet = 3,skip = 1)
meta.df = readRDS('../data/metadata.RDS')
meta.df$`Patient ID` = rownames(meta.df)
merged = merge(meta.df,supdata,by = 'Patient ID')

#generate the glm model
merged$multimf_cg2at = merged$`Regional difference in single molecule C>A frequency (per Mb of evaluable bases)`
merged$patient_type = factor(merged$patient_type, levels = c("Control", "Case"))
merged$patient_binary= 0
merged$patient_binary[merged$patient_type =='Case'] = 1
factor(merged$patient_type, levels = c("Control", "Case"))
lucas_glm <- glm(patient_type ~ multimf_cg2at, data = merged, family = "binomial")

#data from IUCPQ cohort (Gemini score prediction)
IUCPQ_data = read.csv('../outDir/gemini_score_hisat.csv')
IUCPQ_data$patient_type=factor(IUCPQ_data$X, levels = c("Control", "Case"))
IUCPQ_predicted_data = predict(lucas_glm,data.frame(patient_type = IUCPQ_data$patient_type, multimf_cg2at = IUCPQ_data$multimf_cg2at),type="response")
IUCPQ_data$IUCPQ_predicted_data = IUCPQ_predicted_data
IUCPQ_data$`Record ID` = sapply(strsplit(IUCPQ_data$id,'.',fixed = T),'[',5)
IUCPQ_data$ID = gsub('_cts','',IUCPQ_data$`Record ID`)
IUCPQ_data$`Record ID` = as.numeric(sapply(strsplit(IUCPQ_data$`Record ID`,'_',fixed = T),'[',2))
IUCPQ_data$provenance = 'IUCPQ'
IUCPQ_data$provenance[nchar(IUCPQ_data$ID) == 10] = 'GQ'

#IUCPQ_data = IUCPQ_data[IUCPQ_data$provenance == 'IUCPQ',]
#IUCPQ_data = IUCPQ_data[,c(5,11,12)]

#données sujet:
donnees_sujet = readxl::read_xlsx('../data/Données_sujets_cfDNA_2025-01-30.xlsx')
IUCPQ_data = merge(IUCPQ_data,donnees_sujet , by = 'Record ID')
IUCPQ_data$`Histological type`[IUCPQ_data$`Histological type`=='Adenocarcinoma'] = 'Adeno-\ncarcinoma'
IUCPQ_data$`Histological type`[IUCPQ_data$`Histological type`=='Squamous cell carcinoma'] = 'Squamous cell\ncarcinoma'

#plot
gemini_score_plot4 = ggplot() +
  geom_smooth(data = merged, aes(x = multimf_cg2at,y = `GEMINI score[C>A]`),col ='black',se = F,method = "glm", method.args = list(family = "binomial")) + 
  geom_point(data = merged, aes(x = multimf_cg2at,y = `GEMINI score[C>A]`,col = patient_type),pch = 21,size = 3,stroke =2,alpha = 1) + 
  scale_color_manual(values=RColorBrewer::brewer.pal(n = 4,'Set1')) +
  xlab('Regional difference in single molecule C>A frequency (per Mb of evaluable bases)') +
  ggtitle('high risk LUCAS cohort and logistic regression.') +
  theme_bw() + 
  theme(legend.position= 'inside',legend.position.inside = c(0.8, 0.2),legend.background = element_rect(fill="lightblue"))

gemini_score_plot3 = ggplot() +
  geom_smooth(data = merged, aes(x = multimf_cg2at,y = `GEMINI score[C>A]`),col ='black',se = F,method = "glm", method.args = list(family = "binomial")) + 
  geom_rect(aes(xmin = -7.5, xmax =-0.5, ymin =0, ymax=0.55),fill = 'blue',alpha = 0.4) +
  geom_rect(aes(xmin = -0.5, xmax =15, ymin =0.55, ymax=1),fill = 'red',alpha = 0.4) + 
  geom_point(data = IUCPQ_data,aes(x=multimf_cg2at,y=IUCPQ_predicted_data,col = `Histological type`),pch = 4,stroke =4, alpha = 1,inherit.aes =F) +
  geom_text_repel(data = IUCPQ_data,aes(x=multimf_cg2at,y=IUCPQ_predicted_data,label = ID),max.overlaps=200, hjust = -0.1) + 
  geom_text(data =IUCPQ_data,x=10, y=0.75, label="Cases",size =10) +
  geom_text(data =IUCPQ_data,x=-5, y=0.25, label="Control",size =10) +
  scale_color_manual(values=RColorBrewer::brewer.pal(n = 4,'Set1')[1:4]) +
  theme_bw() + 
  xlab('Regional difference in single molecule C>A frequency (per Mb of evaluable bases)') +
  ggtitle('logistic regression and fitted values of pilot experiment') +
  theme(legend.position= 'inside',legend.position.inside = c(0.8, 0.2),legend.background = element_rect(fill="lightblue"))


#save
pdf(file.path(paste0('../figures/Figure_Gemini_score.pdf')),width = 8,height = 10)
#(gemini_score_plot4 | gemini_score_plot3)  + plot_annotation(tag_levels = 'A')
gemini_score_plot3
dev.off()

###boxplots
variables = c('provenance','grade','Histological type','Predominant feature','Histopathological specification','Sex at birth')
gemini_box = list()
for(i in 1:6){
  gemini_box[[i]] = ggplot(IUCPQ_data,aes(x = .data[[variables[i]]], y = IUCPQ_predicted_data, fill = .data[[variables[i]]])) +
    geom_boxplot(outliers = F) +
    geom_jitter(color="black", size=1, alpha=0.9) + 
    geom_hline(yintercept = 0.55,linetype = 'dashed',col = 'red') +
    annotate("text", label = "Cancer detection threshold",x=0.9,y = 0.57, size = 3, colour = "red") +
    ylab('Gemini score') + 
    guides(fill="none") + 
    theme_bw()# +
  #  theme(axis.text.x = element_text(angle = 45, hjust=1))
}

#
pdf(file.path(paste0('../figures/Figure_Gemini_boxplot.pdf')),width = 4,height = 4)
#(gemini_box[[1]]|gemini_box[[2]]|gemini_box[[3]]) / (gemini_box[[4]]|gemini_box[[5]]|gemini_box[[6]]) + plot_layout(axes = "collect")
gemini_box[[3]]
dev.off()

 # geom_point(data = IUCPQ_data,aes(x=multimf_cg2at,y=IUCPQ_predicted_data,col = `Histological type`),pch = 4,stroke =4, alpha = 1,inherit.aes =F) +
 #  geom_text_repel(data = IUCPQ_data,aes(x=multimf_cg2at,y=IUCPQ_predicted_data,label = ID), hjust = -0.1) + 
  scale_color_manual(values=RColorBrewer::brewer.pal(n = 4,'Set1')[1:4]) +
  theme_bw() + 
  xlab('Regional difference in single molecule C>A frequency (per Mb of evaluable bases)') +
  ggtitle('logistic regression and fitted values of pilot experiment') +
  theme(legend.position= 'inside',legend.position.inside = c(0.8, 0.2),legend.background = element_rect(fill="lightblue"))


#save
#pdf(file.path(paste0('../data/Figure_Gemini_score.pdf')),width = 12,height = 10)
#(gemini_score_plot1 | gemini_score_plot4) / (gemini_score_plot2 | gemini_score_plot3) + plot_annotation(tag_levels = 'A')
#dev.off()


