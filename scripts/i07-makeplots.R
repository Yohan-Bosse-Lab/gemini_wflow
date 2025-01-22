
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
IUCPQ_data = read.csv('../outDir/gemini_score.csv')
IUCPQ_data$patient_type=factor(IUCPQ_data$X, levels = c("Control", "Case"))
IUCPQ_predicted_data = predict(lucas_glm,data.frame(patient_type = IUCPQ_data$patient_type, multimf_cg2at = IUCPQ_data$multimf_cg2at),type="response")
IUCPQ_data$IUCPQ_predicted_data = IUCPQ_predicted_data
IUCPQ_data$patient_ID = sapply(strsplit(IUCPQ_data$id,'.',fixed = T),'[',5)
IUCPQ_data$patient_ID = gsub('_cts','',IUCPQ_data$patient_ID)


#plot
gemini_score_plot4 = ggplot() +
  geom_smooth(data = merged, aes(x = multimf_cg2at,y = `GEMINI score[C>A]`),col ='black',se = F,method = "glm", method.args = list(family = "binomial")) + 
  geom_point(data = merged, aes(x = multimf_cg2at,y = `GEMINI score[C>A]`,col = patient_type),pch = 21,size = 3,stroke =2,alpha = 1) + 
  scale_color_manual(values=RColorBrewer::brewer.pal(n = 4,'Set1')) +
  xlab('Regional difference in single molecule C>A frequency (per Mb of evaluable bases)') +
  ggtitle('high risk LUCAS cohort and logistic regression.') +
  theme_bw() + 
  theme(legend.position  = 'inside',legend.position.inside = c(0.8,0.25))

gemini_score_plot3 = ggplot() +
  geom_smooth(data = merged, aes(x = multimf_cg2at,y = `GEMINI score[C>A]`),col ='black',se = F,method = "glm", method.args = list(family = "binomial")) + 
  geom_point(data = IUCPQ_data,aes(x=multimf_cg2at,y=IUCPQ_predicted_data),col = 'goldenrod',pch = 4,stroke =4, inherit.aes =F) +
  geom_text_repel(data = IUCPQ_data,aes(x=multimf_cg2at,y=IUCPQ_predicted_data,label = patient_ID), hjust = -0.1) + 
  scale_color_manual(values=RColorBrewer::brewer.pal(n = 4,'Set1')) +
  theme_bw() + 
  xlab('Regional difference in single molecule C>A frequency (per Mb of evaluable bases)') +
  ggtitle('logistic regression and fitted values of pilot experiment\n(inversing cancer/control fixed bins)') 


#save
pdf(file.path(paste0('../data/Figure_Gemini_score.pdf')),width = 12,height = 6)
(gemini_score_plot4 | gemini_score_plot3)  + plot_annotation(tag_levels = 'A')
dev.off()


#save
#pdf(file.path(paste0('../data/Figure_Gemini_score.pdf')),width = 12,height = 10)
#(gemini_score_plot1 | gemini_score_plot4) / (gemini_score_plot2 | gemini_score_plot3) + plot_annotation(tag_levels = 'A')
#dev.off()


