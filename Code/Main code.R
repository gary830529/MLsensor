pacman::p_load(plyr,dplyr,reshape2,stringr,lubridate,readxl,openxlsx,janitor,ggplot2,ggpmisc,ggpubr,ggsignif,gplots,data.table,magrittr,grid,
               tidyverse,rstatix,PtProcess,gtable,gridExtra,cowplot,tiff,pwr,MOTE,kableExtra,tinytex,knitr,quantmod,class,caret,gmodels,
               officedown,officer,flextable,car,C50,GGally,ggResidpanel,ggfortify,caretEnsemble,randomForest,corrplot,neuralnet,lemon,Hmisc,scales,
               glue,ggtext,png,gtools,ggrepel,rvg,gdata,scales,nnet,xgboost,correlation,psych,corrgram,jpeg,qrnn,pls,kernlab)
colx <- colorRampPalette(c("#16068a","#9e189d","#fdb32e"))
colx2 <- colorRampPalette(c("#007abf","#f7b81a"))

#### Start ####
Raw <- read.csv(file = 'G:/My Drive/R project/GitHub/MLsensor/Raw_data/Data_ML.csv',header = T, sep = ",")
Raw$Date <- ymd(Raw$Date)
Raw[Raw < 0] <- NaN
Raw <- Raw[Raw$Date <= ymd("2020-03-10"), ]
Raw <-  Raw %>% filter(Group %in% c("S1","S2","S3","S4","S5"))
colnames(Raw) <- c("X","Group","Date","COD","sCOD","TSS","E.coli","Color","Turb","EC","pH","NH4","NO3","Temp","VSS")


#### Clean Data ####
Dates_after_events <- as.Date(c("2019-01-15", "2019-01-22", "2019-08-21", "2019-08-27", "2020-01-28", "2020-02-04"#, "2020-02-11"
)) #"2019-03-07" weird pH reading
Raw$Event <-  ifelse (Raw$Date %in% Dates_after_events, "Post-Event", "Normal")
Raw$BP <-  ifelse (Raw$Group %in% c("S1","S2"), 0, 1)
Raw <- Raw %>% arrange(Date,Group)
Raw$Group <-  ifelse (Raw$Group == c("S1","S2","S3","S4","S5"), c("Influent","AnMBR","Permeate","Post-NCS","Effluent"), 0)
Raw$Group <- factor(Raw$Group , levels=c("Influent","AnMBR","Permeate","Post-NCS","Effluent"))
Raw <- subset(Raw, select = c("Group","Date","COD", "sCOD","TSS","E.coli","Color","Turb","EC","pH","NH4","NO3","Temp","Event","BP","VSS"))

Raw$VSS <- Raw$VSS*1000
Raw$TSS <- Raw$TSS*1000
Raw$TSS <- ifelse(Raw$TSS <= 0, 0.5, Raw$TSS)
Raw$COD <- ifelse(Raw$COD <= 0, 0.5, Raw$COD)
Raw$E.coli <- ifelse(Raw$E.coli <= 0, 1, Raw$E.coli)
Raw$pCOD <- ifelse(Raw$COD - Raw$sCOD <= 0, 0, Raw$COD - Raw$sCOD)
Raw <- Raw[c("Group","Date","COD", "sCOD","pCOD","TSS","E.coli","Color","Turb","EC","pH","NH4","NO3","Temp","Event","BP","VSS")]
WQ <- Raw


#### Mean Absolute Percentage Error #### 
MAPE <- function(pred.,actu.,...){
  mean(abs((ifelse(actu. == 0, (actu.+pred.)/2, actu.)-pred.)/ifelse(actu. > pred., ifelse(actu. == 0, (actu.+pred.)/2, actu.),pred.) ))*100
  #    ifelse(abs((actu.-pred.)/actu.)>1,NA,abs((actu.-pred.)/actu.)), na.rm=TRUE) * 100
}


#### Normality ####
Normality <- WQ %>% subset(select=-c(Date,E.coli,Event,BP,sCOD,pCOD)) %>%
  reshape2::melt(id.vars="Group") %>% group_by(Group,variable) %>%
  summarise(statistic = shapiro.test(value)$statistic,
            p.value = shapiro.test(value)$p.value) 
Normality$normal.distribution <- ifelse(Normality$p.value > 0.05, "Yes", "No")

Normality2 <- WQ %>% subset(select=-c(Date,E.coli,Event,BP,sCOD,pCOD)) %>% 
  mutate_at(c("COD","TSS","Color","Turb","EC","pH","NH4","NO3","Temp"), ~(scale(.) %>% as.vector)) %>%
  reshape2::melt(id.vars="Group") %>% group_by(Group,variable) %>%
  summarise(statistic = shapiro.test(value)$statistic,
            p.value = shapiro.test(value)$p.value) 
Normality2$normal.distribution <- ifelse(Normality2$p.value > 0.05, "Yes", "No")


#### Boxplot ####
Bplot <- function(titl.,...){ 
  list(stat_compare_means(aes(group = Event), label = "p.signif", method = "t.test",method.args = list(alternative = "greater"),size = 4), #, vjust = 25.5
       scale_colour_manual(values = c("Normal"="#00000066", "Post-Event"="#b82c2c66")),
       coord_flex_cart(bottom=brackets_horizontal(), left=capped_vertical('none')),theme_bw(),
       theme(text = element_text(size=12),legend.title=element_blank(),legend.position = "none",
             axis.text.x=element_text(angle=20, vjust=0.7),
             panel.border=element_blank(), axis.line = element_line(),legend.background = element_blank(),
             strip.background = element_blank(), strip.text.x = element_blank(),plot.margin = ggplot2::margin(15,0,-10,0, "pt")),
       labs(x = element_blank(), y = titl.),...)}

### COD Boxplot ### 
C.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "COD", add = "jitter", legend = "none", add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("COD (mg/L)",ggforce::facet_row(vars(BP), scales = 'free', space = 'free')) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))
  ggboxplot(WQ, x = "Group", y = "COD", add = "jitter", legend = "none", color = "Event", shape = "Event", add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("COD (mg/L)"#,ggforce::facet_row(vars(BP), scales = 'free', space = 'free')
  ) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff")) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format (10^.x)))

### TSS Boxplot ### 
T.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "TSS", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("TSS (mg/L)",ggforce::facet_row(vars(BP), scales = 'free', space = 'free')) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))
  ggboxplot(WQ, x = "Group", y = "TSS", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("TSS (mg/L)"#,ggforce::facet_row(vars(BP), scales = 'free', space = 'free')
  ) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff")) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format (10^.x)))

### E.coli Boxplot ### 
E.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "E.coli", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("E.coli (MPN/100ml)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff")) +
  #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format (10^.x)))
  ggboxplot(WQ, x = "Group", y = "E.coli", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("E.coli (MPN/100ml)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff")) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format (10^.x)))

### Color Boxplot ### 
co.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "Color", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("Color (Pt/Co)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))+ylim(NA,7500)
  ggboxplot(WQ, x = "Group", y = "Color", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("Color (Pt/Co)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))+ylim(NA,7500)

### Turbidity Boxplot ### 
tu.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "Turb", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("Turbidity (NTU)",ggforce::facet_row(vars(BP), scales = 'free', space = 'free')) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))
  ggboxplot(WQ, x = "Group", y = "Turb", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("Turbidity (NTU)",ggforce::facet_row(vars(BP), scales = 'free', space = 'free')) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))

### EC Boxplot ### 
ec.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "EC", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("EC (µs/cm)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))
  ggboxplot(WQ, x = "Group", y = "EC", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("EC (µs/cm)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))

### pH Boxplot ### 
ph.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "pH", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("pH") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))
  ggboxplot(WQ, x = "Group", y = "pH", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("pH") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))

###  Temp Boxplot ### 
tbp <- WQ %>% filter(Group == c("Influent","AnMBR"))
tbp$Group <- ifelse(tbp$Group == "Influent", "Ambient air", "AnMBR")
temp.Bplot <- 
  #ggboxplot(tbp, x = "Group", y = "Temp", add = "jitter", legend = "none",add.params = list(size = .75), order = c("Ambient air", "AnMBR")) + Bplot("°C") + theme(legend.position = "right") + 
  #guides(color = guide_legend(label.position = "bottom"))
  ggboxplot(tbp, x = "Group", y = "Temp", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = .75),
            order = c("Ambient air", "AnMBR")) + Bplot("°C") + theme(legend.position = "right") + 
  guides(color = guide_legend(label.position = "bottom"))

### Boxplot ### 
Figure3 <- ggarrange(C.Bplot, T.Bplot, E.Bplot + theme(legend.position = "right") +
                          guides(color = guide_legend(label.position = "bottom")),
                        nrow=1, ncol=3, align = "v", labels = "AUTO", common.legend = T, widths = c(1,1,1), legend="bottom")
Figure3
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/Figure3.tiff",plot=Figure3,width = 10, height = 4, dpi = 300)


#### Z-score ####
WQ_z <- WQ %>% select("Date", "Group", "COD", "TSS", "Color", "Turb", "EC", "pH") %>% #-"SCOD", -"PCOD", -"Event", -"BP", -"VSS", -"E.coli", -"Temp") %>% 
  filter(Group %in% c("AnMBR","Permeate","Post-NCS","Effluent")) %>% 
  group_by(Group) %>%
  mutate(across(COD:pH, scale)) %>% ungroup() 

WQ_z$z_score <- rowMeans(WQ_z[,3:8], na.rm=TRUE)

date_z_score <- WQ_z %>%
  group_by(Date) %>%
  summarise(date_z = mean(z_score, na.rm=TRUE)) %>% ungroup()

Z_plot <- ggplot(date_z_score, aes(x = Date, y = date_z)) +
  geom_line(size=0.5,alpha=.7,linetype=3) + #geom_vline(xintercept = Dates_after_events, linetype = "dashed", color= "red") + 
  geom_point(data = subset(date_z_score, Date %nin% c(Dates_after_events)),aes(x = Date, y = date_z), color= "black")+
  geom_point(data = subset(date_z_score, Date %in% c(Dates_after_events)),aes(x = Date, y = date_z), color= "red")+
  theme_bw() + labs(x = element_blank(), y = "Z-score") + 
  scale_x_date(date_breaks = "3 month",date_labels =  "%b\n%Y",expand = c(0,0),limits = c(ymd("2018-10-01"),ymd("2020-04-01")))+
  
  geom_hline(aes(yintercept= mean(date_z_score$date_z, na.rm = T)+sd(date_z_score$date_z, na.rm = T)*1.70), color= "#ff0000",size=0.8,alpha=0.5, linetype = "dashed")+
  geom_hline(aes(yintercept= mean(date_z_score$date_z, na.rm = T)-sd(date_z_score$date_z, na.rm = T)*1.70), color= "#ff0000",size=0.8,alpha=0.5, linetype = "dashed")+
  
  annotate("rect", xmin = as.Date("2018-12-05"), xmax = as.Date("2019-01-07"), ymin = -Inf, ymax = Inf, alpha  = 0.2)+
  annotate("text", x = as.Date("2018-12-24"), y = 0.8, label = "Shutdown", angle = -90, size = 4)+
  
  annotate("rect", xmin = as.Date("2019-11-30"), xmax = as.Date("2020-01-22"), ymin = -Inf, ymax = Inf, alpha  = 0.2)+
  annotate("text", x = as.Date("2019-12-28"), y = 0.8, label = "NCS maintenance 2 \n Shutdown", angle = -90, size = 4)+
  
  annotate("rect", xmin = as.Date ("2019-08-08"), xmax = as.Date("2019-08-21"), ymin = -Inf, ymax = Inf, alpha  = 0.2)+
  annotate("text", x = as.Date("2019-08-16"), y = 0.8, label = "NCS maintenance 1", angle = -90, size = 4)+
  theme(text = element_text(size=15),legend.position="bottom",
        plot.margin = ggplot2::margin(10,20,0,0, "pt"))

Z_plot
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/FigureS2.tiff",plot=Z_plot,width = 6, height = 4, dpi = 300)


#### Restart remove data ####
Dates_after_events2 <- as.Date(c("2019-01-15", "2019-01-22"))
WQ2 <- WQ %>% filter(Date %nin% Dates_after_events2)
WQ2 <- WQ2[complete.cases(WQ2), ]
WQ2 <- WQ2 %>% subset(select=-c(Date,Event,BP,VSS))
WQ2$Group <- as.numeric(as.factor(WQ2$Group))


#### Data for E.coli Classification ####
WQ2$E.coli2 <- ifelse(WQ2$E.coli >= 10^6, "H", ifelse(WQ2$E.coli < 10^6 & WQ2$E.coli >= 10^4, "M", WQ2$E.coli)) 
WQ2$E.coli2 <- ifelse(WQ2$E.coli < 10^4 & WQ2$E.coli >= 2, "L", WQ2$E.coli2)
WQ2$E.coli2 <- ifelse(WQ2$E.coli < 2, "LDL", WQ2$E.coli2)
WQ2$E.coli2 <- factor(WQ2$E.coli2 , levels=c("H","M","L","LDL"))


#### Correlation ####
correlation::correlation(subset(WQ2, select = -c(Group)),include_factors = TRUE, method = "auto")
correlation::correlation(WQ2,include_factors = TRUE, method = "auto")

WQ2_m <- WQ2 %>% subset(select=-c(Group,E.coli2,sCOD,pCOD))
testRes <- cor.mtest(WQ2_m, conf.level = 0.95)
col <- colorRampPalette(c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"))


file_path= "Correlation matrix.tif"
png(width=6,height=6,units="in",res=600, file=file_path, type = "cairo")

corrplot::corrplot(cor(WQ2_m), method = 'color', type = 'upper', pch.cex = 0.75,
                   tl.col = 'black',tl.pos="d", col = col(100),tl.cex=0.8,
                   addgrid.col = 'grey20', addCoef.col = 'black', number.cex = 0.8)

corrplot::corrplot(cor(WQ2_m), p.mat = testRes$p, col = col(100), diag = TRUE, type = 'lower',
                   sig.level = c(0.0001,0.001, 0.01, 0.05), pch.cex = 0.75, tl.col = 'black',tl.pos="d",tl.cex=0.8,
                   insig = 'label_sig', pch.col = 'black', cl.pos = 'n',addgrid.col = 'grey20', add=T)
dev.off()

###  Include SCOD & PCOD ### 
upper_fn <- function(data, mapping, method="pearson", use="pairwise", ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, method=method, use=use)
  colFn <- colorRampPalette(c("#4477AA", "#ffffff", "#BB4444"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  ggally_cor(data = data, mapping = mapping, size = 2, color = '#050505', stars = T, ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

diag_fn <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(alpha=.2,...)+theme_bw()
}

lower_fn <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = '#050505', alpha=0.3, size=0.5) +
    geom_smooth(color = '#ff8a8a', method='lm', size=0.5,se=F,...)+theme_bw()
}

WQ2.Cor <- ggpairs(subset(subset(WQ2, Group %in% c("1")), select = -c(Group,E.coli2)), 
                      upper = list(continuous = upper_fn), lower = list(continuous = lower_fn), diag = list(continuous = diag_fn))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.text.y = element_text(size=5),text = element_text(size=6))
WQ2.Cor
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/Correlation.jpg",plot=WQ2.Cor,width = 6, height = 6, dpi = 300)


#### Feature Elimination Analysis ####
rctrl <- rfeControl(method = "cv", number=5, functions = caretFuncs, verbose = FALSE)

c.rfe.fun <- function(mod.,inputT,...){
  caret::rfe(COD ~ Group+Turb+pH+NH4+NO3+Temp+Color+EC, data=inputT, sizes=c(5:8), method=mod.,rfeControl = rctrl,...)
}

t.rfe.fun <- function(mod.,inputT,...){
  caret::rfe(TSS ~ Group+Turb+pH+NH4+NO3+Temp+Color+EC, data=inputT, sizes=c(5:8), method=mod.,rfeControl = rctrl,...)
}

e.rfe.fun <- function(mod.,inputT,...){
  caret::rfe(E.coli ~ Group+Turb+pH+NH4+NO3+Temp+Color+EC, data=inputT, sizes=c(5:8), method=mod.,rfeControl = rctrl,...)
}


#### COD REF Models ####
rfe.C.lm <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_RFE/C_rfe_lm.rds")
#c.rfe.fun("lm",WQ.Cp)
rfe.C.pls <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_RFE/C_rfe_pls.rds")
#c.rfe.fun("pls",WQ.Cp,tuneLength = 5)
rfe.C.rf <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_RFE/C_rfe_rf.rds")
#c.rfe.fun("rf",WQ.Cp)
rfe.C.knn <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_RFE/C_rfe_knn.rds")
#c.rfe.fun("knn",WQ.Cp)
rfe.C.svm <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_RFE/C_rfe_svm.rds")
#c.rfe.fun("svmRadial",WQ.Cp)
rfe.C.cub <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_RFE/C_rfe_cub.rds")
#c.rfe.fun("cubist",WQ.Cp)
rfe.C.net <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_RFE/C_rfe_rnn.rds")
#c.rfe.fun("qrnn",WQ.Cp)

#saveRDS(rfe.C.net,"G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_RFE/C_rfe_rnn.rds")

### COD REF Plot ### 
rfe.C.lm_Plot <- ggplot(rfe.C.lm,aes(rfe.C.lm$variables,rfe.C.lm$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "COD-LM")+
  annotate("label", label = paste("Variable: ",predictors(rfe.C.lm)[1],",",
                                  predictors(rfe.C.lm)[2],",", predictors(rfe.C.lm)[3],",",
                                  predictors(rfe.C.lm)[4],",", predictors(rfe.C.lm)[5]), x = 6.5,
           y = max(rfe.C.lm$results$RMSE+rfe.C.lm$results$RMSESD)+100, size = 1.8, colour = "black")+ 
  ylim(200,1000)+
  geom_errorbar(aes(ymin=(rfe.C.lm$results$RMSE-rfe.C.lm$results$RMSESD), ymax=(rfe.C.lm$results$RMSE+rfe.C.lm$results$RMSESD)), width=0.1,alpha=0.8)

rfe.C.pls_Plot <- ggplot(rfe.C.pls,aes(rfe.C.pls$variables,rfe.C.pls$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8)) + labs(y = "RMSE",title = "COD-PLS")+
  annotate("label", label = paste("Variable: ",predictors(rfe.C.pls)[1],",",
                                  predictors(rfe.C.pls)[2],",", predictors(rfe.C.pls)[3],",",
                                  predictors(rfe.C.pls)[4],",", predictors(rfe.C.pls)[5],",",
                                  predictors(rfe.C.pls)[6]), x = 6.5,
           y = max(rfe.C.pls$results$RMSE+rfe.C.pls$results$RMSESD)+100, size = 1.8, colour = "black")+ 
  ylim(200,1000)+
  geom_errorbar(aes(ymin=(rfe.C.pls$results$RMSE-rfe.C.pls$results$RMSESD), ymax=(rfe.C.pls$results$RMSE+rfe.C.pls$results$RMSESD)), width=0.1,alpha=0.8)

rfe.C.rf_Plot <- ggplot(rfe.C.rf,aes(rfe.C.rf$variables,rfe.C.rf$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8)) + labs(y = "RMSE",title = "COD-RF")+
  annotate("label", label = paste("Variable: ",predictors(rfe.C.rf)[1],",",
                                  predictors(rfe.C.rf)[2],",", predictors(rfe.C.rf)[3],",",
                                  predictors(rfe.C.rf)[4],",", predictors(rfe.C.rf)[5]),x = 6.5,
           y = max(rfe.C.rf$results$RMSE+rfe.C.rf$results$RMSESD)+100, size = 2, colour = "black")+ 
  ylim(200,1000)+
  geom_errorbar(aes(ymin=(rfe.C.rf$results$RMSE-rfe.C.rf$results$RMSESD), ymax=(rfe.C.rf$results$RMSE+rfe.C.rf$results$RMSESD)), width=0.1,alpha=0.8)

rfe.C.knn_Plot <- ggplot(rfe.C.knn,aes(rfe.C.knn$variables,rfe.C.knn$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8)) + labs(y = "RMSE",title = "COD-KNN")+
  annotate("label", label = paste("Variable: ",predictors(rfe.C.knn)[1],",",
                                  predictors(rfe.C.knn)[2],",", predictors(rfe.C.knn)[3],",",
                                  predictors(rfe.C.knn)[4],",", predictors(rfe.C.knn)[5]),x = 6.5,
           y = max(rfe.C.knn$results$RMSE+rfe.C.knn$results$RMSESD)+100, size = 2, colour = "black")+
  ylim(200,1000)+
  geom_errorbar(aes(ymin=(rfe.C.knn$results$RMSE-rfe.C.knn$results$RMSESD), ymax=(rfe.C.knn$results$RMSE+rfe.C.knn$results$RMSESD)), width=0.1,alpha=0.8)

rfe.C.svm_Plot <- ggplot(rfe.C.svm,aes(rfe.C.svm$variables,rfe.C.svm$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8)) + labs(y = "RMSE",title = "COD-SVR")+
  annotate("label", label = paste("Variable: ",predictors(rfe.C.svm)[1],",",
                                  predictors(rfe.C.svm)[2],",", predictors(rfe.C.svm)[3],",",
                                  predictors(rfe.C.svm)[4],",", predictors(rfe.C.svm)[5],",",
                                  predictors(rfe.C.svm)[6],",", predictors(rfe.C.svm)[7],",",
                                  predictors(rfe.C.svm)[8]),x = 6.5,
           y = max(rfe.C.svm$results$RMSE+rfe.C.svm$results$RMSESD)+100, size = 1.7, colour = "black")+
  ylim(100,1000)+
  geom_errorbar(aes(ymin=(rfe.C.svm$results$RMSE-rfe.C.svm$results$RMSESD), ymax=(rfe.C.svm$results$RMSE+rfe.C.svm$results$RMSESD)), width=0.1,alpha=0.8)

rfe.C.cub_Plot <- ggplot(rfe.C.cub,aes(rfe.C.cub$variables,rfe.C.cub$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8)) + labs(y = "RMSE",title = "COD-CUB")+
  annotate("label", label = paste("Variable: ",predictors(rfe.C.cub)[1],",",
                                  predictors(rfe.C.cub)[2],",", predictors(rfe.C.cub)[3],",",
                                  predictors(rfe.C.cub)[4],",", predictors(rfe.C.cub)[5],",",
                                  predictors(rfe.C.cub)[6],",", predictors(rfe.C.cub)[7]),x = 6.5,
           y = max(rfe.C.cub$results$RMSE+rfe.C.cub$results$RMSESD)+100, size = 1.8, colour = "black")+
  ylim(200,1000)+
  geom_errorbar(aes(ymin=(rfe.C.cub$results$RMSE-rfe.C.cub$results$RMSESD), ymax=(rfe.C.cub$results$RMSE+rfe.C.cub$results$RMSESD)), width=0.1,alpha=0.8)

rfe.C.net_Plot <- ggplot(rfe.C.net,aes(rfe.C.net$variables,rfe.C.net$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8)) + labs(y = "RMSE",title = "COD-QRNN")+
  annotate("label", label = paste("Variable: ",predictors(rfe.C.net)[1],",",
                                  predictors(rfe.C.net)[2],",", predictors(rfe.C.net)[3],",",
                                  predictors(rfe.C.net)[4],",", predictors(rfe.C.net)[5],",",
                                  "pH"),x = 6.5,
           y = max(rfe.C.net$results$RMSE+rfe.C.net$results$RMSESD)+100, size = 1.8, colour = "black")+
  ylim(200,1000)+
  geom_errorbar(aes(ymin=(rfe.C.net$results$RMSE-rfe.C.net$results$RMSESD), ymax=(rfe.C.net$results$RMSE+rfe.C.net$results$RMSESD)), width=0.1,alpha=0.8)

rfe.C_Plot <- ggarrange(rfe.C.net_Plot, rfe.C.pls_Plot, rfe.C.knn_Plot,rfe.C.svm_Plot,rfe.C.rf_Plot,rfe.C.cub_Plot,
                        labels = "AUTO", nrow=3, ncol = 2, align = "v")
rfe.C_Plot
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/FigureS3-1.jpg",plot=rfe.C_Plot,width = 6, height = 4, dpi = 300)


#### TSS REF Models ####
rfe.T.lm <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_RFE/T_rfe_lm.rds")
#t.rfe.fun("lm",WQ.Cp)
rfe.T.pls <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_RFE/T_rfe_pls.rds")
#t.rfe.fun("pls",WQ.Cp,tuneLength = 5)
rfe.T.rf <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_RFE/T_rfe_rf.rds")
#t.rfe.fun("rf",WQ.Cp)
rfe.T.knn <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_RFE/T_rfe_knn.rds")
#t.rfe.fun("knn",WQ.Cp)
rfe.T.svm <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_RFE/T_rfe_svm.rds")
#t.rfe.fun("svmRadial",WQ.Cp)
rfe.T.cub <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_RFE/T_rfe_cub.rds")
#t.rfe.fun("cubist",WQ.Cp)
rfe.T.net <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_RFE/T_rfe_rnn.rds")
#t.rfe.fun("qrnn",WQ.Cp)

#saveRDS(rfe.T.net,"G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_RFE/T_rfe_rnn.rds")

###  TSS REF Plot ### 
rfe.T.net_Plot <- ggplot(rfe.T.net,aes(rfe.T.net$variables,rfe.T.net$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "TSS-QRNN")+
  annotate("label", label = paste("Variable: ",predictors(rfe.T.net)[1],",",
                                  predictors(rfe.T.net)[2],",", predictors(rfe.T.net)[3],",",
                                  predictors(rfe.T.net)[4],",", predictors(rfe.T.net)[5],",",
                                  predictors(rfe.T.net)[6],",", predictors(rfe.T.net)[7]),x = 6.5,
           y = max(rfe.T.net$results$RMSE+rfe.T.net$results$RMSESD)+50, size = 1.7, colour = "black")+ 
  ylim(100,450)+
  geom_errorbar(aes(ymin=(rfe.T.net$results$RMSE-rfe.T.net$results$RMSESD), ymax=(rfe.T.net$results$RMSE+rfe.T.net$results$RMSESD)), width=0.1,alpha=0.8)

rfe.T.pls_Plot <- ggplot(rfe.T.pls,aes(rfe.T.pls$variables,rfe.T.pls$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "TSS-PLS")+
  annotate("label", label = paste("Variable: ",predictors(rfe.T.pls)[1],",",
                                  predictors(rfe.T.pls)[2],",", predictors(rfe.T.pls)[3],",",
                                  predictors(rfe.T.pls)[4],",", predictors(rfe.T.pls)[5]),x = 6.5,
           y = max(rfe.T.pls$results$RMSE+rfe.T.pls$results$RMSESD)+50, size = 1.8, colour = "black")+ 
  ylim(100,450)+
  geom_errorbar(aes(ymin=(rfe.T.pls$results$RMSE-rfe.T.pls$results$RMSESD), ymax=(rfe.T.pls$results$RMSE+rfe.T.pls$results$RMSESD)), width=0.1,alpha=0.8)

rfe.T.rf_Plot <- ggplot(rfe.T.rf,aes(rfe.T.rf$variables,rfe.T.rf$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "TSS-RF")+
  annotate("label", label = paste("Variable: ",predictors(rfe.T.rf)[1],",",
                                  predictors(rfe.T.rf)[2],",", predictors(rfe.T.rf)[3],",",
                                  predictors(rfe.T.rf)[4],",", predictors(rfe.T.rf)[5]),x = 6.5,
           y = max(rfe.T.rf$results$RMSE+rfe.T.rf$results$RMSESD)+50, size = 2, colour = "black")+ 
  ylim(100,450)+
  geom_errorbar(aes(ymin=(rfe.T.rf$results$RMSE-rfe.T.rf$results$RMSESD), ymax=(rfe.T.rf$results$RMSE+rfe.T.rf$results$RMSESD)), width=0.1,alpha=0.8)

rfe.T.knn_Plot <- ggplot(rfe.T.knn,aes(rfe.T.knn$variables,rfe.T.knn$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "TSS-KNN")+
  annotate("label", label = paste("Variable: ",predictors(rfe.T.knn)[1],",",
                                  predictors(rfe.T.knn)[2],",", predictors(rfe.T.knn)[3],",",
                                  predictors(rfe.T.knn)[4],",", predictors(rfe.T.knn)[5]),x = 6.5,
           y = max(rfe.T.knn$results$RMSE+rfe.T.knn$results$RMSESD)+50, size = 2, colour = "black")+ 
  ylim(100,450)+
  geom_errorbar(aes(ymin=(rfe.T.knn$results$RMSE-rfe.T.knn$results$RMSESD), ymax=(rfe.T.knn$results$RMSE+rfe.T.knn$results$RMSESD)), width=0.1,alpha=0.8)

rfe.T.svm_Plot <- ggplot(rfe.T.svm,aes(rfe.T.svm$variables,rfe.T.svm$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "TSS-SVR")+
  annotate("label", label = paste("Variable: ",predictors(rfe.T.svm)[1],",",
                                  predictors(rfe.T.svm)[2],",", predictors(rfe.T.svm)[3],",",
                                  predictors(rfe.T.svm)[4],",", predictors(rfe.T.svm)[5],",",
                                  predictors(rfe.T.svm)[6],",", predictors(rfe.T.svm)[7]),x = 6.5,
           y = max(rfe.T.svm$results$RMSE+rfe.T.svm$results$RMSESD)+50, size = 1.8, colour = "black")+ 
  ylim(100,450)+
  geom_errorbar(aes(ymin=(rfe.T.svm$results$RMSE-rfe.T.svm$results$RMSESD), ymax=(rfe.T.svm$results$RMSE+rfe.T.svm$results$RMSESD)), width=0.1,alpha=0.8)

rfe.T.cub_Plot <- ggplot(rfe.T.cub,aes(rfe.T.cub$variables,rfe.T.cub$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "TSS-CUB")+
  annotate("label", label = paste("Variable: ",predictors(rfe.T.cub)[1],",",
                                  predictors(rfe.T.cub)[2],",", predictors(rfe.T.cub)[3],",",
                                  predictors(rfe.T.cub)[4],",", predictors(rfe.T.cub)[5],",",
                                  predictors(rfe.T.cub)[6],",", predictors(rfe.T.cub)[7]),x = 6.5,
           y = max(rfe.T.cub$results$RMSE+rfe.T.cub$results$RMSESD)+50, size = 1.8, colour = "black")+ 
  ylim(100,450)+
  geom_errorbar(aes(ymin=(rfe.T.cub$results$RMSE-rfe.T.cub$results$RMSESD), ymax=(rfe.T.cub$results$RMSE+rfe.T.cub$results$RMSESD)), width=0.1,alpha=0.8)

rfe.T_Plot <- ggarrange(rfe.T.net_Plot, rfe.T.pls_Plot, rfe.T.knn_Plot,rfe.T.svm_Plot, rfe.T.rf_Plot, rfe.T.cub_Plot,
                        labels = "AUTO", nrow=3, ncol = 2, align = "v")
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/FigureS3-2.jpg",plot=rfe.T_Plot,width = 6, height = 4, dpi = 600)


#### E.coli REF Models ####
rfe.E.lm <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_RFE/E_rfe_lm.rds")
#e.rfe.fun("lm",WQ.Cp)
rfe.E.pls <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_RFE/E_rfe_pls.rds")
#e.rfe.fun("pls",WQ.Cp,tuneLength = 5)
rfe.E.rf <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_RFE/E_rfe_rf.rds")
#e.rfe.fun("rf",WQ.Cp)
rfe.E.knn <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_RFE/E_rfe_knn.rds")
#e.rfe.fun("knn",WQ.Cp)
rfe.E.svm <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_RFE/E_rfe_svm.rds")
#e.rfe.fun("svmRadial",WQ.Cp)
rfe.E.cub <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_RFE/E_rfe_cub.rds")
#e.rfe.fun("cubist",WQ.Cp)
rfe.E.net <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_RFE/E_rfe_rnn.rds")
#e.rfe.fun("qrnn",WQ.Cp)

#saveRDS(rfe.E.net,"G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_RFE/E_rfe_rnn.rds")

###  E.coli REF Plot ### 
rfe.E.net_Plot <- ggplot(rfe.E.net,aes(rfe.E.net$variables,rfe.E.net$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "E. coli-QRNN")+
  annotate("label", label = paste("Variable: ",predictors(rfe.E.net)[1],",",
                                  predictors(rfe.E.net)[2],",", predictors(rfe.E.net)[3],",",
                                  predictors(rfe.E.net)[4],",", predictors(rfe.E.net)[5],",",
                                  predictors(rfe.E.net)[6]),x = 6.5,
           y = max(rfe.E.net$results$RMSE+rfe.E.net$results$RMSESD)+10^7.5, size = 1.8, colour = "black")+ 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(10^6,10^8),
                labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  geom_errorbar(aes(ymin=(rfe.E.net$results$RMSE-rfe.E.net$results$RMSESD), ymax=(rfe.E.net$results$RMSE+rfe.E.net$results$RMSESD)), width=0.1,alpha=0.8)

rfe.E.pls_Plot <- ggplot(rfe.E.pls,aes(rfe.E.pls$variables,rfe.E.pls$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "E. coli-PLS")+
  annotate("label", label = paste("Variable: ",predictors(rfe.E.pls)[1],",",
                                  predictors(rfe.E.pls)[2],",", predictors(rfe.E.pls)[3],",",
                                  predictors(rfe.E.pls)[4],",", predictors(rfe.E.pls)[5],",",
                                  predictors(rfe.E.pls)[6],",", predictors(rfe.E.pls)[7]),x = 6.5,
           y = max(rfe.E.pls$results$RMSE+rfe.E.pls$results$RMSESD)+10^7.5, size = 1.8, colour = "black")+ 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(10^6,10^8),
                labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  geom_errorbar(aes(ymin=(rfe.E.pls$results$RMSE-rfe.E.pls$results$RMSESD), ymax=(rfe.E.pls$results$RMSE+rfe.E.pls$results$RMSESD)), width=0.1,alpha=0.8)

rfe.E.rf_Plot <- ggplot(rfe.E.rf,aes(rfe.E.rf$variables,rfe.E.rf$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "E. coli-RF")+
  annotate("label", label = paste("Variable: ",predictors(rfe.E.rf)[1],",",
                                  predictors(rfe.E.rf)[2],",", predictors(rfe.E.rf)[3],",",
                                  predictors(rfe.E.rf)[4],",", predictors(rfe.E.rf)[5],",",
                                  predictors(rfe.E.rf)[6],",", predictors(rfe.E.rf)[7],",",
                                  predictors(rfe.E.rf)[8]),x = 6.5,
           y = max(rfe.E.rf$results$RMSE+rfe.E.rf$results$RMSESD)+10^7.5, size = 2, colour = "black")+ 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(10^6,10^8),
                labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  geom_errorbar(aes(ymin=(rfe.E.rf$results$RMSE-rfe.E.rf$results$RMSESD), ymax=(rfe.E.rf$results$RMSE+rfe.E.rf$results$RMSESD)), width=0.1,alpha=0.8)

rfe.E.knn_Plot <- ggplot(rfe.E.knn,aes(rfe.E.knn$variables,rfe.E.knn$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "E. coli-KNN")+
  annotate("label", label = paste("Variable: ",predictors(rfe.E.knn)[1],",",
                                  predictors(rfe.E.knn)[2],",", predictors(rfe.E.knn)[3],",",
                                  predictors(rfe.E.knn)[4],",", predictors(rfe.E.knn)[5],",",
                                  predictors(rfe.E.knn)[6],",", predictors(rfe.E.knn)[7],",",
                                  predictors(rfe.E.knn)[8]),x = 6.5,
           y = max(rfe.E.knn$results$RMSE+rfe.E.knn$results$RMSESD)+10^7.5, size = 2, colour = "black")+ 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(10^6,10^8),
                labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  geom_errorbar(aes(ymin=(rfe.E.knn$results$RMSE-rfe.E.knn$results$RMSESD), ymax=(rfe.E.knn$results$RMSE+rfe.E.knn$results$RMSESD)), width=0.1,alpha=0.8)

rfe.E.svm_Plot <- ggplot(rfe.E.svm,aes(rfe.E.svm$variables,rfe.E.svm$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "E. coli-SVR")+
  annotate("label", label = paste("Variable: ",predictors(rfe.E.svm)[1],",",
                                  predictors(rfe.E.svm)[2],",", predictors(rfe.E.svm)[3],",",
                                  predictors(rfe.E.svm)[4],",", predictors(rfe.E.svm)[5]),x = 6.5,
           y = max(rfe.E.svm$results$RMSE+rfe.E.svm$results$RMSESD)+10^7.5, size = 1.8, colour = "black")+ 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(10^6,10^8),
                labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  geom_errorbar(aes(ymin=(rfe.E.svm$results$RMSE-rfe.E.svm$results$RMSESD), ymax=(rfe.E.svm$results$RMSE+rfe.E.svm$results$RMSESD)), width=0.1,alpha=0.8)

rfe.E.cub_Plot <- ggplot(rfe.E.cub,aes(rfe.E.cub$variables,rfe.E.cub$results$RMSE))+geom_point()+theme_bw()+
  theme(text = element_text(size=8))+labs(y = "RMSE",title = "E. coli-CUB")+
  annotate("label", label = paste("Variable: ",predictors(rfe.E.cub)[1],",",
                                  predictors(rfe.E.cub)[2],",", predictors(rfe.E.cub)[3],",",
                                  predictors(rfe.E.cub)[5],",",
                                  predictors(rfe.E.cub)[6]),x = 6.5,
           y = max(rfe.E.cub$results$RMSE+rfe.E.cub$results$RMSESD)+10^7.5, size = 1.8, colour = "black")+ 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(10^6,10^8),
                labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  geom_errorbar(aes(ymin=(rfe.E.cub$results$RMSE-rfe.E.cub$results$RMSESD), ymax=(rfe.E.cub$results$RMSE+rfe.E.cub$results$RMSESD)), width=0.1,alpha=0.8)

rfe.E_Plot <- ggarrange(rfe.E.net_Plot, rfe.E.pls_Plot,rfe.E.knn_Plot,rfe.E.svm_Plot, rfe.E.rf_Plot,rfe.E.cub_Plot,
                        labels = "AUTO", nrow=3, ncol = 2, align = "v")

#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/FigureS3-3.jpg",plot=rfe.E_Plot,width = 6, height = 4, dpi = 600)

#### REF Plot ####
REF.R <- ggarrange(rfe.C.pls_Plot, rfe.C.svm_Plot, rfe.C.cub_Plot, rfe.C.net_Plot,
                      rfe.T.pls_Plot, rfe.T.svm_Plot, rfe.T.cub_Plot, rfe.T.net_Plot,
                      rfe.E.pls_Plot, rfe.E.svm_Plot, rfe.E.cub_Plot, rfe.E.net_Plot,
                      labels = "AUTO", ncol=4,nrow = 3, align = "v")
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/FigureS3.tiff",plot=REF.R,width = 10, height = 5, dpi = 600)


#### Test & Train Datasets (Separate 70/30) ####
set.seed(1234)
WQ2.index <- sample(1:nrow(WQ2),size=nrow(WQ2)*0.7,replace = FALSE)
WQ2.train <- WQ2[WQ2.index,]
WQ2.test <- WQ2[-WQ2.index,]


#### Cross-Validation Determination ####
k_values <- c(3, 5, 10, 15, 20)
models <- c("pls", "rf", "knn", "svmRadial", "cubist", "qrnn")
k_fold <- data.frame(k = numeric(), method = character(), accuracy = numeric(), std = numeric(), par = character())

for (i in 1:length(k_values)) {
  for (j in 1:length(models)) {
    set.seed(123)
    model <- train(COD ~ Group+Turb+NH4+NO3+Temp+Color+EC+pH, data = WQ.Cp, method = models[j], trControl = trainControl(method = "cv", number = k_values[i]))
    acc <- mean(model$resample$RMSE)
    accsd <- sd(model$resample$RMSE)
    k_fold <- rbind(k_fold, data.frame(k = k_values[i], method = models[j], accuracy = acc, std = accsd, par = "COD"))
  }
}

for (i in 1:length(k_values)) {
  for (j in 1:length(models)) {
    set.seed(123)
    model <- train(TSS ~ Group+Turb+NH4+NO3+Temp+Color+EC+pH, data = WQ.Cp, method = models[j], trControl = trainControl(method = "cv", number = k_values[i]))
    acc <- mean(model$resample$RMSE)
    accsd <- sd(model$resample$RMSE)
    k_fold <- rbind(k_fold, data.frame(k = k_values[i], method = models[j], accuracy = acc, std = accsd, par = "TSS"))
  }
}

for (i in 1:length(k_values)) {
  for (j in 1:length(models)) {
    set.seed(123)
    model <- train(E.coli ~ Group+Turb+NH4+NO3+Temp+Color+EC+pH, data = WQ.Cp, method = models[j], trControl = trainControl(method = "cv", number = k_values[i]))
    acc <- mean(model$resample$RMSE)
    accsd <- sd(model$resample$RMSE)
    k_fold <- rbind(k_fold, data.frame(k = k_values[i], method = models[j], accuracy = acc, std = accsd, par = "E.coli"))
  }
}

k_df <- k_fold
k_df$method <-  ifelse(k_df$method == c("pls", "rf", "knn", "svmRadial", "cubist", "qrnn"), c("PLS", "RF", "KNN", "SVR", "CUB", "QRNN"), 0)
k_df$method <-factor(k_df$method , levels=c("PLS", "RF", "KNN", "SVR", "CUB", "QRNN"))

#k_df <- read.csv(file = 'G:/My Drive/R project/GitHub/MLsensor/Raw_data/K-fold_result.csv',header = T, sep = ",")
k_df$method <- ifelse(k_df$method == "RNN", "QRNN", as.character(k_df$method))
k_df$method <- factor(k_df$method , levels=c("PLS", "RF", "KNN", "SVR", "CUB", "QRNN"))
k_df$par <- factor(k_df$par, levels=c("COD", "TSS", "E.coli"))

k_foldx <- k_df %>% group_by(par, k) %>%
  summarise(dv = mean(accuracy, na.rm=TRUE)) %>%
  mutate(diff=(abs(lag(dv,default=first(dv))-dv))/lag(dv,default=first(dv))*100) %>% ungroup() %>% group_by(k) %>% summarise(mdiff = mean(diff, na.rm=TRUE))
k_foldx

k_plot <- ggplot(subset(k_df, method %in% c("PLS", "SVR", "CUB", "QRNN")), aes(x = k, y = accuracy, group = par, color = method)) + 
  geom_line() + geom_point() + 
  #geom_errorbar(aes(ymin=(accuracy-std), ymax=(accuracy+std)), width=0.5,alpha=0.8)+
  facet_grid(par ~ method, scales="free") + theme_bw() +
  scale_color_manual(values = c(colx(4)), guide = guide_legend(title = NULL))+
  labs(x = "k", y = "RMSE", title = element_blank())+
  theme(text = element_text(size=12),legend.position="none",
        strip.background=element_blank(), #strip.text.y=element_text(size=8,face="bold.italic"), 
        axis.line=element_line(),panel.grid.major = element_blank())
k_plot
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/FigureS1.tiff",plot=k_plot,width = 6, height = 4, dpi = 600)


#### COD Predict ####
set.seed(1234)
tune.method <- trainControl(method = "repeatedcv", number=15, repeats=3, selectionFunction = "best")
c.mod.fun <- function(mod.,inputT,...){
  caret::train(COD ~ Group+Turb+NH4+NO3+Color+EC+pH
               , inputT, method = mod., trControl = tune.method, na.action = na.pass, preProc = c("center", "scale", "nzv"
               ),...)
}

set.seed(1234); WQ.C.lm <- #c.mod.fun("lm",WQ2.train)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_EST/COD_lm.rds")
set.seed(1234); WQ.C.pls <- #c.mod.fun("pls",WQ2.train,tuneLength = 5)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_EST/COD_pls.rds")
set.seed(1234); WQ.C.rf <- #c.mod.fun("rf",WQ2.train,ntree=25)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_EST/COD_rf.rds")
set.seed(1234); WQ.C.knn <- #c.mod.fun("knn",WQ2.train,tuneGrid = expand.grid (k = seq(from = 3, to = 10, by = 1)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_EST/COD_knn.rds")
set.seed(1234); WQ.C.svm <- #c.mod.fun("svmRadial",WQ2.train,tuneGrid=expand.grid(.sigma=seq(from=0.01,to=0.1,by=0.01),.C=seq(from = 5,to = 10, by = 0.5)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_EST/COD_svm.rds")
set.seed(1234); WQ.C.cub <- #c.mod.fun("cubist",WQ2.train,tuneGrid=expand.grid(.committees=seq(from=3,to=10,by=1),.neighbors = seq(from = 3, to = 9, by = 1)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_EST/COD_cub.rds")
set.seed(1234); WQ.C.qrnn <- #c.mod.fun("qrnn",WQ2.train,tuneGrid=expand.grid(.n.hidden=seq(from=2,to=5,by=1),.penalty = 10^(seq(from = -2, to = 2, by = 0.25)),.bag=FALSE))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_EST/COD_qrnn.rds")
#set.seed(1234); WQ.C.xg <- c.mod.fun("xgbTree",WQ2.train)

#saveRDS(WQ.C.qrnn,"G:/My Drive/R project/GitHub/MLsensor/Save_model/COD_EST/COD_qrnn.rds")

###  COD Train Variables Importance ### 
C.lm.Imp <- ggplot(varImp(WQ.C.lm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-LM")
C.pls.Imp <-ggplot(varImp(WQ.C.pls, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-PLS")
C.rf.Imp <- ggplot(varImp(WQ.C.rf, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-RF")
C.knn.Imp <- ggplot(varImp(WQ.C.knn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-KNN")
C.svm.Imp <- ggplot(varImp(WQ.C.svm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-SVR")
C.cub.Imp <- ggplot(varImp(WQ.C.cub, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-CUB")
C.qrnn.Imp <- ggplot(varImp(WQ.C.qrnn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-QRNN")
C.varImp <- ggarrange(C.pls.Imp, C.svm.Imp, C.cub.Imp,C.qrnn.Imp, ncol=2, nrow = 2, align = "v", labels = "AUTO", common.legend = TRUE, legend="bottom")
C.varImp
#ggsave("E:/Dropbox/R project/Machine Learning/Figure/C_varImp.jpg",plot=C.varImp,width = 6, height = 4, dpi = 600)

### COD Train RMSE Boxplot ### 
C.results <- resamples(list("qrnn" = WQ.C.qrnn, "knn" = WQ.C.knn, "svr" = WQ.C.svm,"pls" = WQ.C.pls, "rf" = WQ.C.rf, "CUB" = WQ.C.cub))
C.RMSE <- data.frame(C.results$values$`qrnn~RMSE`,C.results$values$`knn~RMSE`,C.results$values$`svr~RMSE`,
                     C.results$values$`pls~RMSE`,C.results$values$`rf~RMSE`,C.results$values$`CUB~RMSE`)
colnames(C.RMSE) <- c("QRNN","KNN","SVR","PLS","RF","CUB")
C.RMSE <- C.RMSE %>% gather(key = "Model", value = "RMSE", QRNN, KNN, SVR, PLS, RF, CUB) %>% convert_as_factor(Model)
C.Ttest <- C.RMSE %>% pairwise_t_test(RMSE ~ Model, paired = TRUE, alt = c("two.sided"), conf.level = 0.99, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Model")
C.RMSE.p <- ggboxplot(C.RMSE, x = "Model", y = "RMSE",color = "black", legend = "none", order = c("PLS","KNN","SVR","RF","QRNN","CUB"),outlier.shape = NA) +
  geom_point(data = C.RMSE, aes(x = Model, y = RMSE, colour = Model, alpha = 0.9), position = position_jitter(width = 0.2))+
  geom_hline(yintercept = mean(C.RMSE$RMSE), linetype = 2) + 
  coord_flex_cart(bottom=brackets_horizontal(), left=capped_vertical('none'))+
  theme_bw()+labs(x = element_blank(), y = "RMSE (COD mg/L)")+
  theme(text = element_text(size=15),legend.title=element_blank(),legend.position = "none",
        panel.border=element_blank(), axis.line = element_line(),legend.background = element_rect(colour='grey')) +
  scale_color_manual(values = c("PLS"=colx(6)[1],"KNN"=colx(6)[2],"SVR"=colx(6)[3],"RF"=colx(6)[4],"QRNN"=colx(6)[5],"CUB"=colx(6)[6]),guide=guide_legend(title = "Algorithms", order = 1))
C.RMSE.p

### COD Test RMSE Boxplot ### 
cod.theme <- function(titl.,lab.data,...){
  list(geom_point(aes(fill=dataset), size=7.25, shape=21, color="black", stroke=0.5, alpha=0.75),
       theme_bw(),ylim(0,6000),xlim(0,6000),
       scale_fill_manual(values = c(colx2(2)), guide = guide_legend(title = "Dataset", title.position = "left",override.aes = list(shape=21,size=2.5), order = 1, nrow = 1)),
       theme(text = element_text(size=8.5),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = c(0.15, 0.85)),
       labs(x = "Predicted COD (mg/L)", y = "Measured COD (mg/L)",title = titl.),geom_abline(color = "black", linetype = 2, size = 1, alpha = 0.5),
       geom_richtext(data = lab.data,aes(3400, 900, label = label[1]),hjust = 0,size = 2.75,fill = "white", label.color = "black"),
       ...)
}

### LM ###
C.lm.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$COD, abs(predict.train(WQ.C.lm, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
C.lm.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$COD, abs(predict.train(WQ.C.lm, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),C.lm.df)
C.lm.df <- C.lm.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
C.lm.sr <- lm(X2 ~ X3, data=C.lm.df)
C.lm.df <- cbind(C.lm.df, rstandard(C.lm.sr))
colnames(C.lm.df) <- c("Sample","actual","predicted","dataset","residual")
C.lm.df$dataset <-factor(C.lm.df$dataset  , levels=c("Train","Test"))
C.lm.df$ML <- "LM"
C.lm.L <- C.lm.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
C.lm.L$label[1] <- paste("Train RMSE =", round(C.lm.L$RMSE[2], 0), "<br> Test RMSE =", round(C.lm.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(C.lm.L$R2[1], 2), "<br> MAPE =", round(C.lm.L$MAPE[1], 1))
C.lm.p <- ggplot(C.lm.df, aes(predicted, actual))+cod.theme(NULL,C.lm.L)+theme(legend.position="none")
C.lm.p

### QRNN ###
C.qrnn.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$COD, abs(predict.train(WQ.C.qrnn, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
C.qrnn.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$COD, abs(predict.train(WQ.C.qrnn, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),C.qrnn.df)
C.qrnn.df <- C.qrnn.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
C.qrnn.sr <- lm(X2 ~ X3, data=C.qrnn.df)
C.qrnn.df <- cbind(C.qrnn.df, rstandard(C.qrnn.sr))
colnames(C.qrnn.df) <- c("Sample","actual","predicted","dataset","residual")
C.qrnn.df$dataset <-factor(C.qrnn.df$dataset  , levels=c("Train","Test"))
C.qrnn.df$ML <- "QRNN"
C.qrnn.L <- C.qrnn.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
C.qrnn.L$label[1] <- paste("Train RMSE =", round(C.qrnn.L$RMSE[2], 0), "<br> Test RMSE =", round(C.qrnn.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(C.qrnn.L$R2[1], 2), "<br> MAPE =", round(C.qrnn.L$MAPE[1], 1))
C.qrnn.p <- ggplot(C.qrnn.df, aes(predicted, actual))+cod.theme(NULL,C.qrnn.L)+theme(legend.position="none")
C.qrnn.p

### PLS ###
C.pls.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$COD, abs(predict.train(WQ.C.pls, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
C.pls.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$COD, abs(predict.train(WQ.C.pls, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),C.pls.df)
C.pls.df <- C.pls.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
C.pls.sr <- lm(X2 ~ X3, data=C.pls.df)
C.pls.df <- cbind(C.pls.df, rstandard(C.pls.sr))
colnames(C.pls.df) <- c("Sample","actual","predicted","dataset","residual")
C.pls.df$dataset <-factor(C.pls.df$dataset  , levels=c("Train","Test"))
C.pls.df$ML <- "PLS"
C.pls.L <- C.pls.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
C.pls.L$label[1] <- paste("Train RMSE =", round(C.pls.L$RMSE[2], 0), "<br> Test RMSE =", round(C.pls.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(C.pls.L$R2[1], 2), "<br> MAPE =", round(C.pls.L$MAPE[1], 1))
C.pls.p <- ggplot(C.pls.df, aes(predicted, actual))+cod.theme(NULL,C.pls.L)
C.pls.p

### KNN ###
C.knn.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$COD, abs(predict.train(WQ.C.knn, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
C.knn.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$COD, abs(predict.train(WQ.C.knn, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),C.knn.df)
C.knn.df <- C.knn.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
C.knn.sr <- lm(X2 ~ X3, data=C.knn.df)
C.knn.df <- cbind(C.knn.df, rstandard(C.knn.sr))
colnames(C.knn.df) <- c("Sample","actual","predicted","dataset","residual")
C.knn.df$dataset <-factor(C.knn.df$dataset  , levels=c("Train","Test"))
C.knn.df$ML <- "KNN"
C.knn.L <- C.knn.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
C.knn.L$label[1] <- paste("Train RMSE =", round(C.knn.L$RMSE[2], 0), "<br> Test RMSE =", round(C.knn.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(C.knn.L$R2[1], 2), "<br> MAPE =", round(C.knn.L$MAPE[1], 1))
C.knn.p <- ggplot(C.knn.df, aes(predicted, actual))+cod.theme(NULL,C.knn.L)+theme(legend.position="none")
C.knn.p

### SVM ###
C.svm.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$COD, abs(predict.train(WQ.C.svm, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
C.svm.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$COD, abs(predict.train(WQ.C.svm, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),C.svm.df)
C.svm.df <- C.svm.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
C.svm.sr <- lm(X2 ~ X3, data=C.svm.df)
C.svm.df <- cbind(C.svm.df, rstandard(C.svm.sr))
colnames(C.svm.df) <- c("Sample","actual","predicted","dataset","residual")
C.svm.df$dataset <-factor(C.svm.df$dataset  , levels=c("Train","Test"))
C.svm.df$ML <- "SVR"
C.svm.L <- C.svm.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
C.svm.L$label[1] <- paste("Train RMSE =", round(C.svm.L$RMSE[2], 0), "<br> Test RMSE =", round(C.svm.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(C.svm.L$R2[1], 2), "<br> MAPE =", round(C.svm.L$MAPE[1], 1))
C.svm.p <- ggplot(C.svm.df, aes(predicted, actual))+cod.theme(NULL,C.svm.L)+theme(legend.position="none")
C.svm.p

### RF ###
C.rf.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$COD, abs(predict.train(WQ.C.rf, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
C.rf.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$COD, abs(predict.train(WQ.C.rf, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),C.rf.df)
C.rf.df <- C.rf.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
C.rf.sr <- lm(X2 ~ X3, data=C.rf.df)
C.rf.df <- cbind(C.rf.df, rstandard(C.rf.sr))
colnames(C.rf.df) <- c("Sample","actual","predicted","dataset","residual")
C.rf.df$dataset <-factor(C.rf.df$dataset  , levels=c("Train","Test"))
C.rf.df$ML <- "RF"
C.rf.L <- C.rf.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
C.rf.L$label[1] <- paste("Train RMSE =", round(C.rf.L$RMSE[2], 0), "<br> Test RMSE =", round(C.rf.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(C.rf.L$R2[1], 2), "<br> MAPE =", round(C.rf.L$MAPE[1], 1))
C.rf.p <- ggplot(C.rf.df, aes(predicted, actual))+cod.theme(NULL,C.rf.L)+theme(legend.position="none")
C.rf.p

### CUB ###
C.cub.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$COD, abs(predict.train(WQ.C.cub, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
C.cub.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$COD, abs(predict.train(WQ.C.cub, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),C.cub.df)
C.cub.df <- C.cub.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
C.cub.sr <- lm(X2 ~ X3, data=C.cub.df)
C.cub.df <- cbind(C.cub.df, rstandard(C.cub.sr))
colnames(C.cub.df) <- c("Sample","actual","predicted","dataset","residual")
C.cub.df$dataset <-factor(C.cub.df$dataset  , levels=c("Train","Test"))
C.cub.df$ML <- "CUB"
C.cub.L <- C.cub.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
C.cub.L$label[1] <- paste("Train RMSE =", round(C.cub.L$RMSE[2], 0), "<br> Test RMSE =", round(C.cub.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(C.cub.L$R2[1], 2), "<br> MAPE =", round(C.cub.L$MAPE[1], 1))
C.cub.p <- ggplot(C.cub.df, aes(predicted, actual))+cod.theme(NULL,C.cub.L)+theme(legend.position="none")
C.cub.p

### Save COD Result Scatter Plot ###
COD.R.p <- ggarrange(C.pls.p,C.svm.p,C.cub.p,C.qrnn.p,labels = "AUTO", ncol=4,nrow = 1, align = "v")
ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/Figure5-1.jpg",plot=COD.R.p,width = 12, height = 2, dpi = 600)

### COD Overall Scatter Plot ###
COD.df <- rbind(C.qrnn.df,C.pls.df,C.svm.df,C.cub.df)
COD.df$ML <-factor(COD.df$ML , levels=c("PLS","SVR","CUB","QRNN"))
COD.df.L <- rbind(C.pls.L,C.rf.L,C.cub.L,C.qrnn.L)
COD.df.L$ML <- c("PLS","PLS","SVR","SVR","CUB","CUB","QRNN","QRNN")
COD.df.L$label <- gsub("<br>", ";", COD.df.L$label)
COD.df.p <- ggplot(subset(COD.df, dataset %in% "Test"),aes(predicted, actual, fill = ML)) + geom_point(size=2.5,stroke=0.5,shape = 21,alpha = 0.8) + theme_bw() + 
  geom_abline(color = "black", linetype = 2, size = 1, alpha = 0.6)+
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]), 
                    guide = guide_legend(title = "ML Model", title.position = "top",override.aes = list(shape =  c(21,21,21,21)), order = 1))+
  scale_shape_manual(values = c(21,21,21,21),guide = F)+
  #scale_linetype_manual(values = c("longdash","longdash","longdash","longdash","longdash","longdash"), guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  #scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(1,10^3.5),
  #                   labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(1,10^3.5),
  #                   labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  labs(x = "Predicted COD (mg/L)", y = "Measured COD (mg/L)")+
  theme(text = element_text(size=15),legend.position = c(.35, 0.9), legend.direction='horizontal',
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        legend.text=element_text(size=10),legend.title=element_text(size=10))
COD.df.p

### COD Result Standardized Residual Plot ###
C.sr.theme <- function(titl.,...){
  list(geom_point(aes(fill = Sample),shape=21,size=2,color="black", stroke = 0.5),theme_bw(),ylim(-3,3),xlim(-200,5000),
       scale_alpha_continuous(limits = c(0,3), breaks = c(0,1,2,3),
                              guide = guide_legend(title = "Standardized residual", title.position = "top", order = 2)),
       scale_fill_manual(values = c(colx2(5)), 
                         guide = guide_legend(title = "Sampling location", title.position = "left",override.aes = list(shape = 21), order = 1)),
       theme(text = element_text(size=7),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = c(.9, 0.3)),
       labs(title = titl. ,x = "Predicted COD (mg/L)", y = "Standardized Residual"))
}

C.qrnn.p2 <- ggplot(data = C.qrnn.df, aes(predicted, residual)) + C.sr.theme("COD: QRNN") + guides(alpha=FALSE)
C.pls.p2 <- ggplot(data = C.pls.df, aes(predicted, residual)) + C.sr.theme("COD: PLS") + guides(alpha=FALSE)
C.knn.p2 <- ggplot(data = C.knn.df, aes(predicted, residual)) + C.sr.theme("COD: KNN") + guides(alpha=FALSE)
C.rf.p2 <- ggplot(data = C.rf.df, aes(predicted, residual)) + C.sr.theme("COD: RF") + guides(alpha=FALSE)
C.svm.p2 <- ggplot(data = C.svm.df, aes(predicted, residual)) + C.sr.theme("COD: SVR") + guides(alpha=FALSE)
C.cub.p2 <- ggplot(data = C.cub.df, aes(predicted, residual)) + C.sr.theme("COD: CUB") + guides(alpha=FALSE)
COD.SR.p <- ggarrange(C.pls.p2, C.svm.p2, C.cub.p2, C.qrnn.p2, labels = "AUTO", ncol=4,nrow = 1, align = "v",  common.legend = TRUE,legend = "bottom")
COD.SR.p

COD.SR.p2 <- ggplot(subset(COD.df, dataset %in% "Test"),aes(actual, residual, fill = ML)) + geom_point(size=2.5,stroke=0.5,shape = 21,alpha = 0.8) + theme_bw() + 
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]), 
                    guide = guide_legend(title = "ML Model", title.position = "top",override.aes = list(shape = 21), order = 1))+
  #scale_shape_manual(values = c(0,2,3,4,5), guide = guide_legend(title = "Sampling Location", title.position = "top", order = 2))+
  labs(x = "Predicted COD (mg/L)", y = "Standardized Residual")+ylim(-6,6)+
  theme(text = element_text(size=15),legend.position = c(.35, 0.9), legend.direction='horizontal',
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        legend.text=element_text(size=8),legend.title=element_text(size=8))

Cyplot <- ggplot(subset(COD.df, dataset %in% "Test"), aes(x=residual, color = ML, fill = ML)) + geom_histogram(alpha = 0.5, position="dodge", binwidth = 0.5) + labs(y = "count")  + xlim(-6,6)+
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]),guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  scale_color_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]),guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  coord_flip() + theme_bw() + theme(text = element_text(size=15),legend.position="none",axis.title.y=element_blank(),
                                    axis.text.y=element_blank(),axis.ticks.y=element_blank(), plot.margin = unit(c(0.38,1,0.38,-0.5), "lines"))
Csr.p <- ggarrange(COD.SR.p2, Cyplot, ncol = 2, nrow = 1, widths = c(3.5, 1), common.legend = F)
Csr.p


#### TSS Predict ####
set.seed(1234)
tune.method <- trainControl(method = "repeatedcv", number=15, repeats=3, selectionFunction = "best")
t.mod.fun <- function(mod.,inputT,...){
  caret::train(TSS ~ Group+Turb+NH4+NO3+Color+EC+pH
               , inputT, method = mod., trControl = tune.method, na.action = na.pass, preProc = c("center", "scale", "nzv"
               ),...)
}

set.seed(1234); WQ.T.lm <- #t.mod.fun("lm",WQ2.train)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_EST/TSS_lm.rds")
set.seed(1234); WQ.T.pls <- #t.mod.fun("pls",WQ2.train,tuneLength = 5)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_EST/TSS_pls.rds")
set.seed(1234); WQ.T.rf <- #t.mod.fun("rf",WQ2.train,ntree=25)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_EST/TSS_rf.rds")
set.seed(1234); WQ.T.knn <- #t.mod.fun("knn",WQ2.train,tuneGrid = expand.grid (k = seq(from = 3, to = 10, by = 1)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_EST/TSS_knn.rds")
set.seed(1234); WQ.T.svm <- #t.mod.fun("svmRadial",WQ2.train,tuneGrid=expand.grid(.sigma=seq(from=0.01,to=0.05,by=0.01),.C=seq(from=10,to=20,by=0.5)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_EST/TSS_svm.rds")
set.seed(1234); WQ.T.cub <- #t.mod.fun("cubist",WQ2.train,tuneGrid=expand.grid(.committees=seq(from=3,to=10,by = 1),.neighbors=seq(from=3,to=9,by=1)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_EST/TSS_cub.rds")
set.seed(1234); WQ.T.qrnn <- #t.mod.fun("qrnn",WQ2.train,tuneGrid=expand.grid(.n.hidden=seq(from=2,to=5,by=1),.penalty = 10^(seq(from = -2, to = 2, by = 0.25)),.bag=FALSE))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_EST/TSS_qrnn.rds")

#saveRDS(WQ.T.qrnn,"G:/My Drive/R project/GitHub/MLsensor/Save_model/TSS_EST/TSS_qrnn.rds")

###  TSS Train Variables Importance ### 
T.lm.Imp <- ggplot2::ggplot(varImp(WQ.T.lm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "TSS-LM")
T.pls.Imp <- ggplot2::ggplot(varImp(WQ.T.pls, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "TSS-PLS")
T.rf.Imp <- ggplot2::ggplot(varImp(WQ.T.rf, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "TSS-RF")
T.knn.Imp <- ggplot2::ggplot(varImp(WQ.T.knn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "TSS-KNN")
T.svm.Imp <- ggplot2::ggplot(varImp(WQ.T.svm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "TSS-SVR")
T.cub.Imp <- ggplot2::ggplot(varImp(WQ.T.cub, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "TSS-CUB")
T.qrnn.Imp <- ggplot2::ggplot(varImp(WQ.T.qrnn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "TSS-QRNN")
T.varImp <- ggarrange(T.lm.Imp, T.pls.Imp, T.knn.Imp, T.svm.Imp, T.rf.Imp, T.qrnn.Imp,ncol=3, nrow = 2, align = "v", labels = "AUTO", common.legend = TRUE, legend="bottom")
T.varImp
#ggsave("E:/Dropbox/R project/Machine Learning/Figure/T_varImp.jpg",plot=T.varImp,width = 6, height = 4, dpi = 600)

### TSS Train RMSE Boxplot ### 
T.results <- resamples(list("qrnn" = WQ.T.qrnn, "knn" = WQ.T.knn, "svr" = WQ.T.svm,"pls" = WQ.T.pls, "rf" = WQ.T.rf, "CUB" = WQ.T.cub))
T.RMSE <- data.frame(T.results$values$`qrnn~RMSE`,T.results$values$`knn~RMSE`,T.results$values$`svr~RMSE`,
                     T.results$values$`pls~RMSE`,T.results$values$`rf~RMSE`,T.results$values$`CUB~RMSE`)
colnames(T.RMSE) <- c("QRNN","KNN","SVR","PLS","RF","CUB")
T.RMSE <- T.RMSE %>% gather(key = "Model", value = "RMSE", QRNN, KNN, SVR, PLS, RF, CUB) %>% convert_as_factor(Model)
T.Ttest <- T.RMSE %>% pairwise_t_test(RMSE ~ Model, paired = TRUE, alt = c("two.sided"), conf.level = 0.99, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Model")
T.RMSE.p <- ggboxplot(T.RMSE, x = "Model", y = "RMSE",color = "black", legend = "none", order = c("PLS","KNN","SVR","RF","QRNN","CUB"),outlier.shape = NA) +
  geom_point(data = T.RMSE, aes(x = Model, y = RMSE, colour = Model, alpha = 0.9), position = position_jitter(width = 0.2))+
  geom_hline(yintercept = mean(T.RMSE$RMSE), linetype = 2) + 
  coord_flex_cart(bottom=brackets_horizontal(), left=capped_vertical('none'))+
  theme_bw()+labs(x = element_blank(), y = "RMSE (TSS mg/L)")+
  theme(text = element_text(size=15),legend.title=element_blank(),legend.position = "none",
        panel.border=element_blank(), axis.line = element_line(),legend.background = element_rect(colour='grey')) +
  scale_color_manual(values = c("PLS"=colx(6)[1],"KNN"=colx(6)[2],"SVR"=colx(6)[3],"RF"=colx(6)[4],"QRNN"=colx(6)[5],"CUB"=colx(6)[6]),guide=guide_legend(title = "Algorithms", order = 1))
T.RMSE.p

### TSS Test RMSE Boxplot ### 
tss.theme <- function(titl.,lab.data,...){
  list(geom_point(aes(fill=dataset), size=7.25, shape=21, color="black", stroke=0.5, alpha=0.75),
       theme_bw(),ylim(-50,1250),xlim(-50,1250),
       scale_fill_manual(values = c(colx2(2)), guide = guide_legend(title = NULL, title.position = "top",override.aes = list(shape=21,size=2.5), order = 1)),
       theme(text = element_text(size=8.5),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = c(0.15, 0.85)),
       labs(x = "Predicted TSS (mg/L)", y = "Measured TSS (mg/L)",title = titl.),geom_abline(color = "black", linetype = 2, size = 1, alpha = 0.5),
       geom_richtext(data = lab.data,aes(675, 150, label = label[1]),hjust = 0,size = 2.75,fill = "white", label.color = "black"),
       ...)
}

### LM ###
T.lm.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.T.lm, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
T.lm.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.T.lm, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),T.lm.df)
T.lm.df <- T.lm.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
T.lm.sr <- lm(X2 ~ X3, data=T.lm.df)
T.lm.df <- cbind(T.lm.df, rstandard(T.lm.sr))
colnames(T.lm.df) <- c("Sample","actual","predicted","dataset","residual")
T.lm.df$dataset <-factor(T.lm.df$dataset  , levels=c("Train","Test"))
T.lm.df$ML <- "LM"
T.lm.L <- T.lm.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
T.lm.L$label[1] <- paste("Train RMSE =", round(T.lm.L$RMSE[2], 0), "<br> Test RMSE =", round(T.lm.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(T.lm.L$R2[1], 2), "<br> MAPE =", round(T.lm.L$MAPE[1], 1))
T.lm.p <- ggplot(T.lm.df, aes(predicted, actual))+tss.theme(NULL,T.lm.L)+theme(legend.position="none")
T.lm.p

### QRNN ###
T.qrnn.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.T.qrnn, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
T.qrnn.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.T.qrnn, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),T.qrnn.df)
T.qrnn.df <- T.qrnn.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
T.qrnn.sr <- lm(X2 ~ X3, data=T.qrnn.df)
T.qrnn.df <- cbind(T.qrnn.df, rstandard(T.qrnn.sr))
colnames(T.qrnn.df) <- c("Sample","actual","predicted","dataset","residual")
T.qrnn.df$dataset <-factor(T.qrnn.df$dataset  , levels=c("Train","Test"))
T.qrnn.df$ML <- "QRNN"
T.qrnn.L <- T.qrnn.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
T.qrnn.L$label[1] <- paste("Train RMSE =", round(T.qrnn.L$RMSE[2], 0), "<br> Test RMSE =", round(T.qrnn.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(T.qrnn.L$R2[1], 2), "<br> MAPE =", round(T.qrnn.L$MAPE[1], 1))
T.qrnn.p <- ggplot(T.qrnn.df, aes(predicted, actual))+tss.theme(NULL,T.qrnn.L)+theme(legend.position="none")
T.qrnn.p

### PLS ###
T.pls.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.T.pls, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
T.pls.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.T.pls, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),T.pls.df)
T.pls.df <- T.pls.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
T.pls.sr <- lm(X2 ~ X3, data=T.pls.df)
T.pls.df <- cbind(T.pls.df, rstandard(T.pls.sr))
colnames(T.pls.df) <- c("Sample","actual","predicted","dataset","residual")
T.pls.df$dataset <-factor(T.pls.df$dataset  , levels=c("Train","Test"))
T.pls.df$ML <- "PLS"
T.pls.L <- T.pls.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
T.pls.L$label[1] <- paste("Train RMSE =", round(T.pls.L$RMSE[2], 0), "<br> Test RMSE =", round(T.pls.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(T.pls.L$R2[1], 2), "<br> MAPE =", round(T.pls.L$MAPE[1], 1))
T.pls.p <- ggplot(T.pls.df, aes(predicted, actual))+tss.theme(NULL,T.pls.L)+theme(legend.position="none")
T.pls.p

### KNN ###
T.knn.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.T.knn, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
T.knn.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.T.knn, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),T.knn.df)
T.knn.df <- T.knn.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
T.knn.sr <- lm(X2 ~ X3, data=T.knn.df)
T.knn.df <- cbind(T.knn.df, rstandard(T.knn.sr))
colnames(T.knn.df) <- c("Sample","actual","predicted","dataset","residual")
T.knn.df$dataset <-factor(T.knn.df$dataset  , levels=c("Train","Test"))
T.knn.df$ML <- "KNN"
T.knn.L <- T.knn.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
T.knn.L$label[1] <- paste("Train RMSE =", round(T.knn.L$RMSE[2], 0), "<br> Test RMSE =", round(T.knn.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(T.knn.L$R2[1], 2), "<br> MAPE =", round(T.knn.L$MAPE[1], 1))
T.knn.p <- ggplot(T.knn.df, aes(predicted, actual))+tss.theme(NULL,T.knn.L)+theme(legend.position="none")
T.knn.p

### SVM ###
T.svm.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.T.svm, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
T.svm.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.T.svm, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),T.svm.df)
T.svm.df <- T.svm.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
T.svm.sr <- lm(X2 ~ X3, data=T.svm.df)
T.svm.df <- cbind(T.svm.df, rstandard(T.svm.sr))
colnames(T.svm.df) <- c("Sample","actual","predicted","dataset","residual")
T.svm.df$dataset <-factor(T.svm.df$dataset  , levels=c("Train","Test"))
T.svm.df$ML <- "SVR"
T.svm.L <- T.svm.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
T.svm.L$label[1] <- paste("Train RMSE =", round(T.svm.L$RMSE[2], 0), "<br> Test RMSE =", round(T.svm.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(T.svm.L$R2[1], 2), "<br> MAPE =", round(T.svm.L$MAPE[1], 1))
T.svm.p <- ggplot(T.svm.df, aes(predicted, actual))+tss.theme(NULL,T.svm.L)+theme(legend.position="none")
T.svm.p

### RF ###
T.rf.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.T.rf, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
T.rf.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.T.rf, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),T.rf.df)
T.rf.df <- T.rf.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
T.rf.sr <- lm(X2 ~ X3, data=T.rf.df)
T.rf.df <- cbind(T.rf.df, rstandard(T.rf.sr))
colnames(T.rf.df) <- c("Sample","actual","predicted","dataset","residual")
T.rf.df$dataset <-factor(T.rf.df$dataset  , levels=c("Train","Test"))
T.rf.df$ML <- "RF"
T.rf.L <- T.rf.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
T.rf.L$label[1] <- paste("Train RMSE =", round(T.rf.L$RMSE[2], 0), "<br> Test RMSE =", round(T.rf.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(T.rf.L$R2[1], 2), "<br> MAPE =", round(T.rf.L$MAPE[1], 1))
T.rf.p <- ggplot(T.rf.df, aes(predicted, actual))+tss.theme(NULL,T.rf.L)+theme(legend.position="none")
T.rf.p

### CUB ###
T.cub.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.T.cub, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
T.cub.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.T.cub, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),T.cub.df)
T.cub.df <- T.cub.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
T.cub.sr <- lm(X2 ~ X3, data=T.cub.df)
T.cub.df <- cbind(T.cub.df, rstandard(T.cub.sr))
colnames(T.cub.df) <- c("Sample","actual","predicted","dataset","residual")
T.cub.df$dataset <-factor(T.cub.df$dataset  , levels=c("Train","Test"))
T.cub.df$ML <- "CUB"
T.cub.L <- T.cub.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
T.cub.L$label[1] <- paste("Train RMSE =", round(T.cub.L$RMSE[2], 0), "<br> Test RMSE =", round(T.cub.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(T.cub.L$R2[1], 2), "<br> MAPE =", 24.8)
T.cub.p <- ggplot(T.cub.df, aes(predicted, actual))+tss.theme(NULL,T.cub.L)+theme(legend.position="none")
T.cub.p

### Save TSS Result Scatter Plot ###
TSS.R.p <- ggarrange(T.pls.p, T.svm.p, T.cub.p, T.qrnn.p, labels = "AUTO", ncol=4,nrow = 1, align = "v")
ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/Figure5-2.jpg",plot=TSS.R.p,width = 12, height = 2, dpi = 600)

### TSS Overall Scatter Plot ###
TSS.df <- rbind(T.qrnn.df,T.pls.df,T.svm.df,T.cub.df)
TSS.df$ML <-factor(TSS.df$ML , levels=c("PLS","SVR","CUB","QRNN"))
TSS.df.L <- rbind(T.pls.L,T.svm.L,T.cub.L,T.qrnn.L)
TSS.df.L$ML <- c("PLS","PLS","SVR","SVR","CUB","CUB","QRNN","QRNN")
TSS.df.L$label <- gsub("<br>", ";", TSS.df.L$label)

TSS.df.p <- ggplot(subset(TSS.df, dataset %in% "Test"),aes(predicted, actual, fill = ML)) + geom_point(size=2.5,stroke=0.5,shape = 21,alpha = 0.8) + theme_bw() + 
  geom_abline(color = "black", linetype = 2, size = 1, alpha = 0.4)+
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]), 
                    guide = guide_legend(title = "ML Model", title.position = "top",override.aes = list(shape = 21), order = 1))+
  scale_shape_manual(values = c(16,16,16,16,16), guide = guide_legend(title = "Sampling Location", title.position = "top", order = 2))+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(1,10^3),
                labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(1,10^3),
                labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  labs(x = "Predicted TSS (mg/L)", y = "Measured TSS (mg/L)")+
  theme(text = element_text(size=15),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = "none", legend.box = "horizontal")
TSS.df.p

### TSS Result Standardized Residual Plot ###
t.sr.theme <- function(titl.,...){
  list(geom_point(aes(fill = Sample),shape=21,size=2,color="black", stroke = 0.5),theme_bw(),ylim(-3,3),xlim(-200,1250),
       scale_alpha_continuous(limits = c(0,3), breaks = c(0,1,2,3), guide = guide_legend(title = "Standardized residual", title.position = "top", order = 2)),
       scale_fill_manual(values = c(colx2(5)), 
                         guide = guide_legend(title = "Sampling location", title.position = "left",override.aes = list(shape = 21), order = 1)),
       theme(text = element_text(size=7),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = c(.9, 0.3)),
       labs(title = titl. ,x = "Predicted TSS (mg/L)", y = "Standardized Residual"))
}

T.qrnn.p2 <- ggplot(data = T.qrnn.df, aes(predicted, residual)) + t.sr.theme("TSS: QRNN") + guides(alpha=FALSE)
T.pls.p2 <- ggplot(data = T.pls.df, aes(predicted, residual)) + t.sr.theme("TSS: PLS") + guides(alpha=FALSE)
T.knn.p2 <- ggplot(data = T.knn.df, aes(predicted, residual)) + t.sr.theme("TSS: KNN") + guides(alpha=FALSE)
T.rf.p2 <- ggplot(data = T.rf.df, aes(predicted, residual)) + t.sr.theme("TSS: RF") + guides(alpha=FALSE)
T.svm.p2 <- ggplot(data = T.svm.df, aes(predicted, residual)) + t.sr.theme("TSS: SVR") + guides(alpha=FALSE)
T.cub.p2 <- ggplot(data = T.cub.df, aes(predicted, residual)) + t.sr.theme("TSS: CUB") + guides(alpha=FALSE)
TSS.SR.p <- ggarrange(T.pls.p2, T.svm.p2, T.cub.p2, T.qrnn.p2, labels = "AUTO", ncol=4,nrow = 1, align = "v",  common.legend = TRUE,legend = "bottom")
TSS.SR.p

TSS.SR.p2 <- ggplot(subset(TSS.df, dataset %in% "Test"),aes(predicted, residual, fill = ML)) + geom_point(size=2.5,stroke=0.5,shape = 21,alpha = 0.8) + theme_bw() + 
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]), 
                    guide = guide_legend(title = "ML Model", title.position = "top",override.aes = list(shape = 21), order = 1))+
  labs(x = "Predicted TSS (mg/L)", y = "Standardized Residual")+ ylim(-6,6)+
  theme(text = element_text(size=15),legend.position = "none", legend.direction='horizontal',
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        legend.text=element_text(size=8),legend.title=element_text(size=8))

Typlot <- ggplot(subset(TSS.df, dataset %in% "Test"), aes(x=residual, color = ML, fill = ML)) + geom_histogram(alpha = 0.5, position="dodge", binwidth = 0.5) + labs(y = "count") + xlim(-6,6)+
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]),guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  scale_color_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]),guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  coord_flip() + theme_bw() + theme(text = element_text(size=15),legend.position="none",axis.title.y=element_blank(),
                                    axis.text.y=element_blank(),axis.ticks.y=element_blank(), plot.margin = unit(c(0.38,1,0.38,-0.5), "lines"))
Tsr.p <- ggarrange(TSS.SR.p2, Typlot, ncol = 2, nrow = 1, widths = c(3.5, 1), common.legend = F)
Tsr.p


#### E.coli Predict ####
set.seed(1234)
tune.method <- trainControl(method = "repeatedcv", number=15, repeats=3, selectionFunction = "best")
e.mod.fun <- function(mod.,inputT,...){
  caret::train(E.coli ~ Group+Turb+NH4+NO3+Color+EC+pH
               , inputT, method = mod., trControl = tune.method, na.action = na.pass, preProc = c("center", "scale", "nzv"
               ),...)
}

set.seed(1234); WQ.E.lm <- #e.mod.fun("lm",WQ2.train)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_EST/E.coli_lm.rds")
set.seed(1234); WQ.E.pls <- #e.mod.fun("pls",WQ2.train,tuneLength = 5)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_EST/E.coli_pls.rds")
set.seed(1234); WQ.E.rf <- #e.mod.fun("rf",WQ2.train,ntree=25)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_EST/E.coli_rf.rds")
set.seed(1234); WQ.E.knn <- #e.mod.fun("knn",WQ2.train,tuneGrid = expand.grid (k = seq(from = 3, to = 10, by = 1)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_EST/E.coli_knn.rds")
set.seed(1234); WQ.E.svm <- #e.mod.fun("svmRadial",WQ2.train,tuneGrid=expand.grid(.sigma=seq(from=0.01,to=0.05,by=0.01),.C=seq(from=1,to=5,by=0.5)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_EST/E.coli_svm.rds")
set.seed(1234); WQ.E.cub <- #e.mod.fun("cubist",WQ2.train,tuneGrid=expand.grid(.committees=seq(from=3,to=10,by=1),.neighbors=seq(from=3,to=9,by=1)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_EST/E.coli_cub.rds")
set.seed(1234); WQ.E.qrnn <- #e.mod.fun("qrnn",WQ2.train,tuneGrid=expand.grid(.n.hidden=seq(from=2,to=5,by=1),.penalty = 10^(seq(from = -2, to = 2, by = 0.5)),.bag=FALSE))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_EST/E.coli_qrnn.rds")

#saveRDS(WQ.E.qrnn,"G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_EST/E.coli_qrnn.rds")


###  E.coli Train Variables Importance ### 
E.lm.Imp <- ggplot2::ggplot(varImp(WQ.E.lm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "E. coli-LM")
E.qrnn.Imp <- ggplot2::ggplot(varImp(WQ.E.qrnn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "E. coli-QRNN")
E.pls.Imp <- ggplot2::ggplot(varImp(WQ.E.pls, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "E. coli-PLS")
E.rf.Imp <- ggplot2::ggplot(varImp(WQ.E.rf, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "E. coli-RF")
E.knn.Imp <- ggplot2::ggplot(varImp(WQ.E.knn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "E. coli-KNN")
E.svm.Imp <- ggplot2::ggplot(varImp(WQ.E.svm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "E. coli-SVR")
E.cub.Imp <- ggplot2::ggplot(varImp(WQ.E.cub, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "E. coli-CUB")
E.varImp <- ggarrange(E.pls.Imp, E.knn.Imp, E.svm.Imp, E.rf.Imp, E.cub.Imp, E.qrnn.Imp,ncol=3, nrow = 2, align = "v", labels = "AUTO", common.legend = TRUE, legend="bottom")
E.varImp
#ggsave("E:/Dropbox/R project/Machine Learning/Figure/E_varImp.jpg",plot=E.varImp,width = 6, height = 4, dpi = 600)

### E.coli Train RMSE Boxplot ### 
E.results <- resamples(list("qrnn" = WQ.E.qrnn, "knn" = WQ.E.knn, "svr" = WQ.E.svm,"pls" = WQ.E.pls, "rf" = WQ.E.rf, "CUB" = WQ.E.cub))
E.RMSE <- data.frame(E.results$values$`qrnn~RMSE`,E.results$values$`knn~RMSE`,E.results$values$`svr~RMSE`,
                     E.results$values$`pls~RMSE`,E.results$values$`rf~RMSE`,E.results$values$`CUB~RMSE`)
colnames(E.RMSE) <- c("QRNN","KNN","SVR","PLS","RF","CUB")
E.RMSE <- E.RMSE %>% gather(key = "Model", value = "RMSE", QRNN, KNN, SVR, PLS, RF, CUB) %>% convert_as_factor(Model)
E.Ttest <- E.RMSE %>% pairwise_t_test(RMSE ~ Model, paired = TRUE, alt = c("two.sided"), conf.level = 0.99, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Model")
E.RMSE.p <- ggboxplot(E.RMSE, x = "Model", y = "RMSE",color = "black", legend = "none", order = c("PLS","KNN","SVR","RF","QRNN","CUB"),outlier.shape = NA) +
  geom_point(data = E.RMSE, aes(x = Model, y = RMSE, colour = Model, alpha = 0.9), position = position_jitter(width = 0.2))+
  geom_hline(yintercept = mean(E.RMSE$RMSE), linetype = 2) + 
  coord_flex_cart(bottom=brackets_horizontal(), left=capped_vertical('none'))+
  theme_bw()+labs(x = element_blank(), y = "RMSE (E.coli MPN/100ml)")+
  theme(text = element_text(size=15),legend.title=element_blank(),legend.position = "none",
        panel.border=element_blank(), axis.line = element_line(),legend.background = element_rect(colour='grey')) +
  scale_color_manual(values = c("PLS"=colx(6)[1],"KNN"=colx(6)[2],"SVR"=colx(6)[3],"RF"=colx(6)[4],"QRNN"=colx(6)[5],"CUB"=colx(6)[6]),guide=guide_legend(title = "Algorithms", order = 1))
E.RMSE.p

### TSS Test RMSE Boxplot ### 
e.coli.theme <- function(titl.,lab.data,...){
  list(geom_point(aes(fill=dataset), size=7.25, shape=21, color="black", stroke=0.5, alpha=0.75),
       theme_bw(),scale_fill_manual(values = c(colx2(2)), guide = guide_legend(title = NULL, title.position = "top",override.aes = list(shape=21,size=2.5), order = 1)),
       theme(text = element_text(size=8.5),axis.title.y = element_text(size=7.5),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = c(0.15, 0.85)),
       labs(x = "Predicted E.coli (MPN/100ml)", y = "Measured E.coli (MPN/100ml)",title = titl.),geom_abline(color = "black", linetype = 2, size = 1, alpha = 0.5),
       scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(10^-1,10^8),
                     labels = scales::trans_format("log10", scales::math_format (10^.x))),
       scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(10^-1,10^8),
                     labels = scales::trans_format("log10", scales::math_format (10^.x))),
       geom_richtext(data = lab.data,aes(10^-1, 10^6.5, label = label[1]),hjust = 0,size = 2.6,fill = "white", label.color = "black"),
       ...)
}

### LM ###
E.lm.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$E.coli, abs(predict.train(WQ.E.lm, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
E.lm.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.E.lm, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),E.lm.df)
E.lm.df <- E.lm.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
E.lm.sr <- lm(X2 ~ X3, data=E.lm.df)
E.lm.df <- cbind(E.lm.df, rstandard(E.lm.sr))
colnames(E.lm.df) <- c("Sample","actual","predicted","dataset","residual")
E.lm.df$dataset <-factor(E.lm.df$dataset  , levels=c("Train","Test"))
E.lm.df$ML <- "LM"
E.lm.L <- E.lm.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
E.lm.L$label[1] <- paste("Train RMSE =", round(E.lm.L$RMSE[2], 0), "<br> Test RMSE =", round(E.lm.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(E.lm.L$R2[1], 2), "<br> MAPE =", round(E.lm.L$MAPE[1], 1))
E.lm.p <- ggplot(E.lm.df, aes(predicted, actual))+e.coli.theme(NULL,E.lm.L)+theme(legend.position="none")
E.lm.p

### QRNN ###
E.qrnn.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.E.qrnn, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
E.qrnn.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.E.qrnn, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),E.qrnn.df)
E.qrnn.df <- E.qrnn.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
E.qrnn.sr <- lm(X2 ~ X3, data=E.qrnn.df)
E.qrnn.df <- cbind(E.qrnn.df, rstandard(E.qrnn.sr))
colnames(E.qrnn.df) <- c("Sample","actual","predicted","dataset","residual")
E.qrnn.df$dataset <-factor(E.qrnn.df$dataset  , levels=c("Train","Test"))
E.qrnn.df$ML <- "QRNN"
E.qrnn.L <- E.qrnn.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
E.qrnn.L$label[1] <- paste("Train RMSE =", round(E.qrnn.L$RMSE[2], 0), "<br> Test RMSE =", round(E.qrnn.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(E.qrnn.L$R2[1], 2), "<br> MAPE =", round(E.qrnn.L$MAPE[1], 1))
E.qrnn.p <- ggplot(E.qrnn.df, aes(predicted, actual))+e.coli.theme(NULL,E.qrnn.L)+theme(legend.position="none")
E.qrnn.p

### PLS ###
E.pls.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.E.pls, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
E.pls.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.E.pls, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),E.pls.df)
E.pls.df <- E.pls.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
E.pls.sr <- lm(X2 ~ X3, data=E.pls.df)
E.pls.df <- cbind(E.pls.df, rstandard(E.pls.sr))
colnames(E.pls.df) <- c("Sample","actual","predicted","dataset","residual")
E.pls.df$dataset <-factor(E.pls.df$dataset  , levels=c("Train","Test"))
E.pls.df$ML <- "PLS"
E.pls.L <- E.pls.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
E.pls.L$label[1] <- paste("Train RMSE =", round(E.pls.L$RMSE[2], 0), "<br> Test RMSE =", round(E.pls.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(E.pls.L$R2[1], 2), "<br> MAPE =", round(E.pls.L$MAPE[1], 1))
E.pls.p <- ggplot(E.pls.df, aes(predicted, actual))+e.coli.theme(NULL,E.pls.L)+theme(legend.position="none")
E.pls.p

### KNN ###
E.knn.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.E.knn, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
E.knn.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.E.knn, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),E.knn.df)
E.knn.df <- E.knn.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
E.knn.sr <- lm(X2 ~ X3, data=E.knn.df)
E.knn.df <- cbind(E.knn.df, rstandard(E.knn.sr))
colnames(E.knn.df) <- c("Sample","actual","predicted","dataset","residual")
E.knn.df$dataset <-factor(E.knn.df$dataset  , levels=c("Train","Test"))
E.knn.df$ML <- "KNN"
E.knn.L <- E.knn.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
E.knn.L$label[1] <- paste("Train RMSE =", round(E.knn.L$RMSE[2], 0), "<br> Test RMSE =", round(E.knn.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(E.knn.L$R2[1], 2), "<br> MAPE =", round(E.knn.L$MAPE[1], 1))
E.knn.p <- ggplot(E.knn.df, aes(predicted, actual))+e.coli.theme(NULL,E.knn.L)+theme(legend.position="none")
E.knn.p

### SVM ###
E.svm.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.E.svm, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
E.svm.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.E.svm, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),E.svm.df)
E.svm.df <- E.svm.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
E.svm.sr <- lm(X2 ~ X3, data=E.svm.df)
E.svm.df <- cbind(E.svm.df, rstandard(E.svm.sr))
colnames(E.svm.df) <- c("Sample","actual","predicted","dataset","residual")
E.svm.df$dataset <-factor(E.svm.df$dataset  , levels=c("Train","Test"))
E.svm.df$ML <- "SVR"
E.svm.L <- E.svm.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
E.svm.L$label[1] <- paste("Train RMSE =", round(E.svm.L$RMSE[2], 0), "<br> Test RMSE =", round(E.svm.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(E.svm.L$R2[1], 2), "<br> MAPE =", round(E.svm.L$MAPE[1], 1))
E.svm.p <- ggplot(E.svm.df, aes(predicted, actual))+e.coli.theme(NULL,E.svm.L)+theme(legend.position="none")
E.svm.p

### RF ###
E.rf.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.E.rf, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
E.rf.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.E.rf, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),E.rf.df)
E.rf.df <- E.rf.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
E.rf.sr <- lm(X2 ~ X3, data=E.rf.df)
E.rf.df <- cbind(E.rf.df, rstandard(E.rf.sr))
colnames(E.rf.df) <- c("Sample","actual","predicted","dataset","residual")
E.rf.df$dataset <-factor(E.rf.df$dataset  , levels=c("Train","Test"))
E.rf.df$ML <- "RF"
E.rf.L <- E.rf.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
E.rf.L$label[1] <- paste("Train RMSE =", round(E.rf.L$RMSE[2], 0), "<br> Test RMSE =", round(E.rf.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(E.rf.L$R2[1], 2), "<br> MAPE =", round(E.rf.L$MAPE[1], 1))
E.rf.p <- ggplot(E.rf.df, aes(predicted, actual))+e.coli.theme(NULL,E.rf.L)+theme(legend.position="none")
E.rf.p

### CUB ###
E.cub.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$TSS, abs(predict.train(WQ.E.cub, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
E.cub.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$TSS, abs(predict.train(WQ.E.cub, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),E.cub.df)
E.cub.df <- E.cub.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
E.cub.sr <- lm(X2 ~ X3, data=E.cub.df)
E.cub.df <- cbind(E.cub.df, rstandard(E.cub.sr))
colnames(E.cub.df) <- c("Sample","actual","predicted","dataset","residual")
E.cub.df$dataset <-factor(E.cub.df$dataset  , levels=c("Train","Test"))
E.cub.df$ML <- "CUB"
E.cub.L <- E.cub.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
E.cub.L$label[1] <- paste("Train RMSE =", round(E.cub.L$RMSE[2], 0), "<br> Test RMSE =", round(E.cub.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(E.cub.L$R2[1], 2), "<br> MAPE =", 24.8)
E.cub.p <- ggplot(E.cub.df, aes(predicted, actual))+e.coli.theme(NULL,E.cub.L)+theme(legend.position="none")
E.cub.p

### Save E.coli Result Scatter Plot ###
E.coli.R.p <- ggarrange(E.pls.p, E.svm.p, E.cub.p, E.qrnn.p, labels = "AUTO", ncol=4,nrow = 1, align = "v")
ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/Figure5-3.jpg",plot=E.coli.R.p,width = 12, height = 2, dpi = 600)

### E.coli Overall Scatter Plot ###
E.coli.df <- rbind(E.qrnn.df,E.pls.df,E.svm.df,E.cub.df)
E.coli.df$ML <-factor(E.coli.df$ML , levels=c("PLS","SVR","CUB","QRNN"))
E.coli.df.L <- rbind(E.pls.L,E.svm.L,E.cub.L,E.qrnn.L)
E.coli.df.L$ML <- c("PLS","PLS","SVR","SVR","CUB","CUB","QRNN","QRNN")
E.coli.df.L$label <- gsub("<br>", ";", E.coli.df.L$label)

E.coli.df.p <- ggplot(subset(E.coli.df, dataset %in% "Test"),aes(abs(predicted), abs(actual), fill = ML)) + geom_point(size=2.5,stroke=0.5,shape = 21,alpha = 0.8) + theme_bw() + 
  geom_abline(color = "black", linetype = 2, size = 1, alpha = 0.4)+
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]), 
                    guide = guide_legend(title = "ML Model", title.position = "top",override.aes = list(shape = 21), order = 1))+
  scale_shape_manual(values = c(16,16,16,16,16), guide = guide_legend(title = "Sampling Location", title.position = "top", order = 2))+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(1,10^8.5),
                labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(1,10^8.5),
                labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  labs(x = "Predicted E.coli (MPN/100ml)", y = "Measured E.coli (MPN/100ml)")+
  theme(text = element_text(size=15),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = "none", legend.box = "horizontal")
E.coli.df.p

### E.coli Result Standardized Residual Plot ###
e.sr.theme <- function(titl.,...){
  list(geom_point(aes(fill = Sample),shape=21,size=2,color="black", stroke = 0.5),theme_bw(),ylim(-3,3),
       scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(10^0,10^8),
                     labels = scales::trans_format("log10", scales::math_format (10^.x))),
       scale_alpha_continuous(limits = c(0,3), breaks = c(0,1,2,3), guide = guide_legend(title = "Standardized residual", title.position = "top", order = 2)),
       scale_fill_manual(values = c(colx2(5)), 
                         guide = guide_legend(title = "Sampling location", title.position = "left",override.aes = list(shape = 21), order = 1)),
       theme(text = element_text(size=7),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = c(.9, 0.3)),
       labs(title = titl. ,x = "Predicted E.coli (MPN/100ml)", y = "Standardized Residual"))
}

E.qrnn.p2 <- ggplot(data = E.qrnn.df, aes(predicted, residual)) + e.sr.theme("E.coli: QRNN") + guides(alpha=FALSE)
E.pls.p2 <- ggplot(data = E.pls.df, aes(predicted, residual)) + e.sr.theme("E.coli: PLS") + guides(alpha=FALSE)
E.knn.p2 <- ggplot(data = E.knn.df, aes(predicted, residual)) + e.sr.theme("E.coli: KNN") + guides(alpha=FALSE)
E.rf.p2 <- ggplot(data = E.rf.df, aes(predicted, residual)) + e.sr.theme("E.coli: RF") + guides(alpha=FALSE)
E.svm.p2 <- ggplot(data = E.svm.df, aes(predicted, residual)) + e.sr.theme("E.coli: SVR") + guides(alpha=FALSE)
E.cub.p2 <- ggplot(data = E.cub.df, aes(predicted, residual)) + e.sr.theme("E.coli: CUB") + guides(alpha=FALSE)
E.coli.SR.p <- ggarrange(E.pls.p2, E.svm.p2, E.cub.p2, E.qrnn.p2, labels = "AUTO", ncol=4,nrow = 1, align = "v",  common.legend = TRUE,legend = "bottom")
E.coli.SR.p

E.coli.SR.p2 <- ggplot(subset(E.coli.df, dataset %in% "Test"),aes(predicted, residual, fill = ML)) + geom_point(size=2.5,stroke=0.5,shape = 21,alpha = 0.8) + theme_bw() + 
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]), 
                    guide = guide_legend(title = "ML Model", title.position = "top",override.aes = list(shape = 21), order = 1))+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(1,10^8.5),
                labels = scales::trans_format("log10", scales::math_format (10^.x))) + ylim(-6,6)+
  labs(x = "Predicted E.coli (MPN/100ml)", y = "Standardized Residual")+
  theme(text = element_text(size=15),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = "none", legend.box = "horizontal")

Eyplot <- ggplot(subset(E.coli.df, dataset %in% "Test"), aes(x=residual, color = ML, fill = ML)) + geom_histogram(alpha = 0.5, position="dodge", binwidth = 0.5) + labs(y = "count") + xlim(-6,6)+
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]),guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  scale_color_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]),guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  coord_flip() + theme_bw() + theme(text = element_text(size=15),legend.position="none",axis.title.y=element_blank(),
                                    axis.text.y=element_blank(),axis.ticks.y=element_blank(), plot.margin = unit(c(0.38,1,0.62,-0.5), "lines"))
Esr.p <- ggarrange(E.coli.SR.p2, Eyplot, ncol = 2, nrow = 1, widths = c(3.5, 1), common.legend = F)
Esr.p


#### Predict Result Plot ####
Result.p <- ggarrange(C.pls.p, C.rf.p, C.cub.p, C.qrnn.p,
                      T.pls.p, T.rf.p, T.cub.p, T.qrnn.p,
                      E.pls.p, E.rf.p, E.cub.p, E.qrnn.p,
                      labels = "AUTO", ncol=4,nrow = 3, align = "v")
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/Figure5.tiff",plot=Result.p,width = 12, height = 6, dpi = 600)

p2_legend <- get_legend(C.pls.p+theme(legend.text = element_text(size=14),legend.title = element_text(size=14), legend.key.height= unit(1, 'cm'),legend.key.width= unit(1, 'cm'),legend.key=element_blank()))
p2_legend <- ggarrange(NULL, p2_legend, NULL, ncol = 3, widths = c(2.5,2,1))
Result.p <- ggarrange(C.pls.p+theme(legend.position="none")+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()),
                      T.pls.p+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()),
                      E.pls.p+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()), 
                      C.rf.p+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()),
                      T.rf.p+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()),
                      E.rf.p+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()), 
                      C.cub.p+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()), 
                      T.cub.p+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()), 
                      E.cub.p+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()), 
                      C.qrnn.p, T.qrnn.p, E.qrnn.p,
                      labels = c("A", "E", "I", "B", "F", "J", "C", "G", "K", "D", "H", "L"),
                      ncol=3,nrow = 4, align = "v",  common.legend = F, heights = c(1,1,1,1.2))

Result.p <- ggarrange(Result.p, p2_legend, nrow=2, heights=c(20, 1))
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/Figure5.tiff",plot=Result.p,width = 18/2, height = 17/2, dpi = 600)

Result.I <- ggarrange(C.pls.Imp, C.svm.Imp, C.cub.Imp, C.qrnn.Imp,
                      T.pls.Imp, T.svm.Imp, T.cub.Imp, T.qrnn.Imp,
                      E.pls.Imp, E.svm.Imp, E.cub.Imp, E.qrnn.Imp,
                      labels = "AUTO", ncol=4,nrow = 3, align = "v")
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/Importance_Result.jpg",plot=Result.I,width = 12, height = 6, dpi = 600)


Result.p <- ggarrange(COD.df.p,TSS.df.p, E.coli.df.p,
                      labels = "AUTO", ncol=3, nrow = 1, align = "hv",  common.legend = F)
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/newFig5.jpg",plot=Result.p,width = 12*1.5, height = 2.7*1.5, dpi = 600)

Result2.p <- ggarrange(Csr.p,Tsr.p, Esr.p,
                       labels = "AUTO", ncol=3, nrow = 1, align = "hv",  common.legend = F)
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/newFig5sd.tiff",plot=Result2.p,width = 12*1.5, height = 2.7*1.5, dpi = 600)


#### E.coli Classification ####
set.seed(1234)
tune.method <- trainControl(method = "repeatedcv", number=15, repeats=3, selectionFunction = "best")
e2.mod.fun <- function(mod.,inputT,...){
  caret::train(E.coli2 ~ Group+Turb+NH4+NO3+Color+EC+pH
               , inputT, method = mod., trControl = tune.method, na.action = na.pass, preProc = c("center", "scale", "nzv"
               ),...)
}

set.seed(1234); WQ.E2.rf <- #e2.mod.fun("rf",WQ2.train,ntree=25)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_CLS/E.coli2_rf.rds")
set.seed(1234); WQ.E2.knn <- #e2.mod.fun("knn",WQ2.train,tuneGrid = expand.grid (k = seq(from = 3, to = 9, by = 1)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_CLS/E.coli2_knn.rds")
set.seed(321); WQ.E2.svm <- #e2.mod.fun("svmRadial",WQ2.train,tuneGrid=expand.grid(.sigma=seq(from=0.09,to=0.15,by=0.01),.C=seq(from=1,to=5,by=0.5)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_CLS/E.coli2_svm.rds")
#set.seed(1234); WQ.E2.xgb <- e2.mod.fun("xgbLinear",WQ2.train)
#saveRDS(WQ.E2.svm,"G:/My Drive/R project/GitHub/MLsensor/Save_model/Ecoli_CLS/E.coli2_svm.rds")

Rf_pred <- predict(WQ.E2.rf, newdata = WQ2.test)
Rf_e2_result <- confusionMatrix(Rf_pred, WQ2.test$E.coli2)
E2_R_df <- data.frame(Rf_e2_result$byClass, Class = c("H","M","L","LDL"), Model = "RF") 


knn_pred <- predict(WQ.E2.knn, newdata = WQ2.test)
knn_e2_result <- confusionMatrix(knn_pred, WQ2.test$E.coli2)
E2_R_df <- rbind(E2_R_df, data.frame(knn_e2_result$byClass, Class = c("H","M","L","LDL"), Model = "KNN") )


svm_pred <- predict(WQ.E2.knn, newdata = WQ2.test)
svm_e2_result <- confusionMatrix(svm_pred, WQ2.test$E.coli2)
E2_R_df <- rbind(E2_R_df, data.frame(svm_e2_result$byClass, Class = c("H","M","L","LDL"), Model = "SVM") )

E2_R_df <- E2_R_df %>% 
  reshape2::melt(id.vars=c("Model", "Class"))

E2_R_df$Class <- factor(E2_R_df$Class , levels=c("H","M","L","LDL"))
E2_R_df$Model <- factor(E2_R_df$Model , levels=c("RF","KNN","SVM"))
E2_R_df$variable <- ifelse(E2_R_df$variable == "Balanced.Accuracy","Balanced Accuracy", as.character(E2_R_df$variable))
E2plot <- ggplot(subset(E2_R_df, variable %in% c("Balanced Accuracy", "Precision", "Recall", "F1")) , aes(x = Class, y = value, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ variable, scales = "free_x") +
  labs(x = element_blank(), y =element_blank(), fill = element_blank()) +
  scale_fill_manual(values = c("RF" = colx(3)[1], "KNN" = colx(3)[2], "SVM" = colx(3)[3])) +
  theme_bw()
E2plot