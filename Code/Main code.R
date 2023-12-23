pacman::p_load(plyr,dplyr,reshape2,stringr,lubridate,readxl,openxlsx,janitor,ggplot2,ggpmisc,ggpubr,ggsignif,gplots,data.table,magrittr,grid,
               tidyverse,rstatix,PtProcess,gtable,gridExtra,cowplot,jpg,tiff,pwr,MOTE,kableExtra,tinytex,knitr,quantmod,class,caret,gmodels,
               officedown,officer,flextable,car,C50,GGally,ggResidpanel,ggfortify,caretEnsemble,randomForest,corrplot,neuralnet,lemon,Hmisc,scales,
               glue,ggtext,png,gtools,ggrepel,rvg,gdata,scales,nnet,xgboost,correlation,psych,corrgram,jpeg)
#### Start ####
Raw <- read.csv(file = 'G:/My Drive/R project/GitHub/MLsensor/Raw_data/Raw_data.csv',header = T, sep = ",")
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

# COD Boxplot #
C.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "COD", add = "jitter", legend = "none", add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("COD (mg/L)",ggforce::facet_row(vars(BP), scales = 'free', space = 'free')) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))
  ggboxplot(WQ, x = "Group", y = "COD", add = "jitter", legend = "none", color = "Event", shape = "Event", add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("COD (mg/L)"#,ggforce::facet_row(vars(BP), scales = 'free', space = 'free')
  ) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff")) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format (10^.x)))

# TSS Boxplot #
T.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "TSS", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("TSS (mg/L)",ggforce::facet_row(vars(BP), scales = 'free', space = 'free')) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))
  ggboxplot(WQ, x = "Group", y = "TSS", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("TSS (mg/L)"#,ggforce::facet_row(vars(BP), scales = 'free', space = 'free')
  ) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff")) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format (10^.x)))

# E.coli Boxplot #
E.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "E.coli", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("E.coli (MPN/100ml)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff")) +
  #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format (10^.x)))
  ggboxplot(WQ, x = "Group", y = "E.coli", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("E.coli (MPN/100ml)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff")) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format (10^.x)))

# Color Boxplot #
co.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "Color", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("Color (Pt/Co)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))+ylim(NA,7500)
  ggboxplot(WQ, x = "Group", y = "Color", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("Color (Pt/Co)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))+ylim(NA,7500)

# Turbidity Boxplot #
tu.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "Turb", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("Turbidity (NTU)",ggforce::facet_row(vars(BP), scales = 'free', space = 'free')) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))
  ggboxplot(WQ, x = "Group", y = "Turb", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("Turbidity (NTU)",ggforce::facet_row(vars(BP), scales = 'free', space = 'free')) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))

# EC Boxplot #
ec.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "EC", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("EC (µs/cm)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))
  ggboxplot(WQ, x = "Group", y = "EC", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("EC (µs/cm)") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))

# pH Boxplot #
ph.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "pH", add = "jitter", legend = "none",add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("pH") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))
  ggboxplot(WQ, x = "Group", y = "pH", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("pH") + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))

# Temp Boxplot #
tbp <- WQ %>% filter(Group == c("Influent","AnMBR"))
tbp$Group <- ifelse(tbp$Group == "Influent", "Ambient air", "AnMBR")
temp.Bplot <- 
  #ggboxplot(tbp, x = "Group", y = "Temp", add = "jitter", legend = "none",add.params = list(size = .75), order = c("Ambient air", "AnMBR")) + Bplot("°C") + theme(legend.position = "right") + 
  #guides(color = guide_legend(label.position = "bottom"))
  ggboxplot(tbp, x = "Group", y = "Temp", add = "jitter", legend = "none", color = "Event", shape = "Event",add.params = list(size = .75),
            order = c("Ambient air", "AnMBR")) + Bplot("°C") + theme(legend.position = "right") + 
  guides(color = guide_legend(label.position = "bottom"))

# Boxplot #
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

# Include SCOD & PCOD #
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

#### COD REF Models ####
rfe.C.lm <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/C_rfe_lm.rds")
#c.rfe.fun("lm",WQ.Cp)
rfe.C.pls <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/C_rfe_pls.rds")
#c.rfe.fun("pls",WQ.Cp,tuneLength = 5)
rfe.C.rf <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/C_rfe_rf.rds")
#c.rfe.fun("rf",WQ.Cp)
rfe.C.knn <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/C_rfe_knn.rds")
#c.rfe.fun("knn",WQ.Cp)
rfe.C.svm <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/C_rfe_svm.rds")
#c.rfe.fun("svmRadial",WQ.Cp)
rfe.C.cub <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/C_rfe_cub.rds")
#c.rfe.fun("cubist",WQ.Cp)
rfe.C.net <- readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/C_rfe_rnn.rds")
#c.rfe.fun("qrnn",WQ.Cp)

#saveRDS(rfe.C.net,"G:/My Drive/R project/GitHub/MLsensor/Save_model/C_rfe_rnn.rds")

# REF Plot #
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




