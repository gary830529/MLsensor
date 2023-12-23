pacman::p_load(plyr,dplyr,reshape2,stringr,lubridate,readxl,openxlsx,janitor,ggplot2,ggpmisc,ggpubr,ggsignif,gplots,data.table,magrittr,grid,
               tidyverse,rstatix,PtProcess,gtable,gridExtra,cowplot,jpg,tiff,pwr,MOTE,kableExtra,tinytex,knitr,quantmod,class,caret,gmodels,
               officedown,officer,flextable,car,C50,GGally,ggResidpanel,ggfortify,caretEnsemble,randomForest,corrplot,neuralnet,lemon,Hmisc,scales,
               glue,ggtext,png,gtools,ggrepel,rvg,gdata,scales,nnet,xgboost,correlation,psych,corrgram,jpeg)



#### Clean Data ####
Raw <- read.csv(file = 'G:/My Drive/R project/GitHub/MLsensor/Raw_data/Raw_data.csv',header = T, sep = ",")
Raw$Date <- ymd(Raw$Date)
Raw[Raw < 0] <- NaN
Raw <- Raw[Raw$Date <= ymd("2020-03-10"), ]
Raw <-  Raw %>% filter(Group %in% c("S1","S2","S3","S4","S5"))
colnames(Raw) <- c("X","Group","Date","COD","sCOD","TSS","E.coli","Color","Turb","EC","pH","NH4","NO3","Temp","VSS")

#### Start ####
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

#### Restart remove data ####
Dates_after_events2 <- as.Date(c("2019-01-15", "2019-01-22"))
WQx <- WQ %>% filter(Date %nin% Dates_after_events2)

WQ.Cp <- WQx[complete.cases(WQx), ]
WQ.Cp <- WQ.Cp %>% subset(select=-c(Date,Event,BP,VSS))
WQ.Cp$Group <- as.numeric(as.factor(WQ.Cp$Group))
WQ.Cp$E.coli2 <- ifelse(WQ.Cp$E.coli >= 10^6, "H", ifelse(WQ.Cp$E.coli < 10^6 & WQ.Cp$E.coli >= 10^4, "M", WQ.Cp$E.coli)) 
WQ.Cp$E.coli2 <- ifelse(WQ.Cp$E.coli < 10^4 & WQ.Cp$E.coli >= 2, "L", WQ.Cp$E.coli2)
WQ.Cp$E.coli2 <- ifelse(WQ.Cp$E.coli < 2, "LDL", WQ.Cp$E.coli2)
WQ.Cp$E.coli2 <- factor(WQ.Cp$E.coli2 , levels=c("H","M","L","LDL"))

#### Mean Absolute Percentage Error #### 
MAPE <- function(pred.,actu.,...){
  mean(abs((ifelse(actu. == 0, (actu.+pred.)/2, actu.)-pred.)/ifelse(actu. > pred., ifelse(actu. == 0, (actu.+pred.)/2, actu.),pred.) ))*100
  #    ifelse(abs((actu.-pred.)/actu.)>1,NA,abs((actu.-pred.)/actu.)), na.rm=TRUE) * 100
}


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

#### COD Boxplot ####
C.Bplot <- 
  #ggboxplot(WQ, x = "Group", y = "COD", add = "jitter", legend = "none", add.params = list(size = 0.75), order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  #Bplot("COD (mg/L)",ggforce::facet_row(vars(BP), scales = 'free', space = 'free')) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff"))
  
  ggboxplot(WQ, x = "Group", y = "COD", add = "jitter", legend = "none", color = "Event", shape = "Event", add.params = list(size = 0.75),
            order = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")) + 
  Bplot("COD (mg/L)"
        #,ggforce::facet_row(vars(BP), scales = 'free', space = 'free')
  ) + scale_x_discrete(labels=c("Influent"="Inf","AnMBR"="AnMBR","Permeate"="Perm","Post-NCS"="P-NCS","Effluent"="Eff")) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format (10^.x)))











