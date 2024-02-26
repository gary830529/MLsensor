pacman::p_load(plyr,dplyr,reshape2,stringr,lubridate,readxl,openxlsx,janitor,ggplot2,ggpmisc,ggpubr,ggsignif,gplots,data.table,magrittr,grid,
               tidyverse,rstatix,PtProcess,gtable,gridExtra,cowplot,tiff,pwr,MOTE,kableExtra,tinytex,knitr,quantmod,class,caret,gmodels,
               officedown,officer,flextable,car,C50,GGally,ggResidpanel,ggfortify,caretEnsemble,randomForest,corrplot,neuralnet,lemon,Hmisc,scales,
               glue,ggtext,png,gtools,ggrepel,rvg,gdata,scales,nnet,xgboost,correlation,psych,corrgram,jpeg,qrnn,pls,kernlab)
colx <- colorRampPalette(c("#16068a","#9e189d","#fdb32e"))
colx2 <- colorRampPalette(c("#264653","#2a9d8e","#e9c46b","#f3a261","#e66f51"))
setwd("G:/My Drive/R project")

#### Start ####
Raw <- read.csv(file = 'G:/My Drive/R project/GitHub/MLsensor/Raw_data/Data_ML.csv',header = T, sep = ",")
#pCOD was calculated by COD - sCOD
Raw$Date <- ymd(Raw$Date)
Raw[Raw < 0] <- NaN
Raw <- Raw[Raw$Date <= ymd("2020-03-10"), ]
Raw <-  Raw %>% filter(Group %in% c("S1","S2","S3","S4","S5"))
colnames(Raw) <- c("X","Group","Date","COD","sCOD","TSS","E.coli","Color","Turb","EC","pH","NH4","NO3","Temp","VSS")

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


Dates_after_events2 <- as.Date(c("2019-01-15", "2019-01-22"))
WQ2 <- WQ %>% filter(Date %nin% Dates_after_events2)
WQ2 <- WQ2[complete.cases(WQ2), ]
WQ2 <- WQ2 %>% subset(select=-c(Date,Event,BP,VSS))
WQ2$Group <- as.numeric(as.factor(WQ2$Group))


set.seed(1234)
WQ2.index <- sample(1:nrow(WQ2),size=nrow(WQ2)*0.7,replace = FALSE)
WQ2.train <- WQ2[WQ2.index,]
WQ2.test <- WQ2[-WQ2.index,]


#### Correlation for COD pCOD sCOD ####
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

WQ2.cor <- WQ2 %>% within(Group <- factor(Group, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
colnames(WQ2.cor) <- c("Group", "COD", "sCOD","pCOD","TSS","E.coli","Color","Trubidity","EC","pH","NH4","NO3","Temperature") 


WQ2.Cor <- ggpairs(subset(WQ2.cor, select = -c(Group,E.coli)), 
                   upper = list(continuous = upper_fn), lower = list(continuous = lower_fn), diag = list(continuous = diag_fn))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.text.y = element_text(size=5),text = element_text(size=6))
WQ2.Cor
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/CorrelationP.jpg",plot=WQ2.Cor,width = 5, height = 5, dpi = 300)

#### Correlation between COD pCOD sCOD ####
upper_fn2 <- function(data, mapping, method="pearson", use="pairwise", ...){
  ggally_cor(data = data, mapping = mapping, size = 2.3, fontface = "bold", ...) + 
    theme_void()
}

diag_fn <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(alpha=.2,...)+theme_bw()
}

lower_fn2 <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha=0.3, size=0.9) +
    geom_smooth(method='lm', size=0.5, se = FALSE,...)+theme_bw()
}

WQ2.COD.cor <- ggpairs(WQ2.cor[1:5],aes(colour = Group, alpha = 0.4), title="COD Correlations",
                      upper = list(continuous = upper_fn2), lower = list(continuous = lower_fn2), diag = list(continuous = diag_fn))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.text.y = element_text(size=5),panel.grid.major = element_blank())+
  scale_fill_manual(values = c(colx2(5)))+
  scale_color_manual(values = c(colx2(5)))
WQ2.COD.cor
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/Correlation2P.jpg",plot=WQ2.COD.cor,width = 5, height = 5, dpi = 300)

#### sCOD Predict ####
set.seed(1234)
tune.method <- trainControl(method = "repeatedcv", number=15, repeats=3, selectionFunction = "best")
sc.mod.fun <- function(mod.,inputT,...){
  caret::train(sCOD ~ Group+Turb+NH4+NO3+Color+EC+pH
               , inputT, method = mod., trControl = tune.method, na.action = na.pass, preProc = c("center", "scale", "nzv"
               ),...)
}

set.seed(1234); WQ.sC.lm <- #sc.mod.fun("lm",WQ2.train)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_lm.rds")
set.seed(1234); WQ.sC.pls <- #sc.mod.fun("pls",WQ2.train,tuneLength = 10)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_pls.rds")
set.seed(1234); WQ.sC.rf <- #sc.mod.fun("rf",WQ2.train,ntree=25)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_rf.rds")
set.seed(1234); WQ.sC.knn <- #sc.mod.fun("knn",WQ2.train,tuneGrid = expand.grid (k = seq(from = 3, to = 10, by = 1)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_knn.rds")
set.seed(1234); WQ.sC.svm <- #sc.mod.fun("svmRadial",WQ2.train,tuneGrid=expand.grid(.sigma=seq(from=0.01,to=0.1,by=0.01),.C=seq(from = 25,to = 35, by = 0.5)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_svm.rds")
set.seed(1234); WQ.sC.cub <- #sc.mod.fun("cubist",WQ2.train,tuneGrid=expand.grid(.committees=seq(from=3,to=10,by=1),.neighbors = seq(from = 3, to = 9, by = 1)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_cub.rds")
set.seed(1234); WQ.sC.qrnn <- #sc.mod.fun("qrnn",WQ2.train,tuneGrid=expand.grid(.n.hidden=seq(from=2,to=5,by=1),.penalty = 10^(seq(from = -2, to = 2, by = 0.25)),.bag=FALSE))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_qrnn.rds")
#set.seed(1234); WQ.C.xg <- sc.mod.fun("xgbTree",WQ2.train)

#saveRDS(WQ.sC.svm,"G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_svm.rds")

###  sCOD Train Variables Importance ### 
sC.lm.Imp <- ggplot(varImp(WQ.sC.lm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "sCOD-LM")
sC.pls.Imp <-ggplot(varImp(WQ.sC.pls, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "sCOD-PLS")
sC.rf.Imp <- ggplot(varImp(WQ.sC.rf, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "sCOD-RF")
sC.knn.Imp <- ggplot(varImp(WQ.sC.knn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "sCOD-KNN")
sC.svm.Imp <- ggplot(varImp(WQ.sC.svm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "sCOD-SVR")
sC.cub.Imp <- ggplot(varImp(WQ.sC.cub, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "sCOD-CUB")
sC.qrnn.Imp <- ggplot(varImp(WQ.sC.qrnn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "sCOD-QRNN")
sC.varImp <- ggarrange(sC.lm.Imp,sC.pls.Imp, sC.svm.Imp,sC.rf.Imp, sC.cub.Imp,sC.qrnn.Imp, ncol=3, nrow = 2, align = "v", labels = "AUTO", common.legend = TRUE, legend="bottom")
sC.varImp
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/SC_varImp.jpg",plot=sC.varImp,width = 6, height = 4, dpi = 300)

### sCOD Train RMSE Boxplot ### 
sC.results <- resamples(list("LM" = WQ.sC.lm, "qrnn" = WQ.sC.qrnn, "knn" = WQ.sC.knn, "svr" = WQ.sC.svm,"pls" = WQ.sC.pls, "rf" = WQ.sC.rf, "CUB" = WQ.sC.cub))
sC.RMSE <- data.frame(sC.results$values$`LM~RMSE`,sC.results$values$`qrnn~RMSE`,sC.results$values$`knn~RMSE`,sC.results$values$`svr~RMSE`,
                     sC.results$values$`pls~RMSE`,sC.results$values$`rf~RMSE`,sC.results$values$`CUB~RMSE`)
colnames(sC.RMSE) <- c("LM","QRNN","KNN","SVR","PLS","RF","CUB")
sC.RMSE <- sC.RMSE %>% gather(key = "Model", value = "RMSE", LM, QRNN, KNN, SVR, PLS, RF, CUB) %>% convert_as_factor(Model)
sC.Ttest <- sC.RMSE %>% pairwise_t_test(RMSE ~ Model, paired = TRUE, alt = c("two.sided"), conf.level = 0.99, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Model")
sC.RMSE.p <- ggboxplot(sC.RMSE, x = "Model", y = "RMSE",color = "black", legend = "none", order = c("LM","PLS","KNN","SVR","RF","QRNN","CUB"),outlier.shape = NA) +
  geom_point(data = sC.RMSE, aes(x = Model, y = RMSE, colour = Model, alpha = 0.9), position = position_jitter(width = 0.2))+
  geom_hline(yintercept = mean(sC.RMSE$RMSE), linetype = 2) + 
  coord_flex_cart(bottom=brackets_horizontal(), left=capped_vertical('none'))+
  theme_bw()+labs(x = element_blank(), y = "RMSE (sCOD mg/L)")+
  theme(text = element_text(size=15),legend.title=element_blank(),legend.position = "none",
        panel.border=element_blank(), axis.line = element_line(),legend.background = element_rect(colour='grey')) +
  scale_color_manual(values = c("LM"=colx(7)[1],"PLS"=colx(7)[2],"KNN"=colx(7)[3],"SVR"=colx(7)[4],"RF"=colx(7)[5],"QRNN"=colx(7)[6],"CUB"=colx(7)[7]),guide=guide_legend(title = "Algorithms", order = 1))
sC.RMSE.p
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/SC_RMSE.jpg",plot=sC.RMSE.p,width = 6, height = 4, dpi = 300)

### sCOD Test RMSE Boxplot ### 
sCOD.theme <- function(titl.,lab.data,...){
  list(geom_point(aes(fill=dataset), size=3.25, shape=21, color="black", stroke=0.5, alpha=0.75),
       theme_bw(),ylim(0,6000),xlim(0,6000),
       scale_fill_manual(values = c(colx2(2)), guide = guide_legend(title = "Dataset", title.position = "top",override.aes = list(shape=21,size=2.5), order = 1, nrow = 2)),
       theme(text = element_text(size=8.5),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = c(0.15, 0.85)),
       labs(x = "Estimated sCOD (mg/L)", y = "Measured sCOD (mg/L)",title = titl.),geom_abline(color = "black", linetype = 2, size = 1, alpha = 0.5),
       geom_richtext(data = lab.data,aes(3400, 1000, label = label[1]),hjust = 0,size = 2.75,fill = "white", label.color = "black"),
       ...)
}

### LM ###
sC.lm.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$sCOD, abs(predict.train(WQ.sC.lm, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
sC.lm.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$sCOD, abs(predict.train(WQ.sC.lm, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),sC.lm.df)
sC.lm.df <- sC.lm.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
sC.lm.sr <- lm(X2 ~ X3, data=sC.lm.df)
sC.lm.df <- cbind(sC.lm.df, rstandard(sC.lm.sr))
colnames(sC.lm.df) <- c("Sample","actual","predicted","dataset","residual")
sC.lm.df$dataset <-factor(sC.lm.df$dataset  , levels=c("Train","Test"))
sC.lm.df$ML <- "LM"
sC.lm.L <- sC.lm.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
sC.lm.L$label[1] <- paste("Train RMSE =", round(sC.lm.L$RMSE[2], 0), "<br> Test RMSE =", round(sC.lm.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(sC.lm.L$R2[1], 2), "<br> MAPE =", round(sC.lm.L$MAPE[1], 1))
sC.lm.p <- ggplot(sC.lm.df, aes(predicted, actual))+sCOD.theme("sCOD-LM",sC.lm.L)+theme(legend.position="none")
sC.lm.p

### QRNN ###
sC.qrnn.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$sCOD, abs(predict.train(WQ.sC.qrnn, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
sC.qrnn.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$sCOD, abs(predict.train(WQ.sC.qrnn, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),sC.qrnn.df)
sC.qrnn.df <- sC.qrnn.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
sC.qrnn.sr <- lm(X2 ~ X3, data=sC.qrnn.df)
sC.qrnn.df <- cbind(sC.qrnn.df, rstandard(sC.qrnn.sr))
colnames(sC.qrnn.df) <- c("Sample","actual","predicted","dataset","residual")
sC.qrnn.df$dataset <-factor(sC.qrnn.df$dataset  , levels=c("Train","Test"))
sC.qrnn.df$ML <- "QRNN"
sC.qrnn.L <- sC.qrnn.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
sC.qrnn.L$label[1] <- paste("Train RMSE =", round(sC.qrnn.L$RMSE[2], 0), "<br> Test RMSE =", round(sC.qrnn.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(sC.qrnn.L$R2[1], 2), "<br> MAPE =", round(sC.qrnn.L$MAPE[1], 1))
sC.qrnn.p <- ggplot(sC.qrnn.df, aes(predicted, actual))+sCOD.theme("sCOD-QRNN",sC.qrnn.L)+theme(legend.position="none")
sC.qrnn.p

### PLS ###
sC.pls.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$sCOD, abs(predict.train(WQ.sC.pls, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
sC.pls.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$sCOD, abs(predict.train(WQ.sC.pls, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),sC.pls.df)
sC.pls.df <- sC.pls.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
sC.pls.sr <- lm(X2 ~ X3, data=sC.pls.df)
sC.pls.df <- cbind(sC.pls.df, rstandard(sC.pls.sr))
colnames(sC.pls.df) <- c("Sample","actual","predicted","dataset","residual")
sC.pls.df$dataset <-factor(sC.pls.df$dataset  , levels=c("Train","Test"))
sC.pls.df$ML <- "PLS"
sC.pls.L <- sC.pls.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
sC.pls.L$label[1] <- paste("Train RMSE =", round(sC.pls.L$RMSE[2], 0), "<br> Test RMSE =", round(sC.pls.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(sC.pls.L$R2[1], 2), "<br> MAPE =", round(sC.pls.L$MAPE[1], 1))
sC.pls.p <- ggplot(sC.pls.df, aes(predicted, actual))+sCOD.theme("sCOD-PLS",sC.pls.L)
sC.pls.p

### KNN ###
sC.knn.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$sCOD, abs(predict.train(WQ.sC.knn, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
sC.knn.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$sCOD, abs(predict.train(WQ.sC.knn, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),sC.knn.df)
sC.knn.df <- sC.knn.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
sC.knn.sr <- lm(X2 ~ X3, data=sC.knn.df)
sC.knn.df <- cbind(sC.knn.df, rstandard(sC.knn.sr))
colnames(sC.knn.df) <- c("Sample","actual","predicted","dataset","residual")
sC.knn.df$dataset <-factor(sC.knn.df$dataset  , levels=c("Train","Test"))
sC.knn.df$ML <- "KNN"
sC.knn.L <- sC.knn.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
sC.knn.L$label[1] <- paste("Train RMSE =", round(sC.knn.L$RMSE[2], 0), "<br> Test RMSE =", round(sC.knn.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(sC.knn.L$R2[1], 2), "<br> MAPE =", round(sC.knn.L$MAPE[1], 1))
sC.knn.p <- ggplot(sC.knn.df, aes(predicted, actual))+sCOD.theme("sCOD-KNN",sC.knn.L)+theme(legend.position="none")
sC.knn.p

### SVM ###
sC.svm.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$sCOD, abs(predict.train(WQ.sC.svm, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
sC.svm.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$sCOD, abs(predict.train(WQ.sC.svm, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),sC.svm.df)
sC.svm.df <- sC.svm.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
sC.svm.sr <- lm(X2 ~ X3, data=sC.svm.df)
sC.svm.df <- cbind(sC.svm.df, rstandard(sC.svm.sr))
colnames(sC.svm.df) <- c("Sample","actual","predicted","dataset","residual")
sC.svm.df$dataset <-factor(sC.svm.df$dataset  , levels=c("Train","Test"))
sC.svm.df$ML <- "SVR"
sC.svm.L <- sC.svm.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
sC.svm.L$label[1] <- paste("Train RMSE =", round(sC.svm.L$RMSE[2], 0), "<br> Test RMSE =", round(sC.svm.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(sC.svm.L$R2[1], 2), "<br> MAPE =", round(sC.svm.L$MAPE[1], 1))
sC.svm.p <- ggplot(sC.svm.df, aes(predicted, actual))+sCOD.theme("sCOD-SVR",sC.svm.L)+theme(legend.position="none")
sC.svm.p

### RF ###
sC.rf.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$sCOD, abs(predict.train(WQ.sC.rf, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
sC.rf.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$sCOD, abs(predict.train(WQ.sC.rf, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),sC.rf.df)
sC.rf.df <- sC.rf.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
sC.rf.sr <- lm(X2 ~ X3, data=sC.rf.df)
sC.rf.df <- cbind(sC.rf.df, rstandard(sC.rf.sr))
colnames(sC.rf.df) <- c("Sample","actual","predicted","dataset","residual")
sC.rf.df$dataset <-factor(sC.rf.df$dataset  , levels=c("Train","Test"))
sC.rf.df$ML <- "RF"
sC.rf.L <- sC.rf.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
sC.rf.L$label[1] <- paste("Train RMSE =", round(sC.rf.L$RMSE[2], 0), "<br> Test RMSE =", round(sC.rf.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(sC.rf.L$R2[1], 2), "<br> MAPE =", round(sC.rf.L$MAPE[1], 1))
sC.rf.p <- ggplot(sC.rf.df, aes(predicted, actual))+sCOD.theme("sCOD-RF",sC.rf.L)+theme(legend.position="none")
sC.rf.p

### CUB ###
sC.cub.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$sCOD, abs(predict.train(WQ.sC.cub, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
sC.cub.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$sCOD, abs(predict.train(WQ.sC.cub, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),sC.cub.df)
sC.cub.df <- sC.cub.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
sC.cub.sr <- lm(X2 ~ X3, data=sC.cub.df)
sC.cub.df <- cbind(sC.cub.df, rstandard(sC.cub.sr))
colnames(sC.cub.df) <- c("Sample","actual","predicted","dataset","residual")
sC.cub.df$dataset <-factor(sC.cub.df$dataset  , levels=c("Train","Test"))
sC.cub.df$ML <- "CUB"
sC.cub.L <- sC.cub.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
sC.cub.L$label[1] <- paste("Train RMSE =", round(sC.cub.L$RMSE[2], 0), "<br> Test RMSE =", round(sC.cub.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(sC.cub.L$R2[1], 2), "<br> MAPE =", round(sC.cub.L$MAPE[1], 1))
sC.cub.p <- ggplot(sC.cub.df, aes(predicted, actual))+sCOD.theme("sCOD-CUB",sC.cub.L)+theme(legend.position="none")
sC.cub.p

### Save sCOD Result Scatter Plot ###
leg <- get_legend(sC.pls.p)

sCOD.R.p <- ggarrange(sC.lm.p,sC.pls.p+theme(legend.position="none"),sC.knn.p,sC.svm.p,
                      sC.rf.p,sC.cub.p,sC.qrnn.p,as_ggplot(leg),
                      labels = c("A","B","C","D","E","F","G",""), ncol=4,nrow = 2, align = "v")
sCOD.R.p
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/sCOD_ScatterP.jpg",plot=sCOD.R.p,width = 12, height = 6, dpi = 300)

### sCOD Overall Scatter Plot ###
sCOD.df <- rbind(sC.qrnn.df,sC.pls.df,sC.svm.df,sC.cub.df)
sCOD.df$ML <-factor(sCOD.df$ML , levels=c("PLS","SVR","CUB","QRNN"))
sCOD.df.L <- rbind(sC.pls.L,sC.rf.L,sC.cub.L,sC.qrnn.L)
sCOD.df.L$ML <- c("PLS","PLS","SVR","SVR","CUB","CUB","QRNN","QRNN")
sCOD.df.L$label <- gsub("<br>", ";", sCOD.df.L$label)
sCOD.df.p <- ggplot(subset(sCOD.df, dataset %in% "Test"),aes(predicted, actual, fill = ML)) + geom_point(size=2.5,stroke=0.5,shape = 21,alpha = 0.8) + theme_bw() + 
  geom_abline(color = "black", linetype = 2, size = 1, alpha = 0.6)+
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]), 
                    guide = guide_legend(title = "ML Model", title.position = "top",override.aes = list(shape =  c(21,21,21,21)), order = 1))+
  scale_shape_manual(values = c(21,21,21,21),guide = F)+
  #scale_linetype_manual(values = c("longdash","longdash","longdash","longdash","longdash","longdash"), guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  #scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(1,10^3.5),
  #                   labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(1,10^3.5),
  #                   labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  labs(x = "Estimated sCOD (mg/L)", y = "Measured sCOD (mg/L)")+
  theme(text = element_text(size=15),legend.position = c(.35, 0.9), legend.direction='horizontal',
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        legend.text=element_text(size=10),legend.title=element_text(size=10))
sCOD.df.p

### COD Result Standardized Residual Plot ###
sC.sr.theme <- function(titl.,...){
  list(geom_point(aes(fill = Sample),shape=21,size=2,color="black", stroke = 0.5),theme_bw(),ylim(-3,3),xlim(-200,5000),
       scale_alpha_continuous(limits = c(0,3), breaks = c(0,1,2,3),
                              guide = guide_legend(title = "Standardized residual", title.position = "top", order = 2)),
       scale_fill_manual(values = c(colx2(5)), 
                         guide = guide_legend(title = "Sampling location", title.position = "left",override.aes = list(shape = 21), order = 1)),
       theme(text = element_text(size=7),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = c(.9, 0.3)),
       labs(title = titl. ,x = "Estimated sCOD (mg/L)", y = "Standardized Residual"))
}

sC.qrnn.p2 <- ggplot(data = subset(sC.qrnn.df, dataset %in% "Test"), aes(predicted, residual)) + sC.sr.theme("sCOD: QRNN") + guides(alpha=FALSE)
sC.pls.p2 <- ggplot(data = subset(sC.pls.df, dataset %in% "Test"), aes(predicted, residual)) + sC.sr.theme("sCOD: PLS") + guides(alpha=FALSE)
sC.knn.p2 <- ggplot(data = subset(sC.knn.df, dataset %in% "Test"), aes(predicted, residual)) + sC.sr.theme("sCOD: KNN") + guides(alpha=FALSE)
sC.rf.p2 <- ggplot(data = subset(sC.rf.df, dataset %in% "Test"), aes(predicted, residual)) + sC.sr.theme("sCOD: RF") + guides(alpha=FALSE)
sC.svm.p2 <- ggplot(data = subset(sC.svm.df, dataset %in% "Test"), aes(predicted, residual)) + sC.sr.theme("sCOD: SVR") + guides(alpha=FALSE)
sC.cub.p2 <- ggplot(data = subset(sC.cub.df, dataset %in% "Test"), aes(predicted, residual)) + sC.sr.theme("sCOD: CUB") + guides(alpha=FALSE)
sCOD.SR.p <- ggarrange(sC.pls.p2, sC.svm.p2, sC.cub.p2, sC.qrnn.p2, labels = "AUTO", ncol=4,nrow = 1, align = "v",  common.legend = TRUE,legend = "bottom")
sCOD.SR.p

sCOD.SR.p2 <- ggplot(subset(sCOD.df, dataset %in% "Test"),aes(actual, residual, fill = ML)) + geom_point(size=2.5,stroke=0.5,shape = 21,alpha = 0.8) + theme_bw() + 
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]), 
                    guide = guide_legend(title = "ML Model", title.position = "top",override.aes = list(shape = 21), order = 1))+
  #scale_shape_manual(values = c(0,2,3,4,5), guide = guide_legend(title = "Sampling Location", title.position = "top", order = 2))+
  labs(x = "Estimated sCOD (mg/L)", y = "Standardized Residual")+ylim(-6,6)+
  theme(text = element_text(size=15),legend.position = c(.35, 0.9), legend.direction='horizontal',
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        legend.text=element_text(size=8),legend.title=element_text(size=8))

sCyplot <- ggplot(subset(sCOD.df, dataset %in% "Test"), aes(x=residual, color = ML, fill = ML)) + geom_histogram(alpha = 0.5, position="dodge", binwidth = 0.5) + labs(y = "count")  + xlim(-6,6)+
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]),guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  scale_color_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]),guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  coord_flip() + theme_bw() + theme(text = element_text(size=15),legend.position="none",axis.title.y=element_blank(),
                                    axis.text.y=element_blank(),axis.ticks.y=element_blank(), plot.margin = unit(c(0.38,1,0.38,-0.5), "lines"))
sCsr.p <- ggarrange(sCOD.SR.p2, sCyplot, ncol = 2, nrow = 1, widths = c(3.5, 1), common.legend = F)
sCsr.p

#### pCOD Predict ####
set.seed(1234)
tune.method <- trainControl(method = "repeatedcv", number=15, repeats=3, selectionFunction = "best")
pc.mod.fun <- function(mod.,inputT,...){
  caret::train(pCOD ~ Group+Turb+NH4+NO3+Color+EC+pH
               , inputT, method = mod., trControl = tune.method, na.action = na.pass, preProc = c("center", "scale", "nzv"
               ),...)
}

set.seed(1234); WQ.pC.lm <- #pc.mod.fun("lm",WQ2.train)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/pCOD_EST/pCOD_lm.rds")
set.seed(1234); WQ.pC.pls <- #pc.mod.fun("pls",WQ2.train,tuneLength = 10)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/pCOD_EST/pCOD_pls.rds")
set.seed(1234); WQ.pC.rf <- #pc.mod.fun("rf",WQ2.train,ntree=25)
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/pCOD_EST/pCOD_rf.rds")
set.seed(1234); WQ.pC.knn <- #pc.mod.fun("knn",WQ2.train,tuneGrid = expand.grid (k = seq(from = 3, to = 10, by = 1)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/pCOD_EST/pCOD_knn.rds")
set.seed(1234); WQ.pC.svm <- #pc.mod.fun("svmRadial",WQ2.train,tuneGrid=expand.grid(.sigma=seq(from=0.01,to=0.1,by=0.01),.C=seq(from = 0.5,to = 5, by = 0.5)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/pCOD_EST/pCOD_svm.rds")
set.seed(1234); WQ.pC.cub <- #pc.mod.fun("cubist",WQ2.train,tuneGrid=expand.grid(.committees=seq(from=3,to=10,by=1),.neighbors = seq(from = 3, to = 9, by = 1)))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/pCOD_EST/pCOD_cub.rds")
set.seed(1234); WQ.pC.qrnn <- #pc.mod.fun("qrnn",WQ2.train,tuneGrid=expand.grid(.n.hidden=seq(from=2,to=5,by=1),.penalty = 10^(seq(from = -2, to = 2, by = 0.25)),.bag=FALSE))
  readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/pCOD_EST/pCOD_qrnn.rds")
#set.seed(1234); WQ.C.xg <- pc.mod.fun("xgbTree",WQ2.train)

#saveRDS(WQ.pC.qrnn,"G:/My Drive/R project/GitHub/MLsensor/Save_model/pCOD_EST/pCOD_qrnn.rds")

###  pCOD Train Variables Importance ### 
pC.lm.Imp <- ggplot(varImp(WQ.pC.lm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "pCOD-LM")
pC.pls.Imp <-ggplot(varImp(WQ.pC.pls, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "pCOD-PLS")
pC.rf.Imp <- ggplot(varImp(WQ.pC.rf, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "pCOD-RF")
pC.knn.Imp <- ggplot(varImp(WQ.pC.knn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "pCOD-KNN")
pC.svm.Imp <- ggplot(varImp(WQ.pC.svm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "pCOD-SVR")
pC.cub.Imp <- ggplot(varImp(WQ.pC.cub, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "pCOD-CUB")
pC.qrnn.Imp <- ggplot(varImp(WQ.pC.qrnn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "pCOD-QRNN")
pC.varImp <- ggarrange(pC.lm.Imp,pC.pls.Imp, pC.svm.Imp, pC.rf.Imp, pC.cub.Imp,pC.qrnn.Imp, ncol=3, nrow = 2, align = "v", labels = "AUTO", common.legend = TRUE, legend="bottom")
pC.varImp
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/PC_varImp.jpg",plot=pC.varImp,width = 6, height = 4, dpi = 300)

### pCOD Train RMSE Boxplot ### 
pC.results <- resamples(list("lm" = WQ.pC.lm,"qrnn" = WQ.pC.qrnn, "knn" = WQ.pC.knn, "svr" = WQ.pC.svm,"pls" = WQ.pC.pls, "rf" = WQ.pC.rf, "CUB" = WQ.pC.cub))
pC.RMSE <- data.frame(pC.results$values$`lm~RMSE`,pC.results$values$`qrnn~RMSE`,pC.results$values$`knn~RMSE`,pC.results$values$`svr~RMSE`,
                      pC.results$values$`pls~RMSE`,pC.results$values$`rf~RMSE`,pC.results$values$`CUB~RMSE`)
colnames(pC.RMSE) <- c("LM","QRNN","KNN","SVR","PLS","RF","CUB")
pC.RMSE <- pC.RMSE %>% gather(key = "Model", value = "RMSE", LM, QRNN, KNN, SVR, PLS, RF, CUB) %>% convert_as_factor(Model)
pC.Ttest <- pC.RMSE %>% pairwise_t_test(RMSE ~ Model, paired = TRUE, alt = c("two.sided"), conf.level = 0.99, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Model")
pC.RMSE.p <- ggboxplot(pC.RMSE, x = "Model", y = "RMSE",color = "black", legend = "none", order = c("LM","PLS","KNN","SVR","RF","QRNN","CUB"),outlier.shape = NA) +
  geom_point(data = pC.RMSE, aes(x = Model, y = RMSE, colour = Model, alpha = 0.9), position = position_jitter(width = 0.2))+
  geom_hline(yintercept = mean(pC.RMSE$RMSE), linetype = 2) + 
  coord_flex_cart(bottom=brackets_horizontal(), left=capped_vertical('none'))+
  theme_bw()+labs(x = element_blank(), y = "RMSE (pCOD mg/L)")+
  theme(text = element_text(size=15),legend.title=element_blank(),legend.position = "none",
        panel.border=element_blank(), axis.line = element_line(),legend.background = element_rect(colour='grey')) +
  scale_color_manual(values = c("LM"=colx(7)[1],"PLS"=colx(7)[2],"KNN"=colx(7)[3],"SVR"=colx(7)[4],"RF"=colx(7)[5],"QRNN"=colx(7)[6],"CUB"=colx(7)[7]),guide=guide_legend(title = "Algorithms", order = 1))
pC.RMSE.p
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/PC_RMSE.jpg",plot=pC.RMSE.p,width = 6, height = 4, dpi = 300)

### pCOD Test RMSE Boxplot ### 
pCOD.theme <- function(titl.,lab.data,...){
  list(geom_point(aes(fill=dataset), size=3.25, shape=21, color="black", stroke=0.5, alpha=0.75),
       theme_bw(),ylim(0,3000),xlim(0,3000),
       scale_fill_manual(values = c(colx2(2)), guide = guide_legend(title = "Dataset", title.position = "left",override.aes = list(shape=21,size=2.5), order = 1, nrow = 1)),
       theme(text = element_text(size=8.5),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = c(0.15, 0.85)),
       labs(x = "Estimated pCOD (mg/L)", y = "Measured pCOD (mg/L)",title = titl.),geom_abline(color = "black", linetype = 2, size = 1, alpha = 0.5),
       geom_richtext(data = lab.data,aes(1750, 500, label = label[1]),hjust = 0,size = 2.75,fill = "white", label.color = "black"),
       ...)
}

### LM ###
pC.lm.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$pCOD, abs(predict.train(WQ.pC.lm, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
pC.lm.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$pCOD, abs(predict.train(WQ.pC.lm, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),pC.lm.df)
pC.lm.df <- pC.lm.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
pC.lm.sr <- lm(X2 ~ X3, data=pC.lm.df)
pC.lm.df <- cbind(pC.lm.df, rstandard(pC.lm.sr))
colnames(pC.lm.df) <- c("Sample","actual","predicted","dataset","residual")
pC.lm.df$dataset <-factor(pC.lm.df$dataset  , levels=c("Train","Test"))
pC.lm.df$ML <- "LM"
pC.lm.L <- pC.lm.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
pC.lm.L$label[1] <- paste("Train RMSE =", round(pC.lm.L$RMSE[2], 0), "<br> Test RMSE =", round(pC.lm.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(pC.lm.L$R2[1], 2), "<br> MAPE =", round(pC.lm.L$MAPE[1], 1))
pC.lm.p <- ggplot(pC.lm.df, aes(predicted, actual))+pCOD.theme("pCOD-LM",pC.lm.L)+theme(legend.position="none")
pC.lm.p

### QRNN ###
pC.qrnn.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$pCOD, abs(predict.train(WQ.pC.qrnn, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
pC.qrnn.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$pCOD, abs(predict.train(WQ.pC.qrnn, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),pC.qrnn.df)
pC.qrnn.df <- pC.qrnn.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
pC.qrnn.sr <- lm(X2 ~ X3, data=pC.qrnn.df)
pC.qrnn.df <- cbind(pC.qrnn.df, rstandard(pC.qrnn.sr))
colnames(pC.qrnn.df) <- c("Sample","actual","predicted","dataset","residual")
pC.qrnn.df$dataset <-factor(pC.qrnn.df$dataset  , levels=c("Train","Test"))
pC.qrnn.df$ML <- "QRNN"
pC.qrnn.L <- pC.qrnn.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
pC.qrnn.L$label[1] <- paste("Train RMSE =", round(pC.qrnn.L$RMSE[2], 0), "<br> Test RMSE =", round(pC.qrnn.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(pC.qrnn.L$R2[1], 2), "<br> MAPE =", round(pC.qrnn.L$MAPE[1], 1))
pC.qrnn.p <- ggplot(pC.qrnn.df, aes(predicted, actual))+pCOD.theme("pCOD-QRNN",pC.qrnn.L)+theme(legend.position="none")
pC.qrnn.p

### PLS ###
pC.pls.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$pCOD, abs(predict.train(WQ.pC.pls, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
pC.pls.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$pCOD, abs(predict.train(WQ.pC.pls, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),pC.pls.df)
pC.pls.df <- pC.pls.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
pC.pls.sr <- lm(X2 ~ X3, data=pC.pls.df)
pC.pls.df <- cbind(pC.pls.df, rstandard(pC.pls.sr))
colnames(pC.pls.df) <- c("Sample","actual","predicted","dataset","residual")
pC.pls.df$dataset <-factor(pC.pls.df$dataset  , levels=c("Train","Test"))
pC.pls.df$ML <- "PLS"
pC.pls.L <- pC.pls.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
pC.pls.L$label[1] <- paste("Train RMSE =", round(pC.pls.L$RMSE[2], 0), "<br> Test RMSE =", round(pC.pls.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(pC.pls.L$R2[1], 2), "<br> MAPE =", round(pC.pls.L$MAPE[1], 1))
pC.pls.p <- ggplot(pC.pls.df, aes(predicted, actual))+pCOD.theme("pCOD-PLS",pC.pls.L)
pC.pls.p

### KNN ###
pC.knn.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$pCOD, abs(predict.train(WQ.pC.knn, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
pC.knn.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$pCOD, abs(predict.train(WQ.pC.knn, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),pC.knn.df)
pC.knn.df <- pC.knn.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
pC.knn.sr <- lm(X2 ~ X3, data=pC.knn.df)
pC.knn.df <- cbind(pC.knn.df, rstandard(pC.knn.sr))
colnames(pC.knn.df) <- c("Sample","actual","predicted","dataset","residual")
pC.knn.df$dataset <-factor(pC.knn.df$dataset  , levels=c("Train","Test"))
pC.knn.df$ML <- "KNN"
pC.knn.L <- pC.knn.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
pC.knn.L$label[1] <- paste("Train RMSE =", round(pC.knn.L$RMSE[2], 0), "<br> Test RMSE =", round(pC.knn.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(pC.knn.L$R2[1], 2), "<br> MAPE =", round(pC.knn.L$MAPE[1], 1))
pC.knn.p <- ggplot(pC.knn.df, aes(predicted, actual))+pCOD.theme("pCOD-KNN",pC.knn.L)+theme(legend.position="none")
pC.knn.p

### SVM ###
pC.svm.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$pCOD, abs(predict.train(WQ.pC.svm, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
pC.svm.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$pCOD, abs(predict.train(WQ.pC.svm, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),pC.svm.df)
pC.svm.df <- pC.svm.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
pC.svm.sr <- lm(X2 ~ X3, data=pC.svm.df)
pC.svm.df <- cbind(pC.svm.df, rstandard(pC.svm.sr))
colnames(pC.svm.df) <- c("Sample","actual","predicted","dataset","residual")
pC.svm.df$dataset <-factor(pC.svm.df$dataset  , levels=c("Train","Test"))
pC.svm.df$ML <- "SVR"
pC.svm.L <- pC.svm.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
pC.svm.L$label[1] <- paste("Train RMSE =", round(pC.svm.L$RMSE[2], 0), "<br> Test RMSE =", round(pC.svm.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(pC.svm.L$R2[1], 2), "<br> MAPE =", round(pC.svm.L$MAPE[1], 1))
pC.svm.p <- ggplot(pC.svm.df, aes(predicted, actual))+pCOD.theme("pCOD-SVR",pC.svm.L)+theme(legend.position="none")
pC.svm.p

### RF ###
pC.rf.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$pCOD, abs(predict.train(WQ.pC.rf, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
pC.rf.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$pCOD, abs(predict.train(WQ.pC.rf, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),pC.rf.df)
pC.rf.df <- pC.rf.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
pC.rf.sr <- lm(X2 ~ X3, data=pC.rf.df)
pC.rf.df <- cbind(pC.rf.df, rstandard(pC.rf.sr))
colnames(pC.rf.df) <- c("Sample","actual","predicted","dataset","residual")
pC.rf.df$dataset <-factor(pC.rf.df$dataset  , levels=c("Train","Test"))
pC.rf.df$ML <- "RF"
pC.rf.L <- pC.rf.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
pC.rf.L$label[1] <- paste("Train RMSE =", round(pC.rf.L$RMSE[2], 0), "<br> Test RMSE =", round(pC.rf.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(pC.rf.L$R2[1], 2), "<br> MAPE =", round(pC.rf.L$MAPE[1], 1))
pC.rf.p <- ggplot(pC.rf.df, aes(predicted, actual))+pCOD.theme("pCOD-RF",pC.rf.L)+theme(legend.position="none")
pC.rf.p

### CUB ###
pC.cub.df <- data.frame(cbind(WQ2.test$Group,WQ2.test$pCOD, abs(predict.train(WQ.pC.cub, newdata = WQ2.test[c(1,7:12)], type="raw"))),"ds"="Test")
pC.cub.df <- rbind(data.frame(cbind(WQ2.train$Group,WQ2.train$pCOD, abs(predict.train(WQ.pC.cub, newdata = WQ2.train[c(1,7:12)], type="raw"))),"ds"="Train"),pC.cub.df)
pC.cub.df <- pC.cub.df %>% within(X1 <- factor(X1, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
pC.cub.sr <- lm(X2 ~ X3, data=pC.cub.df)
pC.cub.df <- cbind(pC.cub.df, rstandard(pC.cub.sr))
colnames(pC.cub.df) <- c("Sample","actual","predicted","dataset","residual")
pC.cub.df$dataset <-factor(pC.cub.df$dataset  , levels=c("Train","Test"))
pC.cub.df$ML <- "CUB"
pC.cub.L <- pC.cub.df %>% group_by(dataset) %>%
  summarise(RMSE = caret::RMSE(predicted,actual),
            R2 = caret::R2(predicted,actual),
            MAPE = MAPE(predicted,actual)) %>% ungroup()
pC.cub.L$label[1] <- paste("Train RMSE =", round(pC.cub.L$RMSE[2], 0), "<br> Test RMSE =", round(pC.cub.L$RMSE[1], 0), "<br> R<sup>2</sup> =", round(pC.cub.L$R2[1], 2), "<br> MAPE =", round(pC.cub.L$MAPE[1], 1))
pC.cub.p <- ggplot(pC.cub.df, aes(predicted, actual))+pCOD.theme("pCOD-CUB",pC.cub.L)+theme(legend.position="none")
pC.cub.p

### Save pCOD Result Scatter Plot ###
pCOD.R.p <- ggarrange(pC.lm.p,pC.pls.p+theme(legend.position="none"),pC.knn.p,pC.svm.p,
                      pC.rf.p,pC.cub.p,pC.qrnn.p,as_ggplot(leg),
                      labels = c("A","B","C","D","E","F","G",""), ncol=4,nrow = 2, align = "v")
pCOD.R.p
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/pCOD_ScatterP.jpg",plot=pCOD.R.p,width = 12, height = 6, dpi = 300)

### pCOD Overall Scatter Plot ###
pCOD.df <- rbind(pC.qrnn.df,pC.pls.df,pC.svm.df,pC.cub.df)
pCOD.df$ML <-factor(pCOD.df$ML , levels=c("PLS","SVR","CUB","QRNN"))
pCOD.df.L <- rbind(pC.pls.L,pC.rf.L,pC.cub.L,pC.qrnn.L)
pCOD.df.L$ML <- c("PLS","PLS","SVR","SVR","CUB","CUB","QRNN","QRNN")
pCOD.df.L$label <- gsub("<br>", ";", pCOD.df.L$label)
pCOD.df.p <- ggplot(subset(pCOD.df, dataset %in% "Test"),aes(predicted, actual, fill = ML)) + geom_point(size=2.5,stroke=0.5,shape = 21,alpha = 0.8) + theme_bw() + 
  geom_abline(color = "black", linetype = 2, size = 1, alpha = 0.6)+
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]), 
                    guide = guide_legend(title = "ML Model", title.position = "top",override.aes = list(shape =  c(21,21,21,21)), order = 1))+
  scale_shape_manual(values = c(21,21,21,21),guide = F)+
  #scale_linetype_manual(values = c("longdash","longdash","longdash","longdash","longdash","longdash"), guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  #scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(1,10^3.5),
  #                   labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),limits=c(1,10^3.5),
  #                   labels = scales::trans_format("log10", scales::math_format (10^.x)))+
  labs(x = "Estimated pCOD (mg/L)", y = "Measured pCOD (mg/L)")+
  theme(text = element_text(size=15),legend.position = c(.35, 0.9), legend.direction='horizontal',
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        legend.text=element_text(size=10),legend.title=element_text(size=10))
pCOD.df.p

### COD Result Standardized Residual Plot ###
pC.sr.theme <- function(titl.,...){
  list(geom_point(aes(fill = Sample),shape=21,size=2,color="black", stroke = 0.5),theme_bw(),ylim(-3,3),xlim(-200,5000),
       scale_alpha_continuous(limits = c(0,3), breaks = c(0,1,2,3),
                              guide = guide_legend(title = "Standardized residual", title.position = "top", order = 2)),
       scale_fill_manual(values = c(colx2(5)), 
                         guide = guide_legend(title = "Sampling location", title.position = "left",override.aes = list(shape = 21), order = 1)),
       theme(text = element_text(size=7),legend.background = element_blank(),legend.box.background = element_blank(),legend.position = c(.9, 0.3)),
       labs(title = titl. ,x = "Estimated sCOD (mg/L)", y = "Standardized Residual"))
}

pC.qrnn.p2 <- ggplot(data = subset(pC.qrnn.df, dataset %in% "Test"), aes(predicted, residual)) + pC.sr.theme("sCOD: QRNN") + guides(alpha=FALSE)
pC.pls.p2 <- ggplot(data = subset(pC.pls.df, dataset %in% "Test"), aes(predicted, residual)) + pC.sr.theme("sCOD: PLS") + guides(alpha=FALSE)
pC.knn.p2 <- ggplot(data = subset(pC.knn.df, dataset %in% "Test"), aes(predicted, residual)) + pC.sr.theme("sCOD: KNN") + guides(alpha=FALSE)
pC.rf.p2 <- ggplot(data = subset(pC.rf.df, dataset %in% "Test"), aes(predicted, residual)) + pC.sr.theme("sCOD: RF") + guides(alpha=FALSE)
pC.svm.p2 <- ggplot(data = subset(pC.svm.df, dataset %in% "Test"), aes(predicted, residual)) + pC.sr.theme("sCOD: SVR") + guides(alpha=FALSE)
pC.cub.p2 <- ggplot(data = subset(pC.cub.df, dataset %in% "Test"), aes(predicted, residual)) + pC.sr.theme("sCOD: CUB") + guides(alpha=FALSE)
pCOD.SR.p <- ggarrange(pC.pls.p2, pC.svm.p2, pC.cub.p2, pC.qrnn.p2, labels = "AUTO", ncol=4,nrow = 1, align = "v",  common.legend = TRUE,legend = "bottom")
pCOD.SR.p

pCOD.SR.p2 <- ggplot(subset(pCOD.df, dataset %in% "Test"),aes(actual, residual, fill = ML)) + geom_point(size=2.5,stroke=0.5,shape = 21,alpha = 0.8) + theme_bw() + 
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]), 
                    guide = guide_legend(title = "ML Model", title.position = "top",override.aes = list(shape = 21), order = 1))+
  #scale_shape_manual(values = c(0,2,3,4,5), guide = guide_legend(title = "Sampling Location", title.position = "top", order = 2))+
  labs(x = "Estimated sCOD (mg/L)", y = "Standardized Residual")+ylim(-6,6)+
  theme(text = element_text(size=15),legend.position = c(.35, 0.9), legend.direction='horizontal',
        legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"),
        legend.text=element_text(size=8),legend.title=element_text(size=8))

pCyplot <- ggplot(subset(pCOD.df, dataset %in% "Test"), aes(x=residual, color = ML, fill = ML)) + geom_histogram(alpha = 0.5, position="dodge", binwidth = 0.5) + labs(y = "count")  + xlim(-6,6)+
  scale_fill_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]),guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  scale_color_manual(values = c("PLS"=colx(4)[1],"SVR"=colx(4)[2],"CUB"=colx(4)[3],"QRNN"=colx(4)[4]),guide = guide_legend(title = "ML Model", title.position = "top", order = 1))+
  coord_flip() + theme_bw() + theme(text = element_text(size=15),legend.position="none",axis.title.y=element_blank(),
                                    axis.text.y=element_blank(),axis.ticks.y=element_blank(), plot.margin = unit(c(0.38,1,0.38,-0.5), "lines"))
pCsr.p <- ggarrange(pCOD.SR.p2, pCyplot, ncol = 2, nrow = 1, widths = c(3.5, 1), common.legend = F)
pCsr.p


#### COD Predict ####
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


###  COD Train Variables Importance ### 
C.lm.Imp <- ggplot(varImp(WQ.C.lm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-LM")
C.pls.Imp <-ggplot(varImp(WQ.C.pls, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-PLS")
C.rf.Imp <- ggplot(varImp(WQ.C.rf, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-RF")
C.knn.Imp <- ggplot(varImp(WQ.C.knn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-KNN")
C.svm.Imp <- ggplot(varImp(WQ.C.svm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-SVR")
C.cub.Imp <- ggplot(varImp(WQ.C.cub, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-CUB")
C.qrnn.Imp <- ggplot(varImp(WQ.C.qrnn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "COD-QRNN")
C.varImp <- ggarrange(C.lm.Imp,C.pls.Imp, C.svm.Imp, C.rf.Imp, C.cub.Imp,C.qrnn.Imp, ncol=3, nrow = 2, align = "v", labels = "AUTO", common.legend = TRUE, legend="bottom")
C.varImp
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/C_varImp.jpg",plot=C.varImp,width = 6, height = 4, dpi = 300)

### COD Train RMSE Boxplot ### 
C.results <- resamples(list("lm" = WQ.C.lm, "qrnn" = WQ.C.qrnn, "knn" = WQ.C.knn, "svr" = WQ.C.svm,"pls" = WQ.C.pls, "rf" = WQ.C.rf, "CUB" = WQ.C.cub))
C.RMSE <- data.frame(C.results$values$`lm~RMSE`,C.results$values$`qrnn~RMSE`,C.results$values$`knn~RMSE`,C.results$values$`svr~RMSE`,
                     C.results$values$`pls~RMSE`,C.results$values$`rf~RMSE`,C.results$values$`CUB~RMSE`)
colnames(C.RMSE) <- c("LM","QRNN","KNN","SVR","PLS","RF","CUB")
C.RMSE <- C.RMSE %>% gather(key = "Model", value = "RMSE", LM,QRNN, KNN, SVR, PLS, RF, CUB) %>% convert_as_factor(Model)
C.Ttest <- C.RMSE %>% pairwise_t_test(RMSE ~ Model, paired = TRUE, alt = c("two.sided"), conf.level = 0.99, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Model")
C.RMSE.p <- ggboxplot(C.RMSE, x = "Model", y = "RMSE",color = "black", legend = "none", order = c("LM","PLS","KNN","SVR","RF","QRNN","CUB"),outlier.shape = NA) +
  geom_point(data = C.RMSE, aes(x = Model, y = RMSE, colour = Model, alpha = 0.9), position = position_jitter(width = 0.2))+
  geom_hline(yintercept = mean(C.RMSE$RMSE), linetype = 2) + 
  coord_flex_cart(bottom=brackets_horizontal(), left=capped_vertical('none'))+
  theme_bw()+labs(x = element_blank(), y = "RMSE (COD mg/L)")+
  theme(text = element_text(size=15),legend.title=element_blank(),legend.position = "none",
        panel.border=element_blank(), axis.line = element_line(),legend.background = element_rect(colour='grey')) +
  scale_color_manual(values = c("LM"=colx(7)[1],"PLS"=colx(7)[2],"KNN"=colx(7)[3],"SVR"=colx(7)[4],"RF"=colx(7)[5],"QRNN"=colx(7)[6],"CUB"=colx(7)[7]),guide=guide_legend(title = "Algorithms", order = 1))
C.RMSE.p
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/C_RMSE.jpg",plot=C.RMSE.p,width = 6, height = 4, dpi = 300)

### COD Test RMSE Boxplot ### 
cod.theme <- function(titl.,lab.data,...){
  list(geom_point(aes(fill=dataset), size=3.25, shape=21, color="black", stroke=0.5, alpha=0.75),
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
C.lm.p <- ggplot(C.lm.df, aes(predicted, actual))+cod.theme("COD-LM",C.lm.L)+theme(legend.position="none")
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
C.qrnn.p <- ggplot(C.qrnn.df, aes(predicted, actual))+cod.theme("COD-QRNN",C.qrnn.L)+theme(legend.position="none")
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
C.pls.p <- ggplot(C.pls.df, aes(predicted, actual))+cod.theme("COD-PLS",C.pls.L)
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
C.knn.p <- ggplot(C.knn.df, aes(predicted, actual))+cod.theme("COD-KNN",C.knn.L)+theme(legend.position="none")
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
C.svm.p <- ggplot(C.svm.df, aes(predicted, actual))+cod.theme("COD-SVR",C.svm.L)+theme(legend.position="none")
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
C.rf.p <- ggplot(C.rf.df, aes(predicted, actual))+cod.theme("COD-RF",C.rf.L)+theme(legend.position="none")
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
C.cub.p <- ggplot(C.cub.df, aes(predicted, actual))+cod.theme("COD-CUB",C.cub.L)+theme(legend.position="none")
C.cub.p

### Save COD Result Scatter Plot ###
COD.R.p <- ggarrange(C.lm.p,C.pls.p+theme(legend.position="none"),C.knn.p,C.svm.p,
                     C.rf.p,C.cub.p,C.qrnn.p,as_ggplot(leg),
                      labels = c("A","B","C","D","E","F","G",""), ncol=4,nrow = 2, align = "v")
COD.R.p
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/COD_ScatterP.jpg",plot=COD.R.p,width = 12, height = 6, dpi = 300)

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