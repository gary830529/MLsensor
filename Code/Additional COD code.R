pacman::p_load(plyr,dplyr,reshape2,stringr,lubridate,readxl,openxlsx,janitor,ggplot2,ggpmisc,ggpubr,ggsignif,gplots,data.table,magrittr,grid,
               tidyverse,rstatix,PtProcess,gtable,gridExtra,cowplot,tiff,pwr,MOTE,kableExtra,tinytex,knitr,quantmod,class,caret,gmodels,
               officedown,officer,flextable,car,C50,GGally,ggResidpanel,ggfortify,caretEnsemble,randomForest,corrplot,neuralnet,lemon,Hmisc,scales,
               glue,ggtext,png,gtools,ggrepel,rvg,gdata,scales,nnet,xgboost,correlation,psych,corrgram,jpeg,qrnn,pls,kernlab)
colx <- colorRampPalette(c("#16068a","#9e189d","#fdb32e"))
colx2 <- colorRampPalette(c("#264653","#2a9d8e","#e9c46b","#f3a261","#e66f51"))

#### Start ####
Raw <- read.csv(file = 'G:/My Drive/R project/GitHub/MLsensor/Raw_data/Data_ML.csv',header = T, sep = ",")
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

Dates_after_events2 <- as.Date(c("2019-01-15", "2019-01-22"))
WQ2 <- WQ %>% filter(Date %nin% Dates_after_events2)
WQ2 <- WQ2[complete.cases(WQ2), ]
WQ2 <- WQ2 %>% subset(select=-c(Date,Event,BP,VSS))
WQ2$Group <- as.numeric(as.factor(WQ2$Group))


set.seed(1234)
WQ2.index <- sample(1:nrow(WQ2),size=nrow(WQ2)*0.7,replace = FALSE)
WQ2.train <- WQ2[WQ2.index,]
WQ2.test <- WQ2[-WQ2.index,]


#### Correlation for COD PCOD SCOD ####
WQ2.cor <- WQ2 %>% within(Group <- factor(Group, labels = c("Influent","AnMBR","Permeate","Post-NCS","Effluent")))
colnames(WQ2.cor) <- c("Group", "COD", "sCOD","pCOD","TSS","E.coli","Color","Trubidity","EC","pH","NH4","NO3","Temperature","E.coli2") 


WQ2.Cor <- ggpairs(subset(WQ2.cor, select = -c(Group,E.coli2)), 
                   upper = list(continuous = upper_fn), lower = list(continuous = lower_fn), diag = list(continuous = diag_fn))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.text.y = element_text(size=5),text = element_text(size=6))
WQ2.Cor

#### Correlation between COD PCOD SCOD ####
WQ2.COD.cor <- ggpairs(WQ2.cor[1:4],aes(colour = Group, alpha = 0.4), title="COD Correlations",
                      upper = list(continuous = upper_fn2), lower = list(continuous = lower_fn2), diag = list(continuous = diag_fn))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5), axis.text.y = element_text(size=5),panel.grid.major = element_blank())+
  scale_fill_manual(values = c(colx2(5)))+
  scale_color_manual(values = c(colx2(5)))
WQ2.COD.cor

#### SCOD Predict ####
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
sC.lm.Imp <- ggplot(varImp(WQ.sC.lm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "SCOD-LM")
sC.pls.Imp <-ggplot(varImp(WQ.sC.pls, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "SCOD-PLS")
sC.rf.Imp <- ggplot(varImp(WQ.sC.rf, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "SCOD-RF")
sC.knn.Imp <- ggplot(varImp(WQ.sC.knn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "SCOD-KNN")
sC.svm.Imp <- ggplot(varImp(WQ.sC.svm, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "SCOD-SVR")
sC.cub.Imp <- ggplot(varImp(WQ.sC.cub, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "SCOD-CUB")
sC.qrnn.Imp <- ggplot(varImp(WQ.sC.qrnn, scale = T))+theme_bw()+labs(x=element_blank(), y="Importance (%)", title = "SCOD-QRNN")
sC.varImp <- ggarrange(sC.pls.Imp, sC.svm.Imp, sC.cub.Imp,sC.qrnn.Imp, ncol=2, nrow = 2, align = "v", labels = "AUTO", common.legend = TRUE, legend="bottom")
sC.varImp
#ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/SC_varImp.jpg",plot=sC.varImp,width = 6, height = 4, dpi = 300)

### sCOD Train RMSE Boxplot ### 
sC.results <- resamples(list("qrnn" = WQ.sC.qrnn, "knn" = WQ.sC.knn, "svr" = WQ.sC.svm,"pls" = WQ.sC.pls, "rf" = WQ.sC.rf, "CUB" = WQ.sC.cub))
sC.RMSE <- data.frame(sC.results$values$`qrnn~RMSE`,sC.results$values$`knn~RMSE`,sC.results$values$`svr~RMSE`,
                     sC.results$values$`pls~RMSE`,sC.results$values$`rf~RMSE`,sC.results$values$`CUB~RMSE`)
colnames(sC.RMSE) <- c("QRNN","KNN","SVR","PLS","RF","CUB")
sC.RMSE <- sC.RMSE %>% gather(key = "Model", value = "RMSE", QRNN, KNN, SVR, PLS, RF, CUB) %>% convert_as_factor(Model)
sC.Ttest <- sC.RMSE %>% pairwise_t_test(RMSE ~ Model, paired = TRUE, alt = c("two.sided"), conf.level = 0.99, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Model")
sC.RMSE.p <- ggboxplot(sC.RMSE, x = "Model", y = "RMSE",color = "black", legend = "none", order = c("PLS","KNN","SVR","RF","QRNN","CUB"),outlier.shape = NA) +
  geom_point(data = sC.RMSE, aes(x = Model, y = RMSE, colour = Model, alpha = 0.9), position = position_jitter(width = 0.2))+
  geom_hline(yintercept = mean(sC.RMSE$RMSE), linetype = 2) + 
  coord_flex_cart(bottom=brackets_horizontal(), left=capped_vertical('none'))+
  theme_bw()+labs(x = element_blank(), y = "RMSE (COD mg/L)")+
  theme(text = element_text(size=15),legend.title=element_blank(),legend.position = "none",
        panel.border=element_blank(), axis.line = element_line(),legend.background = element_rect(colour='grey')) +
  scale_color_manual(values = c("PLS"=colx(6)[1],"KNN"=colx(6)[2],"SVR"=colx(6)[3],"RF"=colx(6)[4],"QRNN"=colx(6)[5],"CUB"=colx(6)[6]),guide=guide_legend(title = "Algorithms", order = 1))
sC.RMSE.p

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
ggsave("G:/My Drive/R project/GitHub/MLsensor/Figure/Figure5-1.jpg",plot=COD.R.p,width = 12, height = 2, dpi = 300)

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

#### PCOD Predict ####
set.seed(1234)
tune.method <- trainControl(method = "repeatedcv", number=15, repeats=3, selectionFunction = "best")
pc.mod.fun <- function(mod.,inputT,...){
  caret::train(sCOD ~ Group+Turb+NH4+NO3+Color+EC+pH
               , inputT, method = mod., trControl = tune.method, na.action = na.pass, preProc = c("center", "scale", "nzv"
               ),...)
}

set.seed(1234); WQ.sC.lm <- #sc.mod.fun("lm",WQ2.train)
  #readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_lm.rds")
set.seed(1234); WQ.sC.pls <- #sc.mod.fun("pls",WQ2.train,tuneLength = 10)
  #readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_pls.rds")
set.seed(1234); WQ.sC.rf <- #sc.mod.fun("rf",WQ2.train,ntree=25)
  #readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_rf.rds")
set.seed(1234); WQ.sC.knn <- #sc.mod.fun("knn",WQ2.train,tuneGrid = expand.grid (k = seq(from = 3, to = 10, by = 1)))
  #readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_knn.rds")
set.seed(1234); WQ.sC.svm <- #sc.mod.fun("svmRadial",WQ2.train,tuneGrid=expand.grid(.sigma=seq(from=0.01,to=0.1,by=0.01),.C=seq(from = 25,to = 35, by = 0.5)))
  #readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_svm.rds")
set.seed(1234); WQ.sC.cub <- #sc.mod.fun("cubist",WQ2.train,tuneGrid=expand.grid(.committees=seq(from=3,to=10,by=1),.neighbors = seq(from = 3, to = 9, by = 1)))
  #readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_cub.rds")
set.seed(1234); WQ.sC.qrnn <- #sc.mod.fun("qrnn",WQ2.train,tuneGrid=expand.grid(.n.hidden=seq(from=2,to=5,by=1),.penalty = 10^(seq(from = -2, to = 2, by = 0.25)),.bag=FALSE))
  #readRDS("G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_qrnn.rds")
#set.seed(1234); WQ.C.xg <- sc.mod.fun("xgbTree",WQ2.train)

#saveRDS(WQ.sC.svm,"G:/My Drive/R project/GitHub/MLsensor/Save_model/sCOD_EST/sCOD_svm.rds")
