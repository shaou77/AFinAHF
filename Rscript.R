library(readr)
library(tidyverse)
library(tableone)
library(pROC)
library(rms)
library(rmda)
library(glmnet)
library(PredictABEL)
library(ResourceSelection)
library(regplot)
library(DynNom)
library(mice)
library(caret)
library(forestplot)
library(ggDCA)

data <- read_csv("data.csv", col_types = cols(af = col_factor(levels = c("0",  "1")), hospital_expire_flag = col_factor(levels = c("0", "1"))))
data <- data[,-c(39:46)]
data <- data %>% replace_na(list(vasoactive=0,dialysis_active=0,invasive=0,noninvasive=0,o2=0,hfnc=0)) %>% mutate(uokg=uo/weight/24)
data$ethnicity <- as.factor(data$ethnicity)
data$mechvent <- as.factor(data$mechvent)
data$vasoactive <- as.factor(data$vasoactive)
data$dialysis_active <- as.factor(data$dialysis_active)
data <- data[,-c(1,2,38,45:48)]
data1 <- data[,-c(2,29,39,40)]
data1t <- data1 %>% filter(los_icu>=2)
missing.percent <- unlist(lapply(data1t, function(x) sum(is.na(x))))/nrow(data1)
missing.percent
vars <- c('gender','age','ethnicity','los_hospital','hospital_expire_flag','los_icu',
          "charlson_comorbidity_index",'sofa','gcs_min','apsiii','hematocrit','mbp','hr','t','spo2','r','wbc','hb','alb','ag','hco3',
          'bun','cr','ca','gs','na','k','pt','ptt','alt','ast','tnt','tb','ntprobnp','dialysis_active','uokg','vasoactive','mechvent')
table1 <- CreateTableOne(vars = vars, strata = c("af"), data = data1t, 
                         testExact = fisher.test, testNonNormal = kruskal.test)
p1<-print(table1, nonnormal=c('age','los_hospital','hospital_expire_flag','los_icu',
                              "charlson_comorbidity_index",'sofa','gcs_min','aspiii','mbp','hr','t','spo2','r','wbc','hb','hematocrit','alb','ag','hco3',
                              'bun','cr','ca','gs','na','k','pt','aptt','alt','ast','tnt','tb','ntprobnp','uokg'),
          argsExact=c('gender','ethnicity','vasoactive','dialysis_active','mechvent'))
tableall <- CreateTableOne(vars = vars, data = data1t)
p2<-print(tableall,nonnormal=c('age','los_hospital','hospital_expire_flag','los_icu',
                               "charlson_comorbidity_index",'sofa','gcs_min','aspiii','mbp','hr','t','spo2','r','wbc','hb','hematocrit','alb','ag','hco3',
                               'bun','cr','ca','gs','na','k','pt','aptt','alt','ast','tnt','tb','ntprobnp','uokg'),
          argsExact=c('gender','ethnicity','vasoactive','dialysis_active','mechvent'))
write.csv(p1,'table1a.csv')
write.csv(p2,'table1b.csv')
data2 <- data1t[,-c(4:7,12,23:25)]
data2 <- data2 %>% replace_na(list(tnt=-1,ntprobnp=-1)) %>% mutate(ntprobnp=ifelse(ntprobnp==-1,0,1),tnt=ifelse(tnt==-1,0,1))
data2$ntprobnp <- as.factor(data2$ntprobnp)
data2$tnt <- as.factor(data2$tnt)
missing.percent2 <- unlist(lapply(data2, function(x) sum(is.na(x))))/nrow(data2)
tempData <- mice(data2,m=5,maxit=5,seed=1024)
data2_imp <- complete(tempData,action = 5)
a <- createDataPartition(data2_imp$af, p = 0.7, list=FALSE)
data_train <- data2_imp[a,]
data_train2 <- data2_imp[a,]
data_test <- data2_imp[-a,]
data_train <- data_train%>%mutate(sex=ifelse(gender=='M',1,0)) 
data_train$sex <- as.numeric(data_train$sex)-1
data_train$ntprobnp <- as.numeric(data_train$ntprobnp)-1
data_train$tnt <- as.numeric(data_train$tnt)-1
data_train$vasoactive <- as.numeric(data_train$vasoactive)-1
data_train$mechvent <- as.numeric(data_train$mechvent)-1
data_train$dialysis_active <- as.numeric(data_train$dialysis_active)-1
x1 <- data_train[,-c(1,2,4)]
y1 <- data_train$af
x1 <- as.matrix(x1)
myfit <- glmnet(x1, y1, family = "binomial",alpha =1 )
plot(myfit, xvar = "lambda", label = TRUE)
myfit2 <- cv.glmnet(x1, y1, family="binomial")
plot(myfit2)
abline(v=log(c(myfit2$lambda.min,myfit2$lambda.1se)),lty="dashed")
myfit2$lambda.1se
coe <- coef(myfit, s = myfit2$lambda.1se)
act_index <- which(coe != 0)
act_coe <- coe[act_index]
row.names(coe)[act_index]
attach(data_train)
dd <- datadist(data_train2)
options(datadist='dd')
fit1 <- lrm(af~age+pt+sofa+vasoactive+apsiii+hr,data = data_train2,x=T,y=T)
fit1
nom1 <- nomogram(fit1,fun=plogis,fun.at=c(0.001,.01,.05,seq(.1,.9,by=.1),.95,.99,.999),lp=F,funlabel='The rate of atrial fibrillation')
plot(nom1)
regplot(fit1,   showP = T,  droplines = F,   rank="sd",    interval="confidence",    title='Nomogram for predicting new-onset AF in patient with AHF',
        points = T) 
cal <- calibrate(fit1,method = 'boot',B=1000)
plot(cal,xlim=c(0,1),ylim=c(0,1))
detach(data_train)
mod <- glm(af~age+pt+sofa+vasoactive+apsiii+hr,data=data_train2,family = 'binomial')
mod2 <- lrm(af~age+pt+sofa+vasoactive+apsiii+hr,data=data_train2,x=T,y=T)
hoslem.test(mod$y,fitted(mod),g=10)
pr.e <- predict(mod,data_test,type = c('response'))
hoslem.test(as.integer(data_test$af)-1,pr.e,g=10)
data_test2 <- data_test
data_test2$af <- as.numeric(data_test2$af)-1
plotCalibration(data = data_test2,cOutcome = 1, predRisk = pr.e, groups = 10, rangeaxis = c(0,1))
fit2 <- lrm(af~age+pt+sofa+vasoactive+apsiii+hr,data = data_test,x=T,y=T)
cal2 <- calibrate(fit2,method = 'boot',B=1000)
plot(cal2,xlim=c(0,1),ylim=c(0,1))

roccurve <- roc(af~pr.e,data=data_test2)
plot.roc(roccurve)
auc(roccurve)
ci.auc(roccurve)
data_test3 <- data_test
data_test3$vasoactive <- as.numeric(data_test3$vasoactive)-1
data_test3$mechvent <-  as.numeric(data_test3$mechvent)-1
data_test3$af <- as.numeric(data_test3$af)-1

dca_model1 <- decision_curve( af~age+pt+sofa+vasoactive+apsiii+hr,data=data_test3,family = binomial(link = logit),thresholds = seq(0.1, by =0.01),
                             confidence.intervals=0.95,study.design='case-control',population.prevalence=0.22)
dca_model2 <- decision_curve( af~apsiii,data=data_test3,family = binomial(link = logit),thresholds = seq(0.1, by =0.01),
                              confidence.intervals=0.95,study.design='case-control',population.prevalence=0.22)
dca_model3 <- decision_curve( af~sofa,data=data_test3,family = binomial(link = logit),thresholds = seq(0.1, by =0.01),
                              confidence.intervals=0.95,study.design='case-control',population.prevalence=0.22)
model_all <- list(dca_model1,dca_model2,dca_model3)
plot_decision_curve(model_all,curve.names =c('Model','APSIII','SOFA'),cost.benefit.axis = T, col=c('red','blue','brown'),confidence.intervals = FALSE,standardize = T)



DNbuilder(mod) 

rs_forest <- read.csv('or.csv',header = FALSE)
forestplot(labeltext = as.matrix(rs_forest[,1:2]), title="Odds Ratio",
           mean = rs_forest$V3,
           lower = rs_forest$V4, 
           upper = rs_forest$V5, 
           is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F),
           zero = 1, 
           boxsize = 0.4, 
           lineheight = unit(8,'mm'),
           colgap = unit(2,'mm'),
           lwd.zero = 2,
           lwd.ci = 2,
           col=fpColors(box='#458B00',summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           xlab="The estimates",
           lwd.xaxis=2,
           lty.ci = "solid",
           graph.pos = 2)
