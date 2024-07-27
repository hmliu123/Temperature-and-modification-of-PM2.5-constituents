###########################################################################################################################################################
#Codes for "Non-optimum ambient temperature and cause-specific hospital admissions for cardiovascular diseases in the context of particulate air pollution"
#Authors for codes: Huimeng Liu, Jian Lei, Yunxing Jiang, Lijun Bai et.al.
#Correspondence to Shaowei Wu, Xi Chen.
##########################################################################################################################################################

######################################################################################################################################
#MODIFICATION BY PM2.5 AND ITS CONSTITUENTS ON THE ASSOCIATIONS BETWEEN NON-OPTIMUM TEMPERATURE AND HOSPITAL ADMISSIONS FOR MAJOR CVDS
##########################################################################################################################################

library(tidyverse);library(rio);library(dlnm);library(splines);
library(tsModel);library(mgcv);library(mixmeta);library(RcppRoll)

#IMPORT TEMPERATURE DATA for 270 cities
tem_270 <- read.csv(file = 'tem_data.csv') %>% 
  mutate(date = as.Date(date))
#IMPORT POLLUTANTS DATA
pol <- READ.CSV("Data/Pol/pol.csv")

#CALCULATE PERCENTILES OF THE TEMPERATURE DISTRIBUTION
quan <- tem_270 %>% 
  group_by(city_code) %>% 
  summarise(q2_5 = quantile(l0avg_temp,0.025,na.rm =T),
            q975 = quantile(l0avg_temp,0.975,na.rm =T),
            q5 = quantile(l0avg_temp,0.05,na.rm =T),
            q95 = quantile(l0avg_temp,0.95,na.rm =T),
            q1 = quantile(l0avg_temp,0.01,na.rm =T),
            q99 = quantile(l0avg_temp,0.99,na.rm =T),
            max = max(l0avg_temp,na.rm =T),
            min = min(l0avg_temp,na.rm =T),
            ranges = max-min)
#DEFINE EXTREME COLD TEMPERATRE AND EXTREME HOT TEMPERATURE
value_cold <- mean(quan$q2_5,na.rm =T)
value_hot<-  mean(quan$q975,na.rm =T)

#BOUND OF THE TEMPERATURE DISTRIBUTION
bound <- c(round(mean(quan$min),1),round(mean(quan$max),1))

#LOAD DATA
m <- "Total_major_CVDs"
city_list <- read.csv(paste0("Result-temperature/City_list/",m,"_city2.csv")) %>% 
  distinct(district)

data_use <- import(paste0("Data/",m,"_compont.csv")) 

#COMBINE TEMPERATURE AND HA
data_use <- left_join(data_use,tem_270, by = c("shiqu","city_code","date"))
#COMBINE pollutants AND HA
data_use <- left_join(data_use,pol,by = c("city_code","date"))
regions <- unique(data_use$district)

#ARRANGE THE DATA AS A LIST OF DATA SETS
data <- lapply(regions,function(x) data_use[data_use$district==x,])
names(data) <- regions

#LOOP FOR PM2.5 AND TIS MAJOR CONSTITUENTS AT LAG01
pollag <- c("PM2.5lag01","BClag01","OMlag01","NO3lag01","SO4lag01","NH4lag01","cllag01")

for(pol_i in pollag) {
  
  #DEFINE MATRIX TO STORE THE RESULTS

  coef<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
  coefm1<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
  coefm2<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
  coefm3<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
  coefm4<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
  
  vcovl <- vector("list",length(regions))
  names(vcovl) <- regions
  vcovl1 <- vector("list",length(regions))
  names(vcovl1) <- regions
  vcovl2 <- vector("list",length(regions))
  names(vcovl2) <- regions
  vcovl3 <- vector("list",length(regions))
  names(vcovl3) <- regions
  vcovl4 <- vector("list",length(regions))
  names(vcovl4) <- regions
  
  MMT1 <- vector()
  MMT2 <- vector()
  MMT3 <- vector()
  MMT4 <- vector()
  
  
  #LOOP FOR CITIES
  for(i in regions) {
    # PRINT
    cat(i,"")
    
    
    #LOAD THE DATA
    sub <- data_use %>% 
      filter(district == i) 
    
    #GENERATRE COVARIATES
    sub$date1 <- 1:length(sub$date)
    sub$rhu <- sub$l0rhu
   
    #CALCULATE pm2.5 AND ITS MAJOR CONSTITUENTS AT LAG01
    sub$PM2.5lag01  = roll_mean(sub$PM2.5, n = 2, align = "right", fill = NA)
    sub$BClag01 = roll_mean(sub$BC, n = 2, align = "right", fill = NA)
    sub$NO3lag01 = roll_mean(sub$NO3, n = 2, align = "right", fill = NA)
    sub$SO4lag01 = roll_mean(sub$SO4, n = 2, align = "right", fill = NA)
    sub$NH4lag01 = roll_mean(sub$NH4, n = 2, align = "right", fill = NA)
    sub$OMlag01 = roll_mean(sub$OM, n =2, align = "right", fill = NA)
    sub$cllag01 = roll_mean(sub$cl, n = 2, align = "right", fill = NA)
    
    #COLUM NUMBER OF TEMPERATURE
    col_index <- which(names(sub) == "l0avg_temp")
    #TRIM THE OUTLIERS  
    p99 <- round(quantile(sub[,c(seq(col_index,col_index+4*lag*1,by = 4))], probs = 0.991,na.rm = T),1)
    p1 <- round(quantile(sub[,c(seq(col_index,col_index+4*lag*1,by = 4))], probs = 0.009,na.rm = T),1)
    sub[,c(seq(col_index,col_index+4*lag*1,by = 4))][sub[,c(seq(col_index,col_index+4*lag*1,by = 4))] > p99] <- NA
    sub[,c(seq(col_index,col_index+4*lag*1,by = 4))][sub[,c(seq(col_index,col_index+4*lag*1,by = 4))] <p1] <- NA
    
    #DEFINE QUARTILE GROUPS OF PM2.5 AND TIS MAJOR CONSTITUENTS
    qk<-quantile(sub[,pol_i],na.rm=T) 
    sub$qgroup[sub[,pol_i]<=qk[2]]<-"q1"
    sub$qgroup[sub[,pol_i]>qk[2] & sub[,pol_i]<=qk[3]]<-"q2"
    sub$qgroup[sub[,pol_i]>qk[3] & sub[,pol_i]<=qk[4]]<-"q3"
    sub$qgroup[sub[,pol_i]>qk[4] & sub[,pol_i]<=qk[5]]<-"q4"
    
    
    #DEFINE THE CROSS-BASIS
    lag <- 28
    arglag <- list(knots=logknots(lag,3))
    argvar <- list(fun="ns",df=3)
    cb <- crossbasis(sub$l0avg_temp, lag = lag, argvar=argvar, arglag=arglag)
    sub["pol_adj"] <- sub[,pol_i]
   
    #RUN THE MODEL
    tryCatch({
      
      fit <- gam(PSN_NO~cb:qgroup+ pol_adj+ns(rhu, df=3) + ns(date1,df=8*length(unique(year)))
                 +as.factor(Holiday)
                 +as.factor(dow),
                 control=list(maxit=100),
                 family=quasipoisson,sub,na.action="na.exclude")
      fit_all <- gam(PSN_NO~cb+ns(rhu, df=3) + ns(date1,df=8*length(unique(year))) 
                     +as.factor(Holiday)
                     +as.factor(dow),
                     control=list(maxit=100),
                     family=quasipoisson,sub,na.action="na.exclude")
      
      red <- crossreduce(cb,fit_all,cen=20)
      coef[i,] <- coef(red)
      vcovl[[i]] <- vcov(red)
      
      #EXTRACT COEFFICIENT OF EACH GROUP
      ## group 1
      indexcoef1<-grep("qgroupq1",names(coef(fit)))
      coef1<-coef(fit)[indexcoef1]
      
      indexvcov1<-grep("qgroupq1",rownames(vcov(fit)))
      vcov1<-vcov(fit)[indexvcov1,indexvcov1]
      
      pred1<-crossreduce(cb, coef=coef1, vcov=vcov1, type="overall",model.link="log",cen= cen)
      coefm1[i,]<-coef(pred1)
      vcovl1[[i]]<-vcov(pred1)
      
      mmt1 <- pred1$predvar[which.min(pred1$RRfit)];
      
      
      ## group 2
      indexcoef2<-grep("qgroupq2",names(coef(fit)))
      coef2<-coef(fit)[indexcoef2]
      
      indexvcov2<-grep("qgroupq2",rownames(vcov(fit)))
      vcov2<-vcov(fit)[indexvcov2,indexvcov2]
      
      pred2<-crossreduce(cb, coef=coef2, vcov=vcov2, type="overall",model.link="log",cen=cen)
      coefm2[i,]<-coef(pred2)
      vcovl2[[i]]<-vcov(pred2)
      mmt2 <- pred2$predvar[which.min(pred2$RRfit)];
      
      
      ## group 3
      indexcoef3<-grep("qgroupq3",names(coef(fit)))
      coef3<-coef(fit)[indexcoef3]
      
      indexvcov3<-grep("qgroupq3",rownames(vcov(fit)))
      vcov3<-vcov(fit)[indexvcov3,indexvcov3]
      
      pred3<-crossreduce(cb, coef=coef3, vcov=vcov3, type="overall",model.link="log",cen=cen)
      coefm3[i,]<-coef(pred3)
      vcovl3[[i]]<-vcov(pred3)
      mmt3 <- pred3$predvar[which.min(pred3$RRfit)];
      
      
      ## group 4
      indexcoef4<-grep("qgroupq4",names(coef(fit)))
      coef4<-coef(fit)[indexcoef4]
      
      indexvcov4<-grep("qgroupq4",rownames(vcov(fit)))
      vcov4<-vcov(fit)[indexvcov4,indexvcov4]
      
      pred4<-crossreduce(cb, coef=coef4, vcov=vcov4, type="overall",model.link="log",cen=cen)
      coefm4[i,]<-coef(pred4)
      vcovl4[[i]]<-vcov(pred4)
      mmt4 <- pred4$predvar[which.min(pred4$RRfit)];
      
      MMT1[i] <- mmt1
      MMT2[i] <- mmt2
      MMT3[i] <- mmt3
      MMT4[i] <- mmt4
      
    }, error=function(e){cat("Error",conditionMessage(e), "\n")})    
  }
 
  
  ### META FOR NATIONAL RESULT
  
  meta_nest1<-mixmeta(coefm1~1,vcovl1,method = "reml")
  meta_nest2<-mixmeta(coefm2~1,vcovl2,method = "reml")
  meta_nest3<-mixmeta(coefm3~1,vcovl3,method = "reml")
  meta_nest4<-mixmeta(coefm4~1,vcovl4,method = "reml")
  
  # PREDICT THE POOLED COEFFICIENTS
  avgtmean <- sapply(data,function(x) mean(x$l0avg_temp,na.rm=T))
  rangetmean <- sapply(data,function(x) diff(range(x$l0avg_temp,na.rm=T)))
  
  datanew <- data.frame(avgtmean=mean(tapply(avgtmean,regions,mean)),
                        rangetmean=mean(tapply(rangetmean,regions,mean))) 
  mvpred1 <- predict(meta_nest1,datanew,vcov=T,format="list")
  mvpred2 <- predict(meta_nest2,datanew,vcov=T,format="list")
  mvpred3 <- predict(meta_nest3,datanew,vcov=T,format="list")
  mvpred4 <- predict(meta_nest4,datanew,vcov=T,format="list") 
 
  #RE-CENTERING
  ranges <- t(sapply(data, function(x) 
    range(x$l0avg_temp,na.rm=T)))
  bound <- round(colMeans(ranges,1),1)
  
  tperc <- seq(bound[1],bound[2],length=50)
  argvar1 <-list(fun="ns",df=3)
  cb.tm1 <- crossbasis(tperc,maxlag=lag,argvar=argvar1,arglag=arglag)
  
  # DEFINE INDICATOR FOR CENTERING PERCENTILE FROM AVERAGE ASSOCIATIONF
  bvar <- do.call("onebasis",c(list(x=tperc),attr(cb.tm1,"argvar")))
  xlag <- 0:lag
  
  #PREDICTION FOR EACH POLLUTANTS' GROUP
  cp1 <- crosspred(bvar,coef=mvpred1$fit,vcov=mvpred1$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt1)
  cp2 <- crosspred(bvar,coef=mvpred2$fit,vcov=mvpred2$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt2)
  cp3 <- crosspred(bvar,coef=mvpred3$fit,vcov=mvpred3$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt3)
  cp4 <- crosspred(bvar,coef=mvpred4$fit,vcov=mvpred4$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt4)
  
  #MINIMUM MORBIDITY TEMPERATURE
  mmt1 <- cp1$predvar[which.min(cp1$allRRfit)]
  mmt2 <- cp2$predvar[which.min(cp2$allRRfit)]
  mmt3 <- cp3$predvar[which.min(cp3$allRRfit)]
  mmt4 <- cp4$predvar[which.min(cp4$allRRfit)]
  #REPREDICT AT MINIMUM MORBIDITY TEMPERATURE
  cp1 <- crosspred(bvar,coef=mvpred1$fit,vcov=mvpred1$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt1)
  cp2 <- crosspred(bvar,coef=mvpred2$fit,vcov=mvpred2$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt2)
  cp3 <- crosspred(bvar,coef=mvpred3$fit,vcov=mvpred3$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt3)
  cp4 <- crosspred(bvar,coef=mvpred4$fit,vcov=mvpred4$vcov,from=bound[1],to=bound[2],model.link="log",by=0.1,cen = mmt4)
}

