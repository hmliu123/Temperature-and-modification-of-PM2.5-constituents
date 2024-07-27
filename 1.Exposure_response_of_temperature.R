###########################################################################################################################################################
#Codes for "Non-optimum ambient temperature and cause-specific hospital admissions for cardiovascular diseases in the context of particulate air pollution"
#Authors for codes: Huimeng Liu, Jian Lei, Yunxing Jiang, Lijun Bai et.al.
#Correspondence to Shaowei Wu, Xi Chen.
##########################################################################################################################################################

######################################################################################################################################
#Overall exposure-response relationships between non-optimum temperature and hospital admissions (HAs) for major CVDs
##########################################################################################################################################
library(tidyverse)
library(lubridate)
library(dlnm) 
library(mvmeta) 
library(splines)
library(mgcv)

#IMPORT TEMPERATURE DATA for 270 cities
tem_270 <- read.csv(file = 'tem_data.csv') %>% 
  mutate(date = as.Date(date))

#MAJOR CVDS
subtype <- c("Total_major_CVDs","CHD","Cardiac arrhythmias","Heart failure","Stroke","ACS","AMI","Angina","Hemorrhagic stroke","Ischemic stroke")
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

for (m in subtype) {
  #PRINT
  cat(m," ")
  #IMPORT HA DATA
  data_use <- import(paste0("/Data/","_compont.csv")) 
  
  #COMBINE TEMPERATURE AND HA
  data_use <- left_join(data_use,tem_270, by = c("shiqu","city_code","date"))
  regions <- unique(data_use$district)
 
  
  #ARRANGE THE DATA AS A LIST OF DATA SETS
  data <- lapply(regions,function(x) data_use[data_use$district==x,])
  names(data) <- regions
  
  
  #SET EMPTY MATRIX TO STORE THE RESULTS
  yall<- matrix(NA,length(regions),3,dimnames=list(regions,paste("b",seq(3),sep=""))) 
  yhot <- matrix(NA,length(regions),5,dimnames=list(regions,paste("b",seq(5),sep="")))
  ycold <- matrix(NA,length(regions),5,dimnames=list(regions,paste("b",seq(5),sep="")))
  
  Sall <- vector("list",length(regions))
  names(Sall) <- regions
  
  Shot <- c()
  Scold <- c()
  
  exlow.rr <- vector()
  exlow.rrlow <- vector()
  exlow.rrhigh <- vector()
  exlow.coef <- vector()
  exlow.se <- vector()
  exhigh.rr <- vector()
  exhigh.rrlow <- vector()
  exhigh.rrhigh <- vector()
  exhigh.coef <- vector()
  exhigh.se <- vector()
  MMT <- vector()
  hightemp <- vector()
  lowtemp <- vector()
  
  #LOOP FOR CITIES
  for(i in regions) {
    # PRINT
    cat(i,"")
    
    # LOAD CITY-SPECIFIC DATA
    sub <- data[[i]]
    
    #GENERATRE COVARIATES
    sub$date1 <- 1:length(sub$date)
    sub$rhu <- sub$l0rhu
    #######################################################################################
    #TRIM THE OUTLIERS 
    p99 <- round(quantile(sub[,c(seq(29,29+4*lag*1,by = 4))], probs = 0.991,na.rm = T),1)
    p1 <- round(quantile(sub[,c(seq(29,29+4*lag*1,by = 4))], probs = 0.009,na.rm = T),1)
    
    sub[,c(seq(29,29+4*lag*1,by = 4))][sub[,c(seq(29,29+4*lag*1,by = 4))] > p99] <- NA
    sub[,c(seq(29,29+4*lag*1,by = 4))][sub[,c(seq(29,29+4*lag*1,by = 4))] <p1] <- NA
    
    #DEFINE THE CROSS-BASIS
    lag <- 28
    temp <- sub[,c(seq(29,29+4*lag*1,by = 4))]
    arglag <- list(knots=logknots(lag,3))
    argvar <- list(fun="ns",df=3) 
    cb <- crossbasis(sub$l0avg_temp, lag = lag, argvar=argvar, arglag=arglag)
    
    #RUN THE MODELS 
    tryCatch({
      mfirst <- gam(PSN_NO~cb+ ns(rhu, df=3) +ns(date1,df=8*length(unique(year)))
                    +as.factor(Holiday)
                    +as.factor(dow),
                    family=quasipoisson,sub,na.action="na.exclude")
      
      # PREDICTION AND REDUCTION TO OVERALL CUMULATIVE EXPOSURE-RESPONSE
      pred <- crosspred(cb, mfirst,by = 0.1,cumul = TRUE);
      mmt <- pred$predvar[which.min(pred$allRRfit)];
      
      crall <- crossreduce(cb,mfirst,cen=20,by = 0.1) #
      crhot <- crossreduce(cb,mfirst,type="var",value=value_hot,cen=20)
      crcold <- crossreduce(cb,mfirst,type="var",value=value_cold,cen=20)
      
      # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
      yall[i,] <- coef(crall)
      Sall[[i]] <- vcov(crall)
      yhot[i,] <- coef(crhot)
      Shot[[i]] <- vcov(crhot)
      ycold[i,] <- coef(crcold)
      Scold[[i]] <- vcov(crcold)
      
      
      exlow.rr[i] <- crall$RRfit[paste(p1)]
      exlow.rrlow[i]<-crall$RRlow[paste(p1)]
      exlow.rrhigh[i]<-crall$RRhigh[paste(p1)]
      
      exhigh.rr[i]<-crall$RRfit[paste(p99)]
      exhigh.rrlow[i]<-crall$RRlow[paste(p99)]
      exhigh.rrhigh[i]<-crall$RRhigh[paste(p99)]
      
      exhigh.se[i] <- crall$se[paste(p99)]
      exlow.se[i] <- crall$se[paste(p1)]
      
      hightemp[i] <- quantile(sub[,c(seq(29,29+4*lag*1,by = 4))], probs = 0.99,na.rm = T)
      lowtemp[i] <- quantile(sub[,c(seq(29,29+4*lag*1,by = 4))], probs = 0.01,na.rm = T)
      
      MMT[i] <- mmt
      
      
    }, error=function(e){cat("Error",conditionMessage(e), "\n")})    
  }
  
  mmt <- round(median(MMT),1)
  HIGH <- round(median(hightemp),1)
  Low <- round(median(lowtemp),1)
  
  
  
  #META
  method <- "reml"
  #OVERALL
  mvall <- mvmeta(yall~1,Sall,method=method) 
  summary(mvall)
  #HOT
  mvhot <- mvmeta(yhot~1,Shot,method=method)
  summary(mvhot)
  #COLD
  mvcold <- mvmeta(ycold~1,Scold,method=method)
  summary(mvcold)
  
  # RE-CENTERING
  # DEFINE RELATED AVERAGE TEMPERATURES
 
  tperc <- seq(bound[1],bound[2],by = 0.1)
  argvar1 <- list(fun="ns",df=3) 
  cb.tm1 <- crossbasis(tperc,lag=lag,argvar=argvar1,arglag=arglag)
  
  # DEFINE INDICATOR FOR CENTERING PERCENTILE FROM AVERAGE ASSOCIATION
  bvar <- do.call("onebasis",c(list(x=tperc),attr(cb.tm1,"argvar")))
  xlag <- 0:lag
  blag <- do.call("onebasis",c(list(x=xlag),attr(cb.tm1,"arglag")))
  
  #PREDICT THE CITY-SPECIFIC CUMULATIVE ASSOCIATIONS
  regall <- lapply(seq(nrow(yall)),function(i) 
    crosspred(bvar,coef=yall[i,],
              vcov=Sall[[i]],model.link="log",cen=23))
  reghot <- lapply(seq(nrow(yhot)),function(i) 
    crosspred(blag,coef=yhot[i,],
              vcov=Shot[[i]],model.link="log",cen=23))
  regcold <- lapply(seq(nrow(ycold)),function(i) 
    crosspred(blag,coef=ycold[i,],
              vcov=Scold[[i]],model.link="log",cen=23))
  
  
  #PREDICT THE OVERALL CUMULATIVE ASSOCIATIONS
  cpall.total <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),
                           model.link="log",by=0.1,from=bound[1],to=bound[2],cen=mmt)
 
  #MINIUM MORBIDITY TEMPERATURE
  mmt <- cpall.total$predvar[which.min(cpall.total$allRRfit)]
  #REPREDICT AT MINIMUM MORBIDITY TEMPERATURE
  cpall.total <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),
                           model.link="log",by=0.1,from=bound[1],to=bound[2],cen=mmt)
  #PREDICT THE LAG PATTERNS OF EXTREME HOT TEMPERATURES
  cphot.total <- crosspred(blag,coef=coef(mvhot),vcov=vcov(mvhot),
                           model.link="log",at=0:lag)
  #PREDICT THE LAG PATTERNS OF EXTREME COLD TEMPERATURES
  cpcold.total <- crosspred(blag,coef=coef(mvcold),vcov=vcov(mvcold),
                            model.link="log",at=0:lag,cen=Low)
  
 }


#PLOT OVERALL EXPOSURE-RESPONSE CURVES
xlab <- paste("Daily Mean Temperature (℃)")
plot(cpall.total,type="n",ylab="Relative Risk",
     ci.arg = list(col = c("#D0E7Ed")),cex.axis = 1.3, cex.lab = 1.3)

lines(cpall.total,col=1,lwd=1)
lines(c(mmt,mmt),c(0,10),lty = 2,lwd = 2,col = "red")
mtext(m,cex=1.2,side = 3,line = 0.5)
mtext(paste0(mmt,"℃"),at = c(mmt+4),line = -0.8,cex=1)

#PLOT LAG PATTERN FOR EXTREME HOT TEMPERATURE
plot(cphot.total,type="n",ylab="Relative Risk",cex.axis = 1.3, cex.lab = 1.5,
     ci.arg = list(col = c("#D0E7Ed")),ylim = c(0.95,1.04))
lines(cphot.total,col=1,lwd=1)
mtext(m,cex=1,line = 0)
#PLOT LAG PATTERN FOR EXTREME HOT TEMPERATURE

plot(cpcold.total,type="n",ylab="Relative Risk",cex.axis = 1.3, cex.lab = 1.5,
     ci.arg = list(col = c("#D0E7Ed")),ylim = c(0.85,1.10))
lines(cpcold.total,col=1,lwd=1)
mtext(m,cex=1,line = 0)


p1 <-round(mean(quan$q1),1)
p2_5 <- round(mean(quan$q25),1)
p5 <- round(mean(quan$q5),1)
p95 <- round(mean(quan$q95),1)
p975 <- round(mean(quan$q975),1)
p99 <-round(mean(quan$q99),1)


#OUT PUT ESTIMATION
results <- round(cbind(cpall.total$allRRfit,cpall.total$allRRlow,cpall.total$allRRhigh)
                 [c(as.character(p1),as.character(p2_5),as.character(p5),
                    as.character(p95),as.character(p975),as.character(p99),
                    as.character(mmt)),],digits=3)

