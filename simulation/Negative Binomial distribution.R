library(mgcv); library(MASS)
library(VGAM)
library(MASS)
library(rmutil)
library(actuar)
library(boot)
library(sandwich)
library(lmtest)
library(robust)
library(robustbase)
library(openxlsx)

n=200
iters=500

glm_coef <- function(data, indices) {
  fit <- glm(y~x, family = poisson(), data = data[indices,])
  return(coef(fit))
}

Simu <- function(filepath, filenumber, beta_1) {
  
  iters <- filenumber
  
  error_iters <- vector() # To store iterations that encounter errors
  
  estall =matrix(999, iters, 6);
  colnames(estall)<-c('intercept','se_intercept','p_intercept','beta1','se_beta1','p_beta1')
  
  estall2=matrix(999, iters, 6);
  colnames(estall2)<-c('intercept','se_intercept','p_intercept','beta1','se_beta1','p_beta1')
  
  estall3=matrix(999, iters, 6);
  colnames(estall3)<-c('intercept','beta1','se_intercept','se_beta1','p_intercept','p_beta1')

  estall4=matrix(999, iters, 6);
  colnames(estall4) = c('intercept','beta1','se_intercept','se_beta1','p_intercept','p_beta1')
  
  # Initialize other matrices similarly...
  pb <- txtProgressBar(min = 0, max = iters, style = 3)
  suppressWarnings(for (III in 1:iters) {
    tryCatch({
      path <- paste0(filepath, III, ".csv")
      data <- read.csv(path)
      datam <- as.matrix(data)
      
      #nb
      datas=data.frame(datam); 
      fitg<-glm.nb(y~x,link="log",data=datas)
      estall[III,1:6]=c(summary(fitg)$coefficients[1,-3], summary(fitg)$coefficients[2,-3])
      
      #poi
      fitg2<-glm(y~x,poisson(link="log"),data=datas)
      estall2[III,1:6]=c(summary(fitg2)$coefficients[1,-3], summary(fitg2)$coefficients[2,-3])
      
      #poission bootstrap
      fitg3coef <- glm_coef(data= datas)
      fitg3results <- boot(data = datas, statistic = glm_coef, R = 1000)
      # Calculate bootstrap standard errors
      fitg3boot_se <- apply(fitg3results$t, 2, sd)
      Zvalue = fitg3coef/fitg3boot_se
      fitg3Pvalue = 2 * (1 - pnorm(abs(fitg3coef/fitg3boot_se)))
      estall3[III,1:6] = c(fitg3coef,fitg3boot_se,fitg3Pvalue)
   
      #poisson sandwich
      # Fit a Poisson GLM
      #model <- glm(y~x, family = poisson(), data = datas)
      estall4[III,1:6]  = as.numeric(coeftest(fitg2, vcov. = vcovHC(fitg2, type = "HC3"))[,c(1,2,4)])
      
      
    }, error = function(e) {
      # If an error occurs, store the iteration number
      error_iters <- c(error_iters, III)
      cat("Error in iteration:", III, "\n")
    })
    setTxtProgressBar(pb, III)
  })
  
  # Optionally print or handle iterations that encountered errors
  if (length(error_iters) > 0) {
    cat("Errors occurred in iterations:", paste(error_iters, collapse = ", "), "\n")
  }
  close(pb)
  # Return all your matrices and error_iters at the end
  return(list(estall,estall2,estall3,estall4, error_iters = error_iters))
}

type1conclusion <- function(list_estall,beta_0,alpha_level = 0.05){
  
  type1result<-function(estall,beta0 = 1, beta1 = 0,alphalevel = 0.05){
    type1=as.numeric(estall[,'p_beta1']<alphalevel)
    zvalue = qnorm(1 - alphalevel / 2)
    cp_beta0=as.numeric(abs(estall[,'intercept']-beta0)<zvalue*estall[,'se_intercept'])
    cp_beta1=as.numeric(abs(estall[,'beta1']-beta1)<zvalue*estall[,'se_beta1'])
    estall=cbind(estall,type1,cp_beta0,cp_beta1)
    result=c(mean(estall[,'intercept']-beta0),sd(estall[,'intercept']),mean(estall[,'se_intercept']), mean((estall[,'intercept']-beta0)^2),
             mean(estall[,'beta1']-beta1),sd(estall[,'beta1']),mean(estall[,'se_beta1']),mean((estall[,'beta1']-beta1)^2),
             mean(type1),mean(cp_beta0),mean(cp_beta1))
    names(result)=c('bias_beta0','sd_beta0','se_beta0','mse_beta0','bias_beta1','sd_beta1','se_beta1','mse_beta1','type1','cp_beta0','cp_beta1')
    return(result)
  }
  clean_error <- function(df){
    df = data.frame(df)
    # Identify columns containing the value 999
    rows_with_999 <- apply(df,1, function(x) any(x == 999))
    # Drop these columns
    df_cleaned <- df[ !rows_with_999 , ]
    return(df_cleaned)
  }
  
  estall<-clean_error(list_estall[[1]])
  estall2<-clean_error(list_estall[[2]])
  estall3<-clean_error(list_estall[[3]])
  estall4<-clean_error(list_estall[[4]])
  
  result_nbr_0<- type1result(estall,beta0 = beta_0 ,alphalevel = alpha_level)
  result_pr_0<- type1result(estall2,beta0 = beta_0,alphalevel = alpha_level)
  result_prboot_0<- type1result(estall3,beta0 = beta_0,alphalevel = alpha_level)
  result_prsandwich_0 <- type1result(estall4,beta0 = beta_0,alphalevel = alpha_level)
  
  result0=data.frame(result_nbr_0,result_pr_0,
                     result_prboot_0,result_prsandwich_0)
  return(result0)
}

Powerconclusion<-function(list_estall,beta_0=0,beta_1,alpha_level = 0.05){
  
  powerresult<-function(est,beta0, beta1,alphalevel = 0.05){
    power=as.numeric(est[,'p_beta1']<alphalevel)
    zvalue = qnorm(1 - alphalevel / 2)
    cp_beta0=as.numeric(abs(est[,'intercept']-beta0)<zvalue*est[,'se_intercept'])
    cp_beta1=as.numeric(abs(est[,'beta1']-beta1)<zvalue*est[,'se_beta1'])
    est=est[,1:6]
    est=cbind(est,power,cp_beta0,cp_beta1)
    result=c(mean(est[,'intercept']-beta0),sd(est[,'intercept']),
             mean(est[,'se_intercept']),mean((est[,'intercept']-beta0)^2),
             mean(est[,'beta1']-beta1),sd(est[,'beta1']),mean(est[,'se_beta1']),
             mean((est[,'beta1']-beta1)^2),
             mean(power),mean(cp_beta0),mean(cp_beta1))
    names(result)=c('bias_beta0','sd_beta0','se_beta0','mse_beta0','bias_beta1','sd_beta1','se_beta1','mse_beta1','power','cp_beta0','cp_beta1')
    return(result)}
  
  clean_error <- function(df){
    df = data.frame(df)
    # Identify columns containing the value 999
    rows_with_999 <- apply(df,1, function(x) any(x == 999))
    # Drop these columns
    df_cleaned <- df[ !rows_with_999 , ]
    return(df_cleaned)
  }
  
  estall<-clean_error(list_estall[[1]])
  estall2<-clean_error(list_estall[[2]])
  estall3<-clean_error(list_estall[[3]])
  estall4<-clean_error(list_estall[[4]])
  
  result_nbr<- powerresult(estall,beta0 = beta_0,beta1 = beta_1,alphalevel = alpha_level )
  result_pr<- powerresult(estall2,beta0 = beta_0,beta1 = beta_1,alphalevel = alpha_level )
  result_prboot<- powerresult(estall3,beta0 = beta_0,beta1 = beta_1,alphalevel = alpha_level )
  result_prsandwich <- powerresult(estall4,beta0 = beta_0,beta1 = beta_1,alphalevel = alpha_level )
  
  result=data.frame(result_nbr,result_pr,
                    result_prboot,result_prsandwich)
  
  return(result)
  
}
#####################################################################################################
#data generation type I error
#NB
filedir=paste0("filepath/NB",beta1,"/")
if (!file.exists(filedir)) {
  dir.create(filedir)
}
datam<-matrix(0,nrow=n,ncol=2)
beta0= 2
beta1= 0
beta=c(beta0, beta1)

# Generate the data
filedir=paste0("filepath/NB",beta1,"/")
if (!file.exists(filedir)) {
  dir.create(filedir)
}
datam<-matrix(0,nrow=n,ncol=2)
for (III in 1:iters){
  set.seed(III+700)
  
  x=c(rep(0,0.5*n),rep(1,0.5*n))
  eta<-beta[1]+ beta[2]* x
  
  datam[,1]<- rnbinom(n,mu=exp(eta),size=4)
  datam[,2]<- x
  
  filename1=paste(filedir, III, ".csv", sep="")	# Generate comma delimited file
  
  data.cost.2=as.matrix(datam)
  colnames(data.cost.2)<-c('y','x')
  write.matrix(data.cost.2, filename1, sep=",")
}
groups=c("NB")

beta1=beta[2]
result_distpower=list() 
for (i in 1:length(groups)) {
  distype=groups[i]
  print(distype)
  
  distfilepathpower=paste0(filepath,distype,beta1,"/")
  Simu_dist_resultpower<-Simu(filepath = distfilepathpower,filenumber =500,beta_1 = beta[2])
  result_distpower[[i]] <- Powerconclusion(Simu_dist_resultpower,beta_0=beta[1],beta_1 =beta[2],alpha_level = 0.05)
  names(result_distpower)[i]=paste0(distype,beta1)
}

result_distpower
lapply(result_distpower)
####################################################################################################
sheetnames = paste0(groups,'0')
# Create a new workbook
wb <- createWorkbook()

# Assuming your list of data frames is called "df_list"
for (i in seq_along(result_disttype1)) {
  # Create a new sheet in the workbook
  addWorksheet(wb, sheetName = sheetnames[i])
  
  # Write the data frame to the sheet, including row names
  writeData(wb, sheet = i, x = result_disttype1[[i]], startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
  
}

# Save the workbook to a file
saveWorkbook(wb, file = "NBtypeIerror.xlsx",overwrite=T)
####################################################################################################
#data generation power
#NB
filedir=paste0("filepath/NB",beta1,"/")
if (!file.exists(filedir)) {
  dir.create(filedir)
}
datam<-matrix(0,nrow=n,ncol=2)
beta0= 2
beta1= 0.1
beta=c(beta0, beta1)

for (III in 1:iters){
  set.seed(III)
  
  x=c(rep(0,0.5*n),rep(1,0.5*n))
  eta<-beta[1]+ beta[2]* x
  
  datam[,1]<- rnbinom(n,mu=exp(eta),size=4)
  datam[,2]<- x
  
  filename1=paste(filedir, III, ".csv", sep="")	# Generate comma delimited file
  
  data.cost.2=as.matrix(datam)
  colnames(data.cost.2)<-c('y','x')
  write.matrix(data.cost.2, filename1, sep=",")
}
groups=c("NB")

beta1=beta[2]
result_distpower=list() 
for (i in 1:length(groups)) {
  distype=groups[i]
  print(distype)
  
  distfilepathpower=paste0(filepath,distype,beta1,"/")
  Simu_dist_resultpower<-Simu(filepath = distfilepathpower,filenumber =500,beta_1 = beta[2])
  result_distpower[[i]] <- Powerconclusion(Simu_dist_resultpower,beta_0=beta[1],beta_1 =beta[2],alpha_level = 0.05)
  names(result_distpower)[i]=paste0(distype,beta1)
}

result_distpower
lapply(result_distpower)
####################################################################################################
sheetnames = paste0(groups,'0')
# Create a new workbook
wb <- createWorkbook()

# Assuming your list of data frames is called "df_list"
for (i in seq_along(result_disttype1)) {
  # Create a new sheet in the workbook
  addWorksheet(wb, sheetName = sheetnames[i])
  
  # Write the data frame to the sheet, including row names
  writeData(wb, sheet = i, x = result_disttype1[[i]], startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE)
  
}

# Save the workbook to a file
saveWorkbook(wb, file = "NBpower.xlsx",overwrite=T)



