####################################################################################################
########## nlminb

p1_1v <- rep(0,500)
p1_2v <- rep(0,500)
p1_3v <- rep(0,500)

p2_1v <- rep(0,500)
p2_2v <- rep(0,500)
p2_3v <- rep(0,500)

p3_1v <- rep(0,500)
p3_2v <- rep(0,500)
p3_3v <- rep(0,500)

p32_1v <- rep(0,500)
p32_2v <- rep(0,500)
p32_3v <- rep(0,500)

ps_1v <- rep(0,500)
ps_2v <- rep(0,500)
ps_3v <- rep(0,500)

ps1_1v <- rep(0,500)
ps1_2v <- rep(0,500)
ps1_3v <- rep(0,500)

ps2_1v <- rep(0,500)
ps2_2v <- rep(0,500)
ps2_3v <- rep(0,500)

#pbj_1v <- rep(0,500)
pbj_2v <- rep(0,500)
pbj_3v <- rep(0,500)

rn <- 1:1000
rseed <- base::sample(rn,500,replace = FALSE)
set.seed(542)
for(i in 1:500){
  set.seed(rseed[i])
  ######
  inv.func <- function(u){
    return(tan(pi*(u-1/2)))
  }
  
  inversion <- function(n,inv.func){
    u <- runif(n)
    inv.func(u)
  }
  
  e_01 <- inversion(300,inv.func)
  ######
  sage <- round(runif(300, 12, 64))
  sage2 <- sage*sage
  err <- e_01
  logtime <- 1.3 + sage*0.15 - sage2 * 0.0016 + err
  stime <- round(10^logtime)
  sstatus <- base::sample(c(rep(1,200),rep(0,100)), 300,
                          replace = FALSE)
  data_i <- data.frame(stime,sage,sage2,sstatus)
  ######data_test <- arrange(data_test,time,desc(status))
  group1 <- filter(data_i, sage < 30) %>%
    arrange(stime) %>% mutate(Group = "1")
  group2 <- filter(data_i, sage >=30 & sage <= 39) %>%
    arrange(stime) %>% mutate(Group = "2")
  group3 <- filter(data_i, sage >=40 & sage <= 49) %>%
    arrange(stime) %>% mutate(Group = "3")
  group4 <- filter(data_i, sage >= 50) %>%
    arrange(stime) %>% mutate(Group = "4")
  ######
  group1[which.min(group1$stime),"sstatus"] <- 0
  group2[which.min(group2$stime),"sstatus"] <- 0
  group3[which.min(group3$stime),"sstatus"] <- 0
  group4[which.min(group4$stime),"sstatus"] <- 0
  ######
  group1 <- gweight(group1,group1$stime,group1$sstatus)
  group2 <- gweight(group2,group2$stime,group2$sstatus)
  group3 <- gweight(group3,group3$stime,group3$sstatus)
  group4 <- gweight(group4,group4$stime,group4$sstatus)
  data_i <- rbind(group1,group2,group3,group4)
  data_i$wp <- data_i$w1/(4*log10(data_i$stime))
  ######
  
  data_i$status <- data_i$sstatus
  data_i$time <- data_i$stime
  data_i_t <- get_synthetic(data = data_i, 
                            time = data_i$time,
                            status = data_i$status, w1 = data_i$w1)
  data_i <- data_i_t
  ######
  data_i_uc <- filter(data_i, sstatus == 1)
  try(lm01 <- lm(log(stime) ~ sage + sage2,data = data_i_uc), silent = TRUE)
  beta_i <- as.vector(lm01$coefficients)
  ######
  
  ###### Synthetic
  try(fit_syn <- lm(log10(Ag1) ~ sage + sage2, data = data_i), silent = TRUE)
  ps_1v[i]<- try(as.vector(fit_syn$coefficients)[1], silent = TRUE)
  ps_2v[i]<- try(as.vector(fit_syn$coefficients)[2], silent = TRUE)
  ps_3v[i]<- try(as.vector(fit_syn$coefficients)[3], silent = TRUE)
  
  p_ms <- c(ps_1v[i],ps_2v[i],ps_3v[i])
  print(p_ms)
  
  ###### New Synthetic 01
  
  try(opt_i <- optimx(beta_i,f_function_SYN,method='nlminb'),silent=TRUE)
  ps1_1v[i]<- opt_i$p1
  ps1_2v[i]<- opt_i$p2
  ps1_3v[i]<- opt_i$p3
  
  p_ms1 <- c(ps1_1v[i],ps1_2v[i],ps1_3v[i])
  print(p_ms1)
  
  ###### NewMethod using W_Prime
  try(opt_i <- optimx(beta_i,f_function_MD01_i,method='nlminb'), silent = TRUE)
  p1_1v[i] <- opt_i$p1
  p1_2v[i] <- opt_i$p2
  p1_3v[i] <- opt_i$p3
  
  p_m1 <- c(p1_1v[i],p1_2v[i],p1_3v[i])
  print(p_m1)
  
  ###### New Synthetic 02
  
  try(opt_i <- optimx(beta_i,f_function_SYN_1,method='nlminb'),silent=TRUE)
  ps2_1v[i]<- opt_i$p1
  ps2_2v[i]<- opt_i$p2
  ps2_3v[i]<- opt_i$p3
  
  p_ms2 <- c(ps2_1v[i],ps2_2v[i],ps2_3v[i])
  print(p_ms2)
  
  ###### Type I Discrete using W1
  try(opt_i <- optimx(beta_i,f_function_MD02_i,method='nlminb'), silent = TRUE)
  p2_1v[i] <- opt_i$p1
  p2_2v[i] <- opt_i$p2
  p2_3v[i] <- opt_i$p3
  
  p_m2 <- c(p2_1v[i],p2_2v[i],p2_3v[i])
  print(p_m2)
  
  ###### Least Square Method by Zhou
  try(opt_i <- optimx(beta_i,f_function_MD03_i,method='nlminb'), silent = TRUE)
  p3_1v[i] <- opt_i$p1
  p3_2v[i] <- opt_i$p2
  p3_3v[i] <- opt_i$p3
  
  p_m3 <- c(p3_1v[i],p3_2v[i],p3_3v[i])
  print(p_m3)
  
  ###### Least Absolute Deviation by Zhou 
  try(opt_i <- optimx(beta_i,f_function_MD032_i,method='nlminb'), silent = TRUE)
  p32_1v[i] <- opt_i$p1
  p32_2v[i] <- opt_i$p2
  p32_3v[i] <- opt_i$p3
  
  p_m32 <- c(p32_1v[i],p32_2v[i],p32_3v[i])
  print(p_m32)
  
  ###### B-J
  lss_fit1 <- try(lss(Surv(log10(stime), sstatus) ~ sage + sage2, data = data_i,
                      mcsize=1, trace = F, cov = T), silent = TRUE)
  
  pbj_2v[i] <- try(lss_fit1$lse[1], silent = TRUE)
  pbj_3v[i] <- try(lss_fit1$lse[2], silent = TRUE)
  
  p_mbj <- try(c(pbj_2v[i],pbj_3v[i]), silent = TRUE)
  try(cat("Number", i, ":", print(p_mbj)), silent = TRUE)
  print("      ")
}


mean(ps_1v, na.rm = TRUE)
mean(ps_2v, na.rm = TRUE)
mean(ps_3v, na.rm = TRUE)

mean(ps1_1v, na.rm = TRUE)
mean(ps1_2v, na.rm = TRUE)
mean(ps1_3v, na.rm = TRUE)

mean(ps2_1v, na.rm = TRUE)
mean(ps2_2v, na.rm = TRUE)
mean(ps2_3v, na.rm = TRUE)

mean(p1_1v)
mean(p1_2v)
mean(p1_3v)


mean(p2_1v)
mean(p2_2v)
mean(p2_3v)


mean(p3_1v)
mean(p3_2v)
mean(p3_3v)


mean(p32_1v)
mean(p32_2v)
mean(p32_3v)



pbj_2v <- as.numeric(pbj_2v)
pbj_3v <- as.numeric(pbj_3v)
mean(pbj_2v, na.rm = TRUE)
mean(pbj_3v, na.rm = TRUE)



datap_nlminb <- data.frame(p1_1v,p1_2v,p1_3v,
                           p2_1v,p2_2v,p2_3v,
                           p3_1v,p3_2v,p3_3v,
                           p32_1v,p32_2v,p32_3v,
                           pbj_2v,pbj_3v,
                           ps_1v,ps_2v,ps_3v,
                           ps1_1v,ps1_2v,ps1_3v,
                           ps2_1v,ps2_2v,ps2_3v)

write.csv(datap_nlminb, "dataD-nlminb-1.csv") 