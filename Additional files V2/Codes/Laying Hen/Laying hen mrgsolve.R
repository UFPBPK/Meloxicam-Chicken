## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(tidyr)       # R-package for tidy messy data
library(tidyverse)   # R-package for tidy messy data
library(truncnorm)   # R package for Truncated normal distribution
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(ggpubr)      # R package for plotting the data

## Set working direction to the data files


library(mrgsolve)  


SolvePBPKOTC <- '
$PARAM @annotated

 
// Oral absorption rate constants
  
Kst     : 2        :                  /h, gastric emptying rate constant
Kint    : 0.2      :                  /h, intestinal transit rate constant.

// IV infusion rate constants
Timeiv  : 0.01     :                  IV injection/infusion time (h)

// Urinary elimination rate constant adjusted by body weight
KurineC : 0.1      :                  L/h/kg, Table 2, Li M et al. 2019

// Metabolic rate constant adjusted by body weight
KmC     : 0.0025   :                 /h/kg, metabolic rate constant

// Dosing, multiple oral gavage
Ka      : 0.5      :                  intestinal absorption rate constant (/h) 


  
  
// Physiological parameters
// Blood flow rates
QCC     : 9.91      :                 Cardiac output (L/h/kg), Table 13, Wang et al. 2021
QLC     : 0.2526    :                 Fraction of flow to the liver, Table 16, Wang et al. 2021
QKC     : 0.2012    :                 Fraction of flow to the kidneys, Table 16, Wang et al. 2021
QFC     : 0.015     :                 Fraction of flow to the fat, Zeng et al. 2019
QMC     : 0.35      :                 Fraction of flow to the muscle, Zeng et al. 2019
QovaryC : 0.1646    :                 Fraction of flow to the ovary, Table 20, Wang et al. 2021


// Tissue  volumes
BW      : 2.2       :                 Body weight(kg)
VLC     : 0.0249    :                 Fractional liver tissue, Table 3, Wang et al. 2021
VKC     : 0.0076    :                 Fractional kidney tissue, Table 3, Wang et al. 2021
VFC     : 0.05      :                 Fractional fat tissue, Zeng et al. 2019
VMC     : 0.4015    :                 Fractional muscle tissue, Table 2, Wang et al. 2021
VovaryC : 0.0191    :                 Fractional ovary tissue, Table 3, Wang et al. 2021
VbloodC : 0.0483    :                 Hematocrit volume, fraction of BW, Table 22, Wang et al. 2021
Htc     : 0.314     :                 Hematocrit volume, fraction of BW, Table 22, Wang et al. 2021

  
// Mass Transfer Parameters (Chemical-specific parameters)
// Partition coefficients(PC, tissue:plasma)
PL      : 0.6116    :                 Liver plasma PC, calculated using AUCtissue:AUCplasma method from Dutch et al. 2021
PK      : 0.6183    :                 Kidney plasma PC, calculated using AUCtissue:AUCplasma method from  Dutch et al. 2021
PM      : 0.0456    :                 Muscle plasma PC, breast muscle, calculated using AUCtissue:AUCplasma method from  Dutch et al. 2021 
Povary  : 0.6116    :                 Richly persused tissue
PF      : 0.0174    :                 Fat plasma PC, calculated using AUCtissue:AUCplasma method from  Dutch et al. 2021
PRest   : 0.479     :                 Rest of the body tissues:plasma PC, Table 2, Li M et al. 2019
PB      : 0.994     :                 Percentage of drug bound to plasma proteins (DrugBank, 2020)


// Egg yolk and white parameters

Ay   : 25                             : apparent maximum follicle weight g
Ky   : 0.08                           : transport constant into yolk
Kw   : 0.025                          : transport constant into albumen
tlag : 1*24                           : h
tsig : 2*24                           : h
s    : 1/24                           : /h
talbumen : 10                         : h
tlay1  : 12                           : tlay1 default
Rwhitefor : (34*0.001)/10             : rate of albumen formation Kg/h 



$MAIN


// Cardiac output and blood flows to tissues(L/h)
double QC = QCC * BW;                  // Cardiac output
double QL = QLC * QC;                  // Liver
double QK = QKC * QC;                  // Kidney
double QF = QFC * QC;                  // Fat
double QM = QMC * QC;                  // Muscle
double Qovary = QovaryC * QC;          // Ovary
double QRest = QC - QK - QL - QF - QM - Qovary; // Rest of the tissues


//  Tissue volumes (L)
double VL = VLC * BW;                  // Liver
double VK = VKC * BW;                  // Kidney
double VF = VFC * BW;                  // Fat
double VM = VMC * BW;                  // Muscle
double Vblood = VbloodC*BW;            // Blood
double Vplasma = VplasmaC * BW;        // Plasma
double Vovary = VovaryC * BW;          // Ovary
double VRest = 1 * BW - VL - VK - VF - VM - Vovary - Vblood; // Rest of the tissues
double VplasmaC = (1-Htc) * VbloodC;   // Plasma volume, fraction of BW, Table 22, Wang et al. 2021
double Kmet = KmC * BW;	               // Metabolic rate constant adjusted by body weight
double Kurine = KurineC * BW;          // # L/h Urinary elimination rate
double Free = 1-PB;                    //Free drug percentage


// Egg yolk and white functions
double tlay   =  tlay1*24;                          // h
double Wy     = (25/(1 + exp(-s*(TIME - (tlay - tsig - tlag)))))*0.001;  // weight of the follicle Kg 
double Cyolk  = Ayolk/Wy;
double Ww = Rwhitefor * (TIME >= tlay-tlag)*(TIME <= tlay - tlag + talbumen);
double Cwhite    = Awhite /AWw;   

      
$CMT ADOSE AST AI Acolon AAO Aplas_free AL AK Aurine Amet AM Aovary Ayolk Awhite AWw AUCCV AUCCL AUCCK AUCCM AUCCF AUCCovary AF ARest Aegg_ex


             
$ODE

dxdt_ADOSE=0;

   // The changing rate of the amount of dose via oral

 
 // Concentration of the chemical in vein and tissue compartments
      double CVL = AL/(VL * PL);      // Concentration in liver / PC of plasma liver
      double CVK = AK/(VK * PK);      // Concentration in Kidney / PC of plasma kidney
      double CVF = AF/(VF * PF);   
      double CVRest = ARest/(VRest * PRest); // Concentration in Rest  / PC of plasma rest tissue
      double CVM = AM/(VM * PM);      // perfusion limited model.
      double CVovary = Aovary/(Vovary * Povary);  // Concentration in ovary/PC of plasma: rest tissue
      double CV = ((QL * CVL  + QK  * CVK  + QF * CVF  + QM * CVM  + Qovary * CVovary + QRest * CVRest)/QC);
      double Cplas_free = Aplas_free/Vplasma;  // Free drug concentration in plasma = amount of chem in plasma / volume of plasma
      double Cplasma = Cplas_free/Free ;
      double CL = AL/VL;              // Concentration in liver
      double CK = AK/VK;              // Concentration in kidney
      double CM = AM/VM;              // Concentration in muscle
      double CF = AF/VF;              // Concentration in fat 
      double Covary = Aovary/Vovary;  // Concentration in ovary 
      double CRest = ARest/VRest;  

  
  // Rate of concentration changes in different compartments    
    double RAST = - Kst * AST;        // MEL oral to the stomach 
    double RAI = Kst * AST - Kint * AI - Ka * AI; // MEL in Intestine 
    double Rcolon = Kint * AI ;       // MEL in Colon 
    double RAO = Ka * AI;
    double Rplas_free = QC * (CV - Cplasma) * Free; // rate of change in amount of chem in plasma compartment of tissues
    double RL = QL * (Cplasma - CVL)* Free + RAO - Kmet * AL ;// rate of change in amount of the chem in liver
    double Rmet = Kmet * AL;          // Rmet the metabolic rate in liver (mg/h)
    double Rurine = Kurine * CVK;
    double RK = QK * (Cplasma - CVK) * Free - Rurine;
    double RM = QM * (Cplasma - CVM) * Free;  // rate of change in amount of the chem in muscle
    double RF = QF * (Cplasma - CVF) * Free;  // rate of change in amount of the chem in fat
    double RRest = QRest * (Cplasma - CVRest) * Free;
    double Ryolk  = Covary * Ky * Wy;          // rate of drug deposition into yolk mg/h        
    double Rwhite = Covary * Kw * AWw;
    double Rovary = Qovary * (Cplasma - CVovary) * Free - Ryolk - Rwhite; // rate of change in amount of the chem in ovary 
    double Regg_ex = Ryolk + Rwhite; // mg/h
   
  
   
   
    dxdt_AST = RAST;  
    dxdt_AI = RAI;
    dxdt_Acolon = Rcolon;    
    dxdt_AAO = RAO;
    dxdt_Aplas_free = Rplas_free;
    dxdt_AL = RL;                    // amount of chemical in liver
    dxdt_Aurine = Rurine;
    dxdt_AK = RK;
    dxdt_AM = RM;
    dxdt_AF = RF;
    dxdt_Aovary = Rovary;
    dxdt_ARest = RRest;
    dxdt_Amet  = Rmet;  
    dxdt_Ayolk = Ryolk;
    dxdt_AWw = Ww;
    dxdt_Awhite   = Rwhite;
    dxdt_Aovary    = Rovary;
    dxdt_Aegg_ex = Regg_ex;
    
    
    
      
     // Equation for the AUC of chemical in tissue compartment 
    dxdt_AUCCV = CV;                // AUC of chemical in vein
    dxdt_AUCCL = CL;                // AUC of chemical in liver compartment 
    dxdt_AUCCK = CK;                // AUC of chemical in kidney compartment  
    dxdt_AUCCM = CM;                // AUC of chemical in muscle compartment
    dxdt_AUCCF = CF;                // AUC of chemical in fat compartment
    dxdt_AUCCovary = Covary;        // AUC of chemical in ovary compartment
    
      
    double Qbal = QC - QL - QK - QM - QF - Qovary - QRest; 
    double Tmass = Aplas_free + AL + AK + Aurine + AM + AF + Aovary + ARest + Amet + Acolon + Aegg_ex + AI + AST;
    double bal   = ADOSE - Tmass;

$TABLE
capture massL = AL;
capture massK = AK;
capture Plasma = CV;
capture Liver  = CL;
capture Kidney = CK;
capture Muscle = CM;
capture Fat  = CF;
capture Ovary  = Covary;
capture Rest = CRest;
capture Bal  = bal;
capture tmass = Tmass;
capture yolk = Cyolk;
capture white = Cwhite;
'



## Build mrgsolve-based PBPK Model

## Load Model
mod <- mcode_cache("pbpk", SolvePBPKOTC) #refer to mcode function in mrgsolve user guide 3.1.2 Inline



pred<- function (pars, tinterval, TDOSE, route, Dose, tlay1 ){
  BW = 2.2
  DOSE = Dose * BW
  
  
  ## Get out of log domain
  pars <- exp(pars)   
  
  
  if (route == "iv") {
    
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "Aplas_free", replicate = FALSE)
    
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
    
    ex <- ev_1+ev_2
    
    
  }
  
  
  
  if (route == "oral") {
    
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval,
                
                addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
    
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval,
                
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    
    ex <- ev_1+ev_2
    
  }
  
  
  
  
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*20, 0.1)
  
  
  
  ## Simulation
  
  out <- mod %>% param(c(tlay1=tlay1))%>%param(pars)%>%
    
    mrgsim_d (data = ex, tgrid = tsamp)
  
  output <- cbind.data.frame(Time = out$time,
                             AL = out$massL,
                             AK = out$massK,
                             CP  = out$Plasma,
                             CL  = out$Liver,
                             CK  = out$Kidney,
                             CM  = out$Muscle,
                             CF  = out$Fat,
                             Cyolk = out$yolk, 
                             Cwhite = out$white, 
                             Bal = out$Bal)
  
  return(output)
  
}

## Input the dataset for plasma calibration and evaluation
data1 <- read.table(file = "oral plasma calibration 1.txt", header=T, skip = 2)
data2 <- read.table(file = "oral plasma calibration 2.txt", header=T, skip = 2)
data3 <- read.table(file = "laying hen multiple dosing plasma calibration.txt", header=T, skip = 3)
data4 <- read.table(file = "oral plasma evaluation 1.txt", header=T, skip = 2)

## Input the dataset for egg yolk and white calibration and evaluation
data5 <- read.table(file = "laying hen multiple dosing egg yolk calibration.txt", header=T, skip = 3)
data6 <- read.table(file = "laying hen multiple dosing egg white calibration.txt", header=T, skip = 3)
data7 <- read.table(file = "laying hen multiple dosing egg yolk evaluation 1.txt", header=T, skip = 3)
data8 <- read.table(file = "oral plasma yolk evaluation 2.txt", header=T, skip = 2)
data9 <- read.table(file = "oral plasma white evaluation 1.txt", header=T, skip = 2)

##Input the dataset for tissue data calibration
data10 <- read.table(file = "Liver study 3 laying hen evaluation.txt", header=T, skip = 2)
data11 <- read.table(file = "Kidney study 3 laying hen evaluation.txt", header=T, skip = 2)
data12 <- read.table(file = "Muscle study 3 laying hen evaluation.txt", header=T, skip = 2)
data13 <- read.table(file = "Fat study 3 laying hen evaluation.txt", header=T, skip = 2)

##Input the dataset for tissue data evaluation
data14 <- read.table(file = "Liver study 4 laying hen evaluation.txt", header=T, skip = 2)
data15 <- read.table(file = "Kidney study 4 laying hen evaluation.txt", header=T, skip = 2)
data16 <- read.table(file = "Muscle study 4 laying hen evaluation.txt", header=T, skip = 2)
data17 <- read.table(file = "Fat study 4 laying hen evaluation.txt", header=T, skip = 2)


##Input the dataset for IV tissue data evaluation
data18 <- read.table(file = "Liver study 5 laying hen evaluation.txt", header=T, skip = 2)
data19 <- read.table(file = "Kidney study 5 laying hen evaluation.txt", header=T, skip = 2)
data20 <- read.table(file = "Muscle study 5 laying hen evaluation.txt", header=T, skip = 2)
data21 <- read.table(file = "Fat study 5 laying hen evaluation.txt", header=T, skip = 2)
data22 <- read.table(file = "Plasma study 5 laying hen evaluation.txt", header=T, skip = 2)

### Input the new dataset for laying hen evaluation
datanew <- read.csv(file = "study 2 group 1.csv", header=T, skip = 2)
datanew1 <- read.csv(file = "study 2 plasma.csv", header=T, skip = 2)
datanew2 <- read.csv(file = "study 2 group 2.csv", header=T, skip = 2)
datanew3 <- read.csv(file = "study 2 group 2 copy.csv", header=T, skip = 2)
datanew4 <- read.csv(file = "study 3 follicle 1.csv", header=T, skip = 2)
datanew5 <- read.csv(file = "study 3 follicle 2.csv", header=T, skip = 2)
datanew6 <- read.csv(file = "study 2 follicle 1.csv", header=T, skip = 2)
datanew7 <- read.csv(file = "study 2 follicle 2.csv", header=T, skip = 2)
datanew8 <- read.csv(file = "study 2 follicle 3.csv", header=T, skip = 2)
datanew9 <- read.csv(file = "study 1.csv", header=T, skip = 2)
datanew10 <- read.csv(file = "study 4 follicle 1.csv", header=T, skip = 2)
datanew11 <- read.csv(file = "study 4 follicle 2.csv", header=T, skip = 2)
datanew12 <- read.csv(file = "study 4 follicle 3.csv", header=T, skip = 2)
datanew13 <- read.csv(file = "study 5 follicle 1.csv", header=T, skip = 2)
datanew14 <- read.csv(file = "study 5 follicle 2.csv", header=T, skip = 2)

### Response letter new dataset
datanew15 <- read.table(file = "new dataset white response letter.txt", header=T, skip = 2)
datanew16 <- read.table(file = "new dataset wyandotte response letter.txt", header=T, skip = 2)


## Transform the dataset format
data1<- cbind.data.frame(Time = data1$Time, CP  = data1$Concentration)
data2<- cbind.data.frame(Time = data2$Time, CP  = data2$Concentration)
data3<- cbind.data.frame(Time = data3$Time, CP  = data3$Concentration)
data4<- cbind.data.frame(Time = data4$Time, CP  = data4$Concentration)

data5<- cbind.data.frame(Time = data5$Time, Cyolk  = data5$Concentration)
data6<- cbind.data.frame(Time = data6$Time, Cwhite  = data6$Concentration)
data7<- cbind.data.frame(Time = data7$Time, Cyolk  = data7$Concentration)
data8<- cbind.data.frame(Time = data8$Time, Cyolk  = data8$Concentration)
data9<- cbind.data.frame(Time = data9$Time, Cwhite  = data9$Concentration)

data10<- cbind.data.frame(Time = data10$Time, CL  = data10$Concentration)
data11<- cbind.data.frame(Time = data11$Time, CK  = data11$Concentration)
data12<- cbind.data.frame(Time = data12$Time, CM  = data12$Concentration)
data13<- cbind.data.frame(Time = data13$Time, CF  = data13$Concentration)

data14<- cbind.data.frame(Time = data14$Time, CL  = data14$Concentration)
data15<- cbind.data.frame(Time = data15$Time, CK  = data15$Concentration)
data16<- cbind.data.frame(Time = data16$Time, CM  = data16$Concentration)
data17<- cbind.data.frame(Time = data17$Time, CF  = data17$Concentration)

data18<- cbind.data.frame(Time = data18$Time, CL  = data18$Concentration)
data19<- cbind.data.frame(Time = data19$Time, CK  = data19$Concentration)
data20<- cbind.data.frame(Time = data20$Time, CM  = data20$Concentration)
data21<- cbind.data.frame(Time = data21$Time, CF  = data21$Concentration)
data22<- cbind.data.frame(Time = data22$Time, CP  = data22$Concentration)


data23<- cbind.data.frame(Time = datanew$Time, Cyolk  = datanew$yolk)
data24<- cbind.data.frame(Time = datanew$Time, Cwhite  = datanew$white)
data25<- cbind.data.frame(Time = datanew1$Time, CP  = datanew1$Concentration)
data26<- cbind.data.frame(Time = datanew2$Time, Cyolk  = datanew2$yolk)
data27<- cbind.data.frame(Time = datanew2$Time, Cwhite  = datanew2$white)
data28<- cbind.data.frame(Time = datanew3$Time, Cyolk  = datanew3$yolk)
data29<- cbind.data.frame(Time = datanew3$Time, Cwhite  = datanew3$white)
data30<- cbind.data.frame(Time = datanew4$Time, Cyolk  = datanew4$follicle)
data31<- cbind.data.frame(Time = datanew5$Time, Cyolk  = datanew5$follicle)
data32<- cbind.data.frame(Time = datanew6$Time, Cyolk  = datanew6$follicle)
data33<- cbind.data.frame(Time = datanew7$Time, Cyolk  = datanew7$follicle)
data34<- cbind.data.frame(Time = datanew8$Time, Cyolk  = datanew8$follicle)
data35<- cbind.data.frame(Time = datanew9$Time, Cyolk  = datanew9$yolk)
data36<- cbind.data.frame(Time = datanew9$Time, Cwhite  = datanew9$white)
data37<- cbind.data.frame(Time = datanew10$Time, Cyolk  = datanew10$follicle)
data38<- cbind.data.frame(Time = datanew11$Time, Cyolk  = datanew11$follicle)
data39<- cbind.data.frame(Time = datanew12$Time, Cyolk  = datanew12$follicle)
data40<- cbind.data.frame(Time = datanew13$Time, Cyolk  = datanew13$follicle)
data41<- cbind.data.frame(Time = datanew14$Time, Cyolk  = datanew14$follicle)

data42<- cbind.data.frame(Time = datanew15$Time, CP  = datanew15$Concentration)
data43<- cbind.data.frame(Time = datanew16$Time, CP  = datanew16$Concentration)


head(data1)
head(data2)
head(data3)
head(data4)
head(data5)
head(data6)
head(data7)
head(data8)
head(data9)
head(data10)
head(data11)
head(data12)
head(data13)
head(data14)
head(data15)
head(data16)
head(data17)
head(data18)
head(data19)
head(data20)
head(data21)
head(data22)






## Create the Cost function
Cost <- function(pars, w) {
  
  
  ## prediction for plasma
  out_A <-  pred (pars,  tinterval = 12, TDOSE = 10, route = 'oral', Dose =1, tlay1 = 12)
  out_B <-  pred (pars,  tinterval = 12, TDOSE = 1, route = 'oral', Dose=1, tlay1 = 12)
  out_C <-  pred (pars,  tinterval = 12, TDOSE = 20, route = 'oral', Dose=1, tlay1 = 12)
  
  
  ## tlength for Cyolk
  t1  <- data5$Time
  tout1 <- list()
  for (i in 1:length(t1)) {
    output1 <- pred(pars, tinterval = 12, TDOSE = 20, route = 'oral', Dose =1, tlay1 = t1[i]/24)
    tout1[[i]] <-output1 %>% filter(Time == t1[i])
  }
  
  yout1 <- do.call(rbind,tout1)%>%select(Time = Time, Cyolk = Cyolk)
  
  ## tlength for Cwhite
  t2  <- data6$Time
  tout2 <- list()
  for (i in 1:length(t2)) {
    output2 <- pred(pars, tinterval = 12, TDOSE = 20, route = 'oral', Dose =1, tlay1 = t2[i]/24)
    tout2[[i]] <-output2 %>% filter(Time == t2[i])
  }
  
  yout2 <- do.call(rbind,tout2)%>%select(Time = Time, Cwhite = Cwhite)  
  
  ## prediction for oral tissue
  out_D <-  pred (pars,  tinterval = 24, TDOSE = 8, route = 'oral', Dose =1, tlay1 = 12)
  
 
  cost<- modCost  (model = out_A, obs = data1, x ="Time", weight = w)
  cost<- modCost  (model = out_B, obs = data2, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_C, obs = data3, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout1, obs = data5, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout2, obs = data6, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_D, obs = data10, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_D, obs = data11, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_D, obs = data12, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_D, obs = data13, x ="Time", weight = w, cost = cost)

  
  return(cost)
  
}

## input parameters
pars <- c(PB = 0.994, 
          Ka = 0.5,
          KurineC = 0.1, 
          PRest =0.479, 
          PK = 0.6183,
          Kint = 0.2,
          KmC = 0.0025,
          Kst = 2,
          Povary = 0.2,
          Kw   = 0.1,
          Ky = 0.09
          )

Cost(log(pars), w="mean")

### Sensitivity analysis
Sns <- sensFun(func = Cost, w = "mean",
               parms = log(pars), varscale = 1)

Sen_1 <- summary(Sns)
Sen_1

## Selected sensitive parameters; 
theta <- pars[abs(Sen_1$Mean) >1.2*mean(abs(Sen_1$Mean))]
theta
pars <-c(
  #PB = 0.994,
  Ka = 0.5,
  KurineC = 0.1, 
  #PRest =0.479, 
  PK = 0.6183,
  Kint = 0.2,
  KmC = 0.0025,
  #Kst = 2,
  Povary = 0.1,
  Kw   = 0.1,
  Ky = 0.09
)


## model fitting

Fit<- modFit(f=Cost, p= log(pars), method ="Nelder-Mead",w="mean", 
             lower = c( log(0.05), log(0.01), log(0.6), log(0.02), log(0.00025), log(0.01), log(0.01), log(0.009)),
             upper = c( log(5),    log(1),    log(1),   log(2),    log(0.025),   log(1),    log(1),    log(0.9)),
             control = nls.lm.control(nprint=1))

summary(Fit)                                 ## Summary of fit 
exp(Fit$par)                                 ## Get the arithmetic value out of the log domain
Cost(Fit$par, w="mean")



# Model plasma calibration results
# Compare the difference between initial and final fitting results
out_Aini <-  pred (log(pars),  tinterval = 12, TDOSE = 10, route = 'oral', Dose =1, tlay1 = 12)
out_Afin <-  pred (Fit$par,  tinterval = 12, TDOSE = 10, route = 'oral', Dose =1, tlay1 = 12)
plot(out_Afin$Time, out_Afin$CP, type = 'l', xlim=c(84, 140), ylim=c(0,9),
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
lines(out_Aini$Time, out_Aini$CP, lty = 2) 
points(data1$Time, data1$CP)

write.csv(out_Afin,'Laying Hen_A calibrate.CSV', row.names = FALSE)


p1 <- ggplot(data=data1, aes(x=Time, y=CP), scale = "Free")+
  geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
  #geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=.1, color="blue", size=0.8) +
  geom_line(data=out_Afin, aes(x=Time, y=CP), size=0.8, colour = "black")+ ggtitle("(A) Plasma, 1 mg/kg oral 10-Dose") +
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 15,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
                                                                        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
                                                                        legend.position = "none")+  #+scale_y_log10()+     
  scale_x_continuous(limits=(c(0,140)))+
  scale_y_continuous(limits=(c(0,9)))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
p1





out_Bini <- pred (log(pars),  tinterval = 12, TDOSE = 1, route = 'oral', Dose=1, tlay1 = 12)
out_Bfin <- pred (Fit$par,  tinterval = 12, TDOSE = 1, route = 'oral', Dose=1, tlay1 = 12)
plot(out_Bfin$Time, out_Bfin$CP, type = 'l', xlim=c(0, 55), ylim=c(0,9),
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
lines(out_Bini$Time, out_Bini$CP, lty = 2)
points(data2$Time, data2$CP)

write.csv(out_Bfin,'Laying Hen_B calibrate.CSV', row.names = FALSE)



p2 <- ggplot(data=data2, aes(x=Time, y=CP), scale = "Free")+
  geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
  #geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=.1, color="blue", size=0.8) +
  geom_line(data=out_Bfin, aes(x=Time, y=CP), size=0.8, colour = "black")+ ggtitle("(B) Plasma, 1 mg/kg oral single dose") +
  
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 15,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
                                                                        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
                                                                        legend.position = "none")+  #+scale_y_log10()+     
  scale_x_continuous(limits=(c(0,60)))+
  scale_y_continuous(limits=(c(0,9)))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
p2




out_Cini <- pred (log(pars),  tinterval = 12, TDOSE = 20, route = 'oral', Dose=1, tlay1 = 12)
out_Cfin <- pred (Fit$par,  tinterval = 12, TDOSE = 20, route = 'oral', Dose=1, tlay1 = 12)
plot(out_Cfin$Time, out_Cfin$CP, type = 'l', xlim=c(200, 300), ylim=c(0,9),
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
lines(out_Cini$Time, out_Cini$CP, lty = 2)
points(data3$Time, data3$CP)




write.csv(out_Cfin,'Laying Hen_C calibrate.CSV', row.names = FALSE)



p3 <- ggplot(data=data3, aes(x=Time, y=CP), scale = "Free")+
  geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
  #geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=.1, color="blue", size=0.8) +
  geom_line(data=out_Cfin, aes(x=Time, y=CP), size=0.8, colour = "black")+ ggtitle("(C) Plasma, 1 mg/kg oral 20-Dose") +
  
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 15,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
                                                                        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
                                                                        legend.position = "none")+  #+scale_y_log10()+     
  scale_x_continuous(limits=(c(130,300)))+
  scale_y_continuous(limits=(c(0,9)))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
p3



#Model plasma evaluation
out_D <-  pred (Fit$par,  tinterval = 24, TDOSE = 8, route = 'oral', Dose =1, tlay1 = 12)
plot(out_D$Time, out_D$CP, type = 'l', xlim=c(50, 300), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data4$Time, data4$CP)


#Evaluation for study 2 plasma
plot(out_D$Time, out_D$CP, type = 'l', xlim=c(0, 300), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data25$Time, data25$CP)


write.csv(out_D,'Laying Hen_A evaluate.CSV', row.names = FALSE)



pa <- ggplot(data=data4, aes(x=Time, y=CP), scale = "Free")+
  geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
  #geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=.1, color="blue", size=0.8) +
  geom_line(data=out_D, aes(x=Time, y=CP), size=0.8, colour = "black")+ ggtitle("(A) Plasma, 1 mg/kg oral 8-Dose") +
  
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 15,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
                                                                        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
                                                                        legend.position = "none")+  #+scale_y_log10()+     
  scale_x_continuous(limits=(c(50,300)))+
  scale_y_continuous(limits=(c(0,7)))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
pa


out_D <-  pred (Fit$par,  tinterval = 24, TDOSE = 8, route = 'oral', Dose =1, tlay1 = 12)
plot(out_D$Time, out_D$CP, type = 'l', xlim=c(-10, 300), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data4$Time, data4$CP)





#Model egg yolk and egg white calibration results
# Compare the difference between initial and final fitting results

t <- c((1:20)*24)
toutini <- list()
toutfin <- list()

for (i in 1:length(t)) {
  outputini <- pred(log(pars), tinterval = 12, TDOSE = 20, route = 'oral', Dose =1, tlay1 = t[i]/24)
  outputfin <- pred(Fit$par, tinterval = 12, TDOSE = 20, route = 'oral', Dose =1, tlay1 = t[i]/24)
  toutini[[i]] <-outputini %>% filter(Time == t[i])
  toutfin[[i]] <-outputfin %>% filter(Time == t[i])
}

youtini <- do.call(rbind,toutini)%>%select(Time = Time, Cyolk = Cyolk, Cwhite = Cwhite)
youtfin <- do.call(rbind,toutfin)%>%select(Time = Time, Cyolk = Cyolk, Cwhite = Cwhite)


# Differences between initial and final egg yolk calibration results
plot(youtfin$Time, youtfin$Cyolk, xlab = "time", ylab = "Cyolk" ,
     xlim = c(0,450), ylim = c(0,0.7), type = 'l')
lines(youtini$Time, youtini$Cyolk, lty = 2)
points(data5$Time, data5$Cyolk)



write.csv(youtfin,'Laying Hen_D & E yolk calibrate.CSV', row.names = FALSE)


plot(youtfin$Time, youtfin$Cyolk, xlab = "time", ylab = "Cyolk" ,
     xlim = c(0,450), ylim = c(0,0.16), type = 'l')
points(data35$Time, data35$Cyolk)


plot(youtfin$Time, youtfin$Cwhite, xlab = "time", ylab = "Cwhite" ,
     xlim = c(0,450), ylim = c(0,0.05), type = 'l')
points(data36$Time, data36$Cwhite)


plot(youtfin$Time, youtfin$Cyolk, xlab = "time", ylab = "Cyolk" ,
     xlim = c(0,450), ylim = c(0,0.25), type = 'l')
points(data37$Time, data37$Cyolk)


plot(youtfin$Time, youtfin$Cyolk, xlab = "time", ylab = "Cyolk" ,
     xlim = c(0,450), ylim = c(0,0.25), type = 'l')
points(data38$Time, data38$Cyolk)


plot(youtfin$Time, youtfin$Cyolk, xlab = "time", ylab = "Cyolk" ,
     xlim = c(0,450), ylim = c(0,0.27), type = 'l')
points(data39$Time, data39$Cyolk)


p4 <- ggplot(data=data5, aes(x=Time, y=Cyolk), scale = "Free")+
  geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
  #geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=.1, color="blue", size=0.8) +
  geom_line(data=youtfin, aes(x=Time, y=Cyolk), size=0.8, colour = "black")+ ggtitle("(D) Egg yolk, 1 mg/kg oral 20-Dose") +
  
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 15,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
                                                                        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
                                                                        legend.position = "none")+  #+scale_y_log10()+     
  scale_x_continuous(limits=(c(0,450)))+
  scale_y_continuous(limits=(c(0,0.23)))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
p4





# Differences between initial and final egg white calibration results
plot(youtfin$Time, youtfin$Cwhite, xlab = "time", ylab = "Cwhite" ,
     xlim = c(0,450), ylim = c(0,0.16), 
     type = 'l')
lines(youtini$Time, youtini$Cwhite, lty = 2)
points(data6$Time, data6$Cwhite)


write.csv(youtfin,'Laying Hen_D&E yolk calibrate.CSV', row.names = FALSE)



p5 <- ggplot(data=data6, aes(x=Time, y=Cwhite), scale = "Free")+
  geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
  #geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=.1, color="blue", size=0.8) +
  geom_line(data=youtfin, aes(x=Time, y=Cwhite), size=0.8, colour = "black")+ ggtitle("(E) Egg white, 1 mg/kg oral 20-Dose") +
  
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 15,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
                                                                        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
                                                                        legend.position = "none")+  #+scale_y_log10()+     
  scale_x_continuous(limits=(c(0,450)))+
  scale_y_continuous(limits=(c(0,0.05)))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
p5




#plot the graph combine all of the Cyolk and Cwhite final results
plot(youtfin$Time, youtfin$Cyolk, xlab = "time", ylab = "Cyolk" ,
     xlim = c(0,450), ylim = c(0,0.2), type = 'l')
lines(youtfin$Time, youtfin$Cwhite, xlab = "time", ylab = "Cwhite" , lty = "dashed")
points(data5$Time, data5$Cyolk)
points(data6$Time, data6$Cwhite)


#Model egg yolk and egg white evaluation
# plot egg yolk evaluation 1
t <- c((1:20)*24)
tout3 <- list()

for (i in 1:length(t)) {
  outputf <- pred(Fit$par, tinterval = 12, TDOSE = 10, route = 'oral', Dose =1, tlay1 = t[i]/24)
  tout3[[i]] <-outputf %>% filter(Time == t[i])
}

yout3 <- do.call(rbind,tout3)%>%select(Time = Time, Cyolk = Cyolk)


plot(yout3$Time, yout3$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l')
points(data7$Time, data7$Cyolk)



write.csv(yout3,'Laying Hen_B yolk evaluate.CSV', row.names = FALSE)



pb <- ggplot(data=data7, aes(x=Time, y=Cyolk), scale = "Free")+
  geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
  #geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=.1, color="blue", size=0.8) +
  geom_line(data=yout3, aes(x=Time, y=Cyolk), size=0.8, colour = "black")+ ggtitle("(B) Egg yolk, 1 mg/kg oral 10-Dose") +
  
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 15,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
                                                                        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
                                                                        legend.position = "none")+  #+scale_y_log10()+     
  scale_x_continuous(limits=(c(0,450)))+
  scale_y_continuous(limits=(c(0,0.16)))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
pb


# plot egg yolk evaluation 2
t <- c((1:20)*24)
tout4 <- list()

for (i in 1:length(t)) {
  output4 <- pred(Fit$par, tinterval = 24, TDOSE = 8, route = 'oral', Dose =1, tlay1 = t[i]/24)
  tout4[[i]] <-output4 %>% filter(Time == t[i])
}

yout4 <- do.call(rbind,tout4)%>%select(Time = Time, Cyolk = Cyolk)

plot(yout4$Time, yout4$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l')
points(data8$Time, data8$Cyolk)


# Evaluation for study 2 group 1 egg yolk
plot(yout4$Time, yout4$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l', ylim = c(0,0.09))
points(data23$Time, data23$Cyolk)

write.csv(yout4,'study 2 group 1 yolk evaluation.CSV', row.names = FALSE)

# Evaluaion for study 2 group 2 egg yolk
plot(yout4$Time, yout4$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l')
points(data26$Time, data26$Cyolk)

# Evaluation for study 
plot(yout4$Time, yout4$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l')
points(data28$Time, data28$Cyolk)


plot(yout4$Time, yout4$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l')
points(data30$Time, data30$Cyolk)



plot(yout4$Time, yout4$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l')
points(data31$Time, data31$Cyolk)



plot(yout4$Time, yout4$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l')
points(data32$Time, data32$Cyolk)




plot(yout4$Time, yout4$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l')
points(data33$Time, data33$Cyolk)



plot(yout4$Time, yout4$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l')
points(data34$Time, data34$Cyolk)






write.csv(yout4,'Laying Hen_C yolk evaluate.CSV', row.names = FALSE)





pc <- ggplot(data=data8, aes(x=Time, y=Cyolk), scale = "Free")+
  geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
  #geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=.1, color="blue", size=0.8) +
  geom_line(data=yout4, aes(x=Time, y=Cyolk), size=0.8, colour = "black")+ ggtitle("(C) Egg yolk, 1 mg/kg oral 8-Dose") +
  
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 15,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
                                                                        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
                                                                        legend.position = "none")+  #+scale_y_log10()+     
  scale_x_continuous(limits=(c(0,450)))+
  scale_y_continuous(limits=(c(0,0.10)))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
pc



# plot egg white evaluation 1
t <- c((1:20)*24)
tout5 <- list()

for (i in 1:length(t)) {
  output5 <- pred(Fit$par, tinterval = 24, TDOSE = 8, route = 'oral', Dose =1, tlay1 = t[i]/24)
  tout5[[i]] <-output5 %>% filter(Time == t[i])
}

yout5 <- do.call(rbind,tout5)%>%select(Time = Time, Cwhite = Cwhite)

plot(yout5$Time, yout5$Cwhite, xlab = "time", ylab = "Cwhite" , type = 'l', ylim = c(0, 0.017))
points(data9$Time, data9$Cwhite)


# Evaluation for study 2 group 1 egg white
plot(yout5$Time, yout5$Cwhite, xlab = "time", ylab = "Cwhite" , type = 'l', xlim = c(0, 300), ylim = c(0, 0.04))
points(data24$Time, data24$Cwhite)

write.csv(yout5,'study2 group1 egg white evaluate.CSV', row.names = FALSE)



# Evaluation for study 2 group 2 egg white
plot(yout5$Time, yout5$Cwhite, xlab = "time", ylab = "Cwhite" , type = 'l', xlim = c(0, 300), ylim = c(0, 0.015))
points(data27$Time, data27$Cwhite)



plot(yout5$Time, yout5$Cwhite, xlab = "time", ylab = "Cwhite" , type = 'l', xlim = c(0, 300), ylim = c(0, 0.015))
points(data29$Time, data29$Cwhite)




write.csv(yout5,'Laying Hen_D white evaluate.CSV', row.names = FALSE)





pd <- ggplot(data=data9, aes(x=Time, y=Cwhite), scale = "Free")+
  geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
  #geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=.1, color="blue", size=0.8) +
  geom_line(data=yout5, aes(x=Time, y=Cwhite), size=0.8, colour = "black")+ ggtitle("(D) Egg yolk, 1 mg/kg oral 8-Dose") +
  
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 15,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
                                                                        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
                                                                        legend.position = "none")+  #+scale_y_log10()+     
  scale_x_continuous(limits=(c(0,450)))+
  scale_y_continuous(limits=(c(0,0.025)))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
pd





## combine evaluation of Cyolk and Cwhite together
plot(yout4$Time, yout4$Cyolk, xlab = "time", ylab = "Cyolk" ,
     xlim = c(0,450), type = 'l', col = "red")
lines(yout5$Time, yout5$Cwhite, xlab = "time", ylab = "Cwhite" , lty = "dashed")
points(data8$Time, data8$Cyolk, col = "red")
points(data9$Time, data9$Cwhite)



#Model egg yolk and egg white evaluation for study 5 
# plot egg yolk evaluation 1
t <- c((1:20)*24)
tout3 <- list()

for (i in 1:length(t)) {
  outputf <- pred(Fit$par, tinterval = 12, TDOSE = 10, route = 'oral', Dose =1, tlay1 = t[i]/24)
  tout3[[i]] <-outputf %>% filter(Time == t[i])
}

yout3 <- do.call(rbind,tout3)%>%select(Time = Time, Cyolk = Cyolk, Cwhite = Cwhite)


plot(yout3$Time, yout3$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l')
points(data7$Time, data7$Cyolk)




##Model Liver calibration study 3
plot(out_D$Time, out_D$CL, type = 'l', xlim=c(50, 300), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data10$Time, data10$CL)


#Model Kidney calibration study 3
plot(out_D$Time, out_D$CK, type = 'l', xlim=c(50, 300),ylim=c(0, 2), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data11$Time, data11$CK)

#Model Muscle calibration study 3
plot(out_D$Time, out_D$CM, type = 'l', xlim=c(50, 300), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data12$Time, data12$CM)

#Model Fat calibration study 3
plot(out_D$Time, out_D$CF, type = 'l', xlim=c(50, 300), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data13$Time, data13$CF)


##Model Liver evaluation study 4
plot(out_Cfin$Time, out_Cfin$CL, type = 'l', xlim=c(50, 300), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data14$Time, data14$CL)


#Model Kidney evaluation study 4
plot(out_Cfin$Time, out_Cfin$CK, type = 'l', xlim=c(50, 300),ylim=c(0, 2), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data15$Time, data15$CK)

#Model Muscle evaluation study 4
plot(out_Cfin$Time, out_Cfin$CM, type = 'l', xlim=c(50, 300), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data16$Time, data16$CM)

#Model Fat evaluation study 4
plot(out_Cfin$Time, out_Cfin$CF, type = 'l', xlim=c(50, 300), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data17$Time, data17$CF)


write.csv(out_Cfin,'Laying Hen_study 4 tissue evaluate.CSV', row.names = FALSE)

##Model Liver evaluation study 5
out_E <-  pred(Fit$par,  tinterval = 1, TDOSE = 1, route = 'iv', Dose =1, tlay1 = 12)
plot(out_E$Time, out_E$CL, type = 'l', xlim=c(0, 30), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data18$Time, data18$CL)


#Model Kidney evaluation study 5
plot(out_E$Time, out_E$CK, type = 'l', xlim=c(0, 30),#ylim=c(0, 0.1), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data19$Time, data19$CK)

#Model Muscle evaluation study 5
plot(out_E$Time, out_E$CM, type = 'l', xlim=c(0, 30), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data20$Time, data20$CM)

#Model Fat evaluation study 5
plot(out_E$Time, out_E$CF, type = 'l', xlim=c(0, 30), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data21$Time, data21$CF)

#Model Plasma evaluation study 5
plot(out_E$Time, out_E$CP, type = 'l', xlim=c(0, 30), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data22$Time, data22$CP)


write.csv(out_E,'Laying Hen study 5 evaluate.CSV', row.names = FALSE)

##Model New dataset response letter
out_E <-  pred(Fit$par,  tinterval = 1, TDOSE = 1, route = 'oral', Dose =1, tlay1 = 12)
plot(out_E$Time, out_E$CP, type = 'l', xlim=c(0, 30), ylim=c(0, 10), 
     xlab = "Time(h)", ylab = "Concentration(mg/Kg)")
points(data42$Time, data42$CP)
points(data43$Time, data43$CP)
write.csv(out_E,'newdataset from response letter.CSV', row.names = FALSE)



#Model egg follicle evaluation for study 5 iv injection 
# plot egg follicle evaluation 
t <- c((1:20)*24)
tout6 <- list()

for (i in 1:length(t)) {
  outputf <- pred(Fit$par, tinterval = 1, TDOSE = 1, route = 'iv', Dose =1, tlay1 = t[i]/24)
  tout6[[i]] <-outputf %>% filter(Time == t[i])
}

yout6 <- do.call(rbind,tout6)%>%select(Time = Time, Cyolk = Cyolk, Cwhite = Cwhite)


plot(yout6$Time, yout6$Cyolk, xlab = "time", ylab = "Cyolk" , xlim = c(0,200), ylim = c(0,0.10),
     type = 'l' )
points(data40$Time, data40$Cyolk)


plot(yout6$Time, yout6$Cyolk, xlab = "time", ylab = "Cyolk" , xlim = c(0,200), ylim = c(0,0.10),
     type = 'l' )
points(data41$Time, data41$Cyolk)


## Check mass balance
out_A <-  pred (log(pars),  tinterval = 12, TDOSE = 10, route = 'oral', Dose =1, tlay1 = 12)%>%filter(Time >0)
Bal = out_A$Bal
plot(out_A$Time, Bal, col = 'red')


## Regression analysis for laying hen model calibration
Fit_laying<- lm(obs ~ mod, data = Cost(Fit$par, w="mean")$residuals)
summary(Fit_laying)

write.csv(Cost(Fit$par, w="mean")$residuals,'Regression analysis laying hen calibration.CSV', row.names = FALSE)


p6 <- ggplot (Cost(Fit$par, w="mean")$residuals, aes(obs, mod)) +
  geom_point   (size = 2) +
  geom_abline  (intercept = 0,
                slope     = 1,
                color     ="black",size = 1) +
  annotation_logticks() +
  scale_y_continuous(limits = c(0,7.5))+
  scale_x_continuous(limits = c(0,7.5)) +
  theme_bw() + labs (x = "Observed value (mg/L)", y = "Predicted value (mg/L)")+
  theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
        = element_text(size = 15,  color = "black", face="bold", family="A"),
        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
        legend.position = "none")+  #+scale_y_log10()+ 
  ggtitle("(F) Laying hen calibration")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())

p6



## Regression analysis for all of laying hen model evaluation except for IV dataset
#Generate Cost function for all of the evaluation results except for IV dataset


Cost1 <- function(pars, w) {
  
  cost<- modCost  (model = out_D, obs = data4, x ="Time", weight = w)
  cost<- modCost  (model = yout3, obs = data7, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout4, obs = data8, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout5, obs = data9, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_Cfin, obs = data14, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_Cfin, obs = data15, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_Cfin, obs = data16, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_Cfin, obs = data17, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = out_E, obs = data18, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = out_E, obs = data19, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = out_E, obs = data20, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = out_E, obs = data21, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = out_E, obs = data22, x ="Time", weight = w, cost = cost)
  
  cost<- modCost  (model = yout4, obs = data23, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout5, obs = data24, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_D, obs = data25, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout4, obs = data26, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout5, obs = data27, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout4, obs = data28, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout5, obs = data29, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout4, obs = data30, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout4, obs = data31, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout4, obs = data32, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout4, obs = data33, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout4, obs = data34, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = youtfin, obs = data35, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = youtfin, obs = data36, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = youtfin, obs = data37, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = youtfin, obs = data38, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = youtfin, obs = data39, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout6, obs = data40, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = yout6, obs = data41, x ="Time", weight = w, cost = cost)
  
  
  
  
  return(cost)
  
}

Cost1(Fit$par, w="mean")

## Regression analysis for all laying hen model evaluation for IV dataset

Fit_laying_eva1<- lm(obs ~ mod, data = Cost1(Fit$par, w="mean")$residuals)
summary(Fit_laying_eva1)
eva1 <- ggplot (Cost4(Fit$par, w="mean")$residuals, aes(obs, mod)) +
  geom_point   (size = 2) +
  geom_abline  (intercept = 0,
                slope     = 1,
                color     ="black",size = 1) +
  annotation_logticks() +
  scale_y_continuous(limits = c(0,15))+
  scale_x_continuous(limits = c(0,15)) +
  theme_bw() + labs (x = "Observed value (mg/L)", y = "Predicted value (mg/L)")+
  theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
        = element_text(size = 15,  color = "black", face="bold", family="A"),
        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
        legend.position = "none")+  #+scale_y_log10()+ 
  ggtitle("(E) Laying hen evaluation 1")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())

eva1

write.csv(Cost1(Fit$par, w="mean")$residuals,'Regression analysis laying hen evaluation 1.CSV', row.names = FALSE)


## Regression analysis for all of laying hen model evaluation only for IV dataset
#Generate Cost function for all of the evaluation results only for IV dataset,18-22



Cost2 <- function(pars, w) {
  
  # cost<- modCost  (model = out_D, obs = data4, x ="Time", weight = w)
  # cost<- modCost  (model = yout3, obs = data7, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout4, obs = data8, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout5, obs = data9, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = out_Cfin, obs = data14, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = out_Cfin, obs = data15, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = out_Cfin, obs = data16, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = out_Cfin, obs = data17, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_E, obs = data18, x ="Time", weight = w)
  cost<- modCost  (model = out_E, obs = data19, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_E, obs = data20, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_E, obs = data21, x ="Time", weight = w, cost = cost)
  cost<- modCost  (model = out_E, obs = data22, x ="Time", weight = w, cost = cost)
  
  # cost<- modCost  (model = yout4, obs = data23, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout5, obs = data24, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = out_D, obs = data25, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout4, obs = data26, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout5, obs = data27, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout4, obs = data28, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout5, obs = data29, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout4, obs = data30, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout4, obs = data31, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout4, obs = data32, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout4, obs = data33, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout4, obs = data34, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = youtfin, obs = data35, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = youtfin, obs = data36, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = youtfin, obs = data37, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = youtfin, obs = data38, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = youtfin, obs = data39, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout6, obs = data40, x ="Time", weight = w, cost = cost)
  # cost<- modCost  (model = yout6, obs = data41, x ="Time", weight = w, cost = cost)
  
  
  
  
  return(cost)
  
}

Cost2(Fit$par, w="mean")
summary(Cost2)
## Regression analysis for all laying hen model evaluation except for IV dataset

Fit_laying_eva2<- lm(obs ~ mod, data = Cost2(Fit$par, w="mean")$residuals)
summary(Fit_laying_eva2)
eva2 <- ggplot (Cost2(Fit$par, w="mean")$residuals, aes(obs, mod)) +
  geom_point   (size = 2) +
  geom_abline  (intercept = 0,
                slope     = 1,
                color     ="black",size = 1) +
  annotation_logticks() +
  scale_y_continuous(limits = c(0,15))+
  scale_x_continuous(limits = c(0,15)) +
  theme_bw() + labs (x = "Observed value (mg/L)", y = "Predicted value (mg/L)")+
  theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
        = element_text(size = 15,  color = "black", face="bold", family="A"),
        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
        legend.position = "none")+  #+scale_y_log10()+ 
  ggtitle("(E) Laying hen evaluation 2")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())

eva2

write.csv(Cost2(Fit$par, w="mean")$residuals,'Regression analysis laying hen evaluation 2.CSV', row.names = FALSE)



## Regression analysis for all new dataset from response letter


Cost3 <- function(pars, w) {

  cost<- modCost  (model = out_E, obs = data42, x ="Time", weight = w)
  cost<- modCost  (model = out_E, obs = data43, x ="Time", weight = w, cost = cost)
  
  return(cost)
  
}
Cost3(Fit$par, w="mean")
summary(Cost3)

Fit_laying_eva3<- lm(obs ~ mod, data = Cost3(Fit$par, w="mean")$residuals)
summary(Fit_laying_eva2)
eva3 <- ggplot (Cost3(Fit$par, w="mean")$residuals, aes(obs, mod)) +
  geom_point   (size = 2) +
  geom_abline  (intercept = 0,
                slope     = 1,
                color     ="black",size = 1) +
  annotation_logticks() +
  scale_y_continuous(limits = c(0,15))+
  scale_x_continuous(limits = c(0,15)) +
  theme_bw() + labs (x = "Observed value (mg/L)", y = "Predicted value (mg/L)")+
  theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
        = element_text(size = 15,  color = "black", face="bold", family="A"),
        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
        legend.position = "none")+  #+scale_y_log10()+ 
  ggtitle("(E) Laying hen evaluation 3")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())

eva3








