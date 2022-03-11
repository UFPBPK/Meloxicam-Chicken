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


SolvePBPK <- '
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
QCC     : 10.17     :                 Cardiac output (L/h/kg), Table 13, Wang et al. 2021
QLC     : 0.2526    :                 Fraction of flow to the liver, Table 16, Wang et al. 2021
QKC     : 0.2012    :                 Fraction of flow to the kidneys, Table 16, Wang et al. 2021
QFC     : 0.015     :                 Fraction of flow to the fat, Zeng et al. 2019
QMC     : 0.35      :                 Fraction of flow to the muscle, Zeng et al. 2019


// Tissue  volumes
BW      : 2.2       :                 Body weight(kg)
VLC     : 0.0204    :                 Fractional liver tissue, Table 2, Wang et al. 2021
VKC     : 0.0057    :                 Fractional kidney tissue, Table 2, Wang et al. 2021
VFC     : 0.05      :                 Fractional fat tissue, Zeng et al. 2019
VMC     : 0.4015    :                 Fractional muscle tissue, Table 2, Wang et al. 2021
VbloodC : 0.0483    :                 Blood volume, fraction of BW , Table 2, Wang et al. 2021 
Htc     : 0.307     :                 Hematocrit volume, fraction of BW, Table 22, Wang et al. 2021

  
// Mass Transfer Parameters (Chemical-specific parameters)
// Partition coefficients(PC, tissue:plasma)
PL      : 0.6116    :                 Liver plasma PC, calculated using AUCtissue:AUCplasma method from Dutch et al. 2021
PK      : 0.6183    :                 Kidney plasma PC, calculated using AUCtissue:AUCplasma method from  Dutch et al. 2021
PM      : 0.0456    :                 Muscle plasma PC, breast muscle, calculated using AUCtissue:AUCplasma method from  Dutch et al. 2021 
PF      : 0.0174    :                 Fat plasma PC, calculated using AUCtissue:AUCplasma method from  Dutch et al. 2021
PRest   : 0.479     :                 Rest of the body tissues:plasma PC, Table 2, Li M et al. 2019
PB      : 0.994     :                 Percentage of drug bound to plasma proteins (DrugBank, 2020)

$MAIN


// Cardiac output and blood flows to tissues(L/h)
double QC = QCC * BW;                  // Cardiac output
double QL = QLC * QC;                  // Liver
double QK = QKC * QC;                  // Kidney
double QF = QFC * QC;                  // Fat
double QM = QMC * QC;                  // Muscle
double QRest = QC - QK - QL - QF - QM; // Rest of the tissues


//  Tissue volumes (L)
double VL = VLC * BW;                  // Liver
double VK = VKC * BW;                  // Kidney
double VF = VFC * BW;                  // Fat
double VM = VMC * BW;                  // Muscle
double Vblood = VbloodC*BW;            // Blood
double Vplasma = VplasmaC * BW;        // Plasma
double VRest = 1 * BW - VL - VK - VF - VM - Vblood; // Rest of the tissues
double VplasmaC = (1-Htc) * VbloodC;   // Plasma volume, fraction of BW, Table 22, Wang et al. 2021
double Kmet = KmC * BW;	               // Metabolic rate constant adjusted by body weight
double Kurine = KurineC * BW;          // # L/h Urinary elimination rate

double Free = 1-PB;                    //Free drug percentage
      

$CMT AIV ADOSE AST AI Acolon AAO Aplas_free AL AK Aurine Amet AM  AUCCV  AUCCL AUCCK  AUCCM  AUCCF AF ARest 


             
$ODE

dxdt_ADOSE=0;

   // The changing rate of the amount of dose via oral

 
 // Concentration of the chemical in vein and tissue compartments
      double CVL = AL/(VL * PL);      // Concentration in liver / PC of plasma liver
      double CVK = AK/(VK * PK);      // Concentration in Kidney / PC of plasma kidney
      double CVF = AF/(VF * PF);   
      double CVRest = ARest/(VRest * PRest); // Concentration in Rest  / PC of plasma rest tissue
      double CVM = AM/(VM * PM);      // perfusion limited model.
      double CV = ((QL * CVL  + QK  * CVK  + QF * CVF  + QM * CVM  + QRest * CVRest)/QC);
      double Cplas_free = Aplas_free/Vplasma;  // Free drug concentration in plasma = amount of chem in plasma / volume of plasma
      double Cplasma = Cplas_free/Free ;
      double CL = AL/VL;              // Concentration in liver
      double CK = AK/VK;              // Concentration in kidney
      double CM = AM/VM;              // Concentration in muscle
      double CF = AF/VF;              // Concentration in fat for perfusion limited model
      double CRest = ARest/VRest;  
  
  
  // Rate of concentration changes in different compartments    
    double RAST = - Kst * AST;        // MEL oral to the stomach 
    double RAI = Kst * AST - Kint * AI - Ka * AI; // MEL in Intestine 
    double Rcolon = Kint * AI ;       // MEL in Intestine 
    double RAO = Ka * AI;
    double Rplas_free = QC * (CV - Cplasma) * Free; // rate of change in amount of chem in plasma compartment of tissues
    double RL = QL * (Cplasma - CVL)* Free + RAO - Kmet * AL ;// rate of change in amount of the chem in liver
    double Rmet = Kmet * AL;          // Rmet the metabolic rate in liver (mg/h)
    double Rurine = Kurine * CVK;
    double RK = QK * (Cplasma - CVK) * Free - Rurine;
    double RM = QM * (Cplasma - CVM) * Free;  // rate of change in amount of the chem in muscle
    double RF = QF * (Cplasma - CVF) * Free;  // rate of change in amount of the chem in muscle
    double RRest = QRest * (Cplasma - CVRest) * Free;
             
   
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
    dxdt_ARest = RRest;
    dxdt_Amet  = Rmet;  
      
      
     // Equation for the AUC of chemical in tissue compartment 
    dxdt_AUCCV = CV;                // AUC of chemical in vein
    dxdt_AUCCL = CL;                // AUC of chemical in liver compartment 
    dxdt_AUCCK = CK;                // AUC of chemical in kidney compartment  
    dxdt_AUCCM = CM;                // AUC of chemical in muscle compartment
    dxdt_AUCCF = CF;                // AUC of chemical in fat compartment
      
    double Qbal = QC - QL - QK - QM - QF - QRest; 
    double Tmass = Aplas_free + AL + AK + Aurine + AM + AF + ARest + Amet + Acolon;
    double bal   = ADOSE  - Tmass;

$TABLE
capture Plasma = CV;
capture Liver  = CL;
capture Kidney = CK;
capture Muscle = CM;
capture Fat  = CF;
capture Rest = CRest;
capture Bal  = bal;
capture tmass = Tmass;
'



## Build mrgsolve-based PBPK Model

## Load Model
mod <- mcode_cache("pbpk", SolvePBPK) #refer to mcode function in mrgsolve user guide 3.1.2 Inline


pred<- function (pars){
  BW         = 2.2
  TDOSE      = 1
  tinterval  = 24
  DOSE       = 0.5*BW

 
  ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "Aplas_free", replicate = FALSE)
  
  ev_2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, tinf = 0.01,
              
              addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
  
  ex <- ev_1 +ev_2
  


  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*20, 0.01)



## Simulation

  out <- mod %>% param(pars)%>%
  
        mrgsim_d (data = ex, tgrid = tsamp)

  output <- cbind.data.frame(Time = out$time,
                             CL   = out$Liver,
                             CK   = out$Kidney,
                             CP   = out$Plasma)
  
  return(output)

}

# Set working directory
# Read in observed data for calibration
data <- read.table(file="for calibration.txt", header=T, skip = 1)
# Read in observed data for evaluation
data1 <- read.table(file="for evaluation.txt", header=T, skip = 1)


# Re-organize the observe data
data_obs <- cbind.data.frame (Time = data$Time,
                  CP = data$Concentration) 

data_obs1 <- cbind.data.frame(Time = data1$Time,
                   CP = data1$Concentration) 




## Cost function
Cost <- function (pars) {
  
  #prediction function with exposure scenario
  out <- pred(pars)
  cost <- modCost(out, data_obs, x = "Time", weight = "mean")
  
  return(cost)
}


# input the pars value
pars = c(KurineC=0.1, 
         PRest = 0.479
         #KmC = 0.0025 
         #Kint = 0.2
)
Cost(pars)


# # Model fitting
Fit <- modFit(f=Cost, p= pars, method ="Marq",
       control = nls.lm.control(nprint=1), lower = c(0.01, 0.0479), upper = c(1, 4.79))



summary(Fit)

## Transform the pars
pars[names(Fit$par)]<-Fit$par
newpars<-pars
outnew <- pred(newpars)

##  Generate a the model with new parameters fitted
write.csv(outnew,'Broiler Chicken.CSV', row.names = FALSE)




## ggplot for calibration
windowsFonts(A = windowsFont("Times New Roman"))
p1 <- ggplot(data=data_obs, aes(x=Time, y=CP), scale = "Free")+
  geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
  #geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=.1, color="blue", size=0.8) +
  geom_line(data=outnew, aes(x=Time, y=CP), size=0.8, colour = "black")+ ggtitle("(A) Plasma, 0.5 mg/kg IV dose") +

  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
                                                                       = element_text(size = 15,  color = "black", face="bold", family="A"),
                                                                       axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
                                                                       strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
                                                                       #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
                                                                       legend.position = "none")+  #+scale_y_log10()+     
  scale_x_continuous(limits=(c(0,24)))+
  scale_y_continuous(limits=(c(0,11)))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
p1


## Regression analysis for chicken model for calibration
Evalu_plasma<- lm(obs ~ mod, data = Cost(Fit$par)$residuals)
summary(Evalu_plasma)
p2 <- ggplot (Cost(Fit$par)$residuals, aes(obs, mod)) +
  geom_point   (size = 2) +
  geom_abline  (intercept = 0,
                slope     = 1,
                color     ="black",size = 1) +
  annotation_logticks() +
  scale_y_continuous(limits = c(0,10))+
  scale_x_continuous(limits = c(0,10)) +
  theme_bw() + labs (x = "Observed value (mg/L)", y = "Predicted value (mg/L)")+
  theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
        = element_text(size = 15,  color = "black", face="bold", family="A"),
        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
        legend.position = "none")+  #+scale_y_log10()+ 
   ggtitle("(B) Broiler chicken calibration")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line())
p2



## ggplot for evaluation
p3 <- ggplot(data=data_obs1, aes(x=Time, y=CP), scale = "Free")+
  geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
  #geom_errorbar(aes(ymin=Plasma-Plasma_sd, ymax=Plasma+Plasma_sd), width=.1, color="blue", size=0.8) +
  geom_line(data=outnew, aes(x=Time, y=CP), size=0.8, colour = "black")+ ggtitle("(C) Plasma, 0.5 mg/kg IV dose") +
  
  ylab("Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 15,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
                                                                        #strip.text.y = element_text(size=20, face="bold", color="black", family="A"),
                                                                        legend.position = "none")+    
  scale_x_continuous(limits=(c(0,24)))+
  scale_y_continuous(limits=(c(0,11)))+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
p3


## Regression analysis for chicken model for evaluation
# Calculate the model cost
out_1 <- pred(pars)
Cost1 <- modCost(out_1, data_obs1, x = "Time", weight = "mean") 
Cost1
Evalu_plasma<- lm(obs ~ mod, data = Cost1$residuals)
summary(Evalu_plasma)
p4 <- ggplot (Cost1$residuals, aes(obs, mod)) +
  geom_point   (size = 2) +
  geom_abline  (intercept = 0,
                slope     = 1,
                color     ="black",size = 1) +
  annotation_logticks() +
  scale_y_continuous(limits = c(0,10))+
  scale_x_continuous(limits = c(0,10)) +
  ggtitle("(D) Broiler chicken calibration")+
  theme_bw() + labs (x = "Observed value (mg/L)", y = "Predicted value (mg/L)")+
  theme(plot.title = element_text(size=15, face="bold", family="A"), axis.text    
        = element_text(size = 15,  color = "black", face="bold", family="A"),
        axis.title   = element_text(size = 18, face = "bold", color = "black", family="A"),
        strip.text.x = element_text(size=1, face="bold", color="black", family="A"),
        legend.position = "none")+  #+scale_y_log10()+ 
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line())
p4

windows()
figure <- ggarrange( p1, p2, p3, p4)
figure


## Generate a chart with observed and modeled values
write.csv(Cost(Fit$par)$residuals,'regression plot for calibration.CSV', row.names = FALSE)
write.csv(Cost1$residuals,'regression plot for evaluation.CSV', row.names = FALSE)
























