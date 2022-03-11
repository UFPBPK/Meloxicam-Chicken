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


SolvePBPK <- '
$PARAM @annotated

 
// Oral absorption rate constants
  
Kst     : 2        :                  /h, gastric emptying rate constant
Kint    : 1.886    :                  /h, intestinal transit rate constant, fitted.

// IV infusion rate constants
Timeiv  : 0.01     :                  IV injection/infusion time (h)

// Urinary elimination rate constant adjusted by body weight
KurineC : 0.011    :                  L/h/kg, fitted.

// Metabolic rate constant adjusted by body weight
KmC     : 0.0234   :                 /h/kg, metabolic rate constant, fitted.

// Dosing, multiple oral gavage
Ka      : 0.568    :                 intestinal absorption rate constant (/h), fitted. 


  
  
// Physiological parameters
// Blood flow rates
QCC      : 9.91      :                 Cardiac output (L/h/kg), Table 13, Wang et al. 2021
QLCa     : 0.2526    :                 Fraction of flow to the liver, Table 16, Wang et al. 2021
QKCa     : 0.2012    :                 Fraction of flow to the kidneys, Table 16, Wang et al. 2021
QFCa     : 0.015     :                 Fraction of flow to the fat, Zeng et al. 2019
QMCa     : 0.35      :                 Fraction of flow to the muscle, Zeng et al. 2019
QovaryCa : 0.1646    :                 Fraction of flow to the ovary, Table 20, Wang et al. 2021
QRestCa  : 0.0166    :                 Calculated


// Tissue  volumes
BW       : 2.2       :                 Body weight(kg)
VLCa     : 0.0249    :                 Fractional liver tissue, Table 3, Wang et al. 2021
VKCa     : 0.0076    :                 Fractional kidney tissue, Table 3, Wang et al. 2021
VFCa     : 0.05      :                 Fractional fat tissue, Zeng et al. 2019
VMCa     : 0.4015    :                 Fractional muscle tissue, Table 2, Wang et al. 2021
VovaryCa : 0.0191    :                 Fractional ovary tissue, Table 3, Wang et al. 2021
VbloodCa : 0.063     :                 Blood volume, fraction of BW, Table 3, Wang et al. 2021
VRestCa  : 0.4354    :                 Calculated
Htc      : 0.314     :                 Hematocrit volume, fraction of BW, Table 22, Wang et al. 2021

  
// Mass Transfer Parameters (Chemical-specific parameters)
// Partition coefficients(PC, tissue:plasma)
PL      : 0.6116    :                 Liver plasma PC, calculated using AUCtissue:AUCplasma method from Dutch et al. 2021
PK      : 0.915     :                 Fitted
PM      : 0.0456    :                 Muscle plasma PC, breast muscle, calculated using AUCtissue:AUCplasma method from  Dutch et al. 2021 
Povary  : 0.0108    :                 Fitted.
PF      : 0.0174    :                 Fat plasma PC, calculated using AUCtissue:AUCplasma method from  Dutch et al. 2021
PRest   : 0.479     :                 Rest of the body tissues:plasma PC, Table 2, Li M et al. 2019
PB      : 0.994     :                 Percentage of drug bound to plasma proteins (DrugBank, 2020)


// Egg yolk and white parameters

Ay   : 25                             : apparent maximum follicle weight g
Ky   : 0.0682                         : transport constant into yolk, fitted.
Kw   : 0.0640                         : transport constant into albumen, fitted.
tlag : 1*24                           : h
tsig : 2*24                           : h
s    : 1/24                           : /h
talbumen : 10                         : h
tlay1  : 12                           : tlay1 default
Rwhitefor : (34*0.001)/10             : rate of albumen formation Kg/h 



$MAIN


double sumQ = QLCa + QKCa + QMCa + QFCa + QovaryCa + QRestCa; // sum up cardiac output fraction
 double QLC = QLCa/sumQ;                     // adjusted blood flow rate fraction to liver
 double QKC = QKCa/sumQ;                     // adjusted blood flow rate fraction to kidney
 double QMC = QMCa/sumQ;                     // adjusted blood flow rate fraction to muscle
 double QFC = QFCa/sumQ;                     // adjusted blood flow rate fraction to fat
 double QovaryC = QovaryCa/sumQ;             // adjusted blood flow rate fraction to ovary
 double QRestC = QRestCa/sumQ;               // adjusted blood flow rate fraction to rest of body
  
double sumV = VLCa + VKCa + VMCa + VFCa + VbloodCa + VovaryCa + VRestCa ; //sum up the tissue volumes
  
  double VLC = VLCa/sumV;                      //adjusted fraction of tissue volume of liver
  double VKC = VKCa/sumV;                      //adjusted fraction of tissue volume of kidney
  double VMC = VMCa/sumV;                      //adjusted fraction of tissue volume of muscle
  double VFC = VFCa/sumV;                      //adjusted fraction of tissue volume of fat
  double VovaryC = VovaryCa/sumV;              //adjusted fraction of tissue volume of ovary
  double VbloodC = VbloodCa/sumV;              //adjusted fraction of tissue volume of blood
  double VRestC = VRestCa/sumV;                //adjusted fraction of tissue volume of rest of body



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
    double Tmass = Aplas_free + AL + AK + Aurine + AM + AF + Aovary + ARest + Amet + Acolon + Aegg_ex;
    double bal   = ADOSE - Tmass-AI-AST;

$TABLE
capture Plasma = CV;
capture Liver  = CL;
capture Kidney = CK;
capture Muscle = CM;
capture Fat  = CF;
capture Ovary  = Covary;
capture Rest = CRest;
capture Yolk = Cyolk;
capture White = Cwhite;
capture Bal  = bal;
capture tmass = Tmass;
'





## Build mrgsolve-based PBPK Model

mod <- mcode_cache("pbpk", SolvePBPK) #refer to mcode function in mrgsolve user guide 3.1.2 Inline


## Physiological global parameter

BW.mean        = 2.2               # kg; accommodate to dataset


# IV infusion rate constants
Timeiv  = 0.01                     #                  IV injection/infusion time (h)


# Oral absorption rate constants

Kst.mean     = 2                   #                  /h, gastric emptying rate constant
Kint.mean    = 1.886               #                  /h, intestinal transit rate constant, fitted.

# Urinary elimination rate constant adjusted by body weight
KurineC.mean = 0.0108              #                  L/h/kg, fitted.

# Metabolic rate constant adjusted by body weight
KmC.mean     = 0.0234              #                 /h/kg, metabolic rate constant, fitted.

# Dosing, multiple oral gavage
Ka.mean      = 0.568               #                intestinal absorption rate constant (/h), fitted. 


## Cardiac Output and Blood flow
QCC.mean = 9.91                    # SD(1.81)# Cardiac output (L/h/kg), Table 13, Wang et al. 2021
QLC.mean = 0.2526                  # Fraction of flow to the liver, Table 16, Wang et al. 2021
QKC.mean = 0.2012                  # SD(0.1244) Fraction of flow to the kidneys, Table 16, Wang et al. 2021
QFC.mean = 0.015                   # Fraction of flow to the fat, Zeng et al. 2019
QMC.mean = 0.35                    #  Fraction of flow to the muscle, Zeng et al. 2019
QovaryC.mean = 0.1646              # SD(4.74/100) Fraction of flow to the ovary, Table 20, Wang et al. 2021
QrestC.mean=1-(QLC.mean+QKC.mean+QMC.mean+QFC.mean+QovaryC.mean)       # Blood flow in the rest of the body


####################### Tissue Volume ######################################################
############################################################################################

VLC.mean       = 0.0249             ## SD(0.49/100) Fractional liver tissue, Table 3, Wang et al. 2021
VKC.mean       = 0.0076             ## SD (0.13/100) Fractional kidney tissue, Table 3, Wang et al. 2021
VFC.mean       = 0.05               ## Fractional fat tissue, Zeng et al. 2019
VMC.mean       = 0.4015             ## SD (0.081/100) Fractional muscle tissue, Table 2, Wang et al. 2021
VovaryC.mean   = 0.0191             ## SD (0.28/100) Fractional ovary tissue, Table 3, Wang et al. 2021
VBloodC.mean   = 0.063              ## SD (1.33/100) Blood volume, fraction of BW , Table 3, Wang et al. 2021        
Htc.mean =0.314                     ## SD (5.1/100) # Hematocrit volume, fraction of BW, Table 22, Wang et al. 2021
VplasmaC = (1-Htc.mean) * VBloodC.mean # Plasma volume, fraction of BW, Table 22, Wang et al. 2021
VRestC.mean=1-(VLC.mean+VKC.mean+VMC.mean+VFC.mean+VovaryC.mean+VBloodC.mean)

############# Mass Transfer Parameters (Chemical-specific parameters) ########################
############### Partition coefficients(PC, tissue:plasma) ####################################

PL.mean = 0.6116                    ## Liver plasma PC, calculated using AUCtissue:AUCplasma method from Dutch et al. 2021
PK.mean = 0.915                     ## Fitted
PM.mean = 0.0456                    ## Muscle plasma PC, breast muscle, calculated using AUCtissue:AUCplasma method from  Dutch et al. 2021 
Povary.mean = 0.0108                ## Fitted
PF.mean = 0.0174                    ## Fat plasma PC, calculated using AUCtissue:AUCplasma method from  Dutch et al. 2021
Prest.mean = 0.479                  ## Rest of the body tissues:plasma PC, Table 2, Li M et al. 2019
PB.mean = 0.994                     ## Percentage of drug bound to plasma proteins (DrugBank, 2020)
Free.mean = 1-PB.mean               ## Unitless; fraction of chemical free in blood after plasma protein binding


######################## Standard deviations for each parameters ######################################
#######################################################################################################

PC.CV=0.2                      ## Coefficient of variation for partition coefficients sd/mean Li 2017
K.CV=0.3                       ## Coefficient of variation for absorption/excretion parameters sd/mean Li 2017
P.CV=0.3                       ## Coefficient of variation for physiological parameters sd/mean Li 2017

###############################################################################################################
#####
QCC.sd = 1.81                  ## Standard deviation for cardiac output Table 13, Wang et al. 2021             
QLC.sd = P.CV*QLC.mean         ## Standard deviation for fractions of blood flow to liver 
QKC.sd = P.CV*QKC.mean         ## (negative number)Standard deviation for fractions of blood flow to Kidney Table 16, Wang et al. 2021
QMC.sd = P.CV*QMC.mean         ## Standard deviation for fractions of blood flow to Muscle 
QFC.sd = P.CV*QFC.mean         ## Standard deviation for fractions of blood flow to fat
QovaryC.sd = 4.74/100          ## Standard deviation for fraction of flow to the ovary, Table 20, Wang et al. 2021
QrestC.sd = P.CV*QrestC.mean   ## Standard deviation for fractions of blood flow to rest of the body


BW.sd = P.CV*BW.mean           ## Standard deviation for bodyweight 
VLC.sd = 0.49/100              ## Standard deviation for fractions of volume of liver Table 3, Wang et al. 2021
VKC.sd = 0.13/100              ## Standard deviation for fractions of volume of kidney Table 2, Wang et al. 2021
VMC.sd = 0.081/100             ## Standard deviation for fractions of volume of muscle Table 2, Wang et al. 2021
VFC.sd = P.CV*VFC.mean         ## Standard deviation for fractions of volume of fat
VovaryC.sd   = 0.28/100        ## Standard deviation for fractional ovary tissue, Table 3, Wang et al. 2021
VRestC.sd = P.CV*VRestC.mean   ## Standard deviation for fractions of volume of rest of the body



PL.sd = PC.CV*PL.mean          ## Standard deviation for liver:Plasma Partition coefficient
PK.sd = PC.CV*PK.mean          ## Standard deviation for kidney:Plasma Partition coefficient
PM.sd =PC.CV*PM.mean           ## Standard deviation for muscle:Plasma Partition coefficient
PF.sd = PC.CV*PF.mean          ## Standard deviation for fat:Plasma Partition coefficient
Povary.sd = PC.CV*Povary.mean  ## Standard deviation for ovary:Plasma Partition coefficient
Prest.sd = PC.CV*Prest.mean    ## Standard deviation for rest of the body:Plasma Partition coefficient


Kst.sd = K.CV*Kst.mean           ## Standard deviation for gastric emptying rate
Kint.sd = K.CV*Kint.mean         ## Standard deviation for intestinal transit rate
KurineC.sd = K.CV*KurineC.mean   ## Standard deviation for urinary excretion rate
KmC.sd = K.CV*KmC.mean           ## Standard deviation for metabolic rate
Ka.sd = K.CV*Ka.mean             ## Standard deviation for intestinal absorption rate
PB.sd = K.CV*PB.mean             ##Standard deviation for percentage of drug bound to plasma proteins 




m.log.PL= log(PL.mean^2/(PL.sd^2+PL.mean^2)^0.5) 

sd.log.PL = (log(1+PL.sd^2/PL.mean^2))^0.5 

m.log.PK = log(PK.mean^2/(PK.sd^2+PK.mean^2)^0.5) 

sd.log.PK = (log(1+PK.sd^2/PK.mean^2))^0.5 

m.log.PM = log(PM.mean^2/(PM.sd^2+PM.mean^2)^0.5)

sd.log.PM = (log(1+PM.sd^2/PM.mean^2))^0.5

m.log.PF= log(PF.mean^2/(PF.sd^2+PF.mean^2)^0.5)

sd.log.PF = (log(1+PF.sd^2/PF.mean^2))^0.5

m.log.Povary = log(Povary.mean^2/(Povary.sd^2+PM.mean^2)^0.5)

sd.log.Povary = (log(1+Povary.sd^2/Povary.mean^2))^0.5

m.log.Prest= log(Prest.mean^2/(Prest.sd^2+Prest.mean^2)^0.5)

sd.log.Prest = (log(1+Prest.sd^2/Prest.mean^2))^0.5





m.log.Kst= log(Kst.mean^2/(Kst.sd^2+Kst.mean^2)^0.5)

sd.log.Kst = (log(1+Kst.sd^2/Kst.mean^2))^0.5

m.log.Kint= log(Kint.mean^2/(Kint.sd^2+Kint.mean^2)^0.5)

sd.log.Kint = (log(1+Kint.sd^2/Kint.mean^2))^0.5

m.log.KurineC= log(KurineC.mean^2/(KurineC.sd^2+KurineC.mean^2)^0.5) 

sd.log.KurineC = (log(1+KurineC.sd^2/KurineC.mean^2))^0.5

m.log.KmC= log(KmC.mean^2/(KmC.sd^2+KmC.mean^2)^0.5)

sd.log.KmC = (log(1+KmC.sd^2/KmC.mean^2))^0.5

m.log.Ka= log(Ka.mean^2/(Ka.sd^2+Ka.mean^2)^0.5) 

sd.log.Ka = (log(1+Ka.sd^2/Ka.mean^2))^0.5 

m.log.PB= log(PB.mean^2/(PB.sd^2+PB.mean^2)^0.5) 

sd.log.PB = (log(1+PB.sd^2/PB.mean^2))^0.5 



# set random seed so that the simulation result is reproducible, because randomly generated data is same if you set same random seed.
N=1000
set.seed(324) 
idata <- 
  data_frame(ID=1:N) %>% 
  mutate(
    QCC = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QCC.mean, sd = QCC.sd),
      b = qnorm(0.975, mean = QCC.mean, sd = QCC.sd),
      mean = QCC.mean,
      sd = QCC.sd
    ),  # Cardiac output index (L/h/kg)
    # Fracion of blood flow to organs (unitless) 
    QLCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QLC.mean, sd = QLC.sd),
      b = qnorm(0.975, mean = QLC.mean, sd = QLC.sd),
      mean = QLC.mean,
      sd = QLC.sd
    ),
    # Fraction of blood flow to the kidneys (2016 Lin)
    QKCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QKC.mean, sd =  QKC.sd),
      b = qnorm(0.975, mean = QKC.mean, sd =  QKC.sd),
      mean = QKC.mean,
      sd =  QKC.sd
    ),
    #QMC = 0.180			# Fraction of blood flow to the muscle (2016 Lin)
    QMCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QMC.mean, sd =  QMC.sd),
      b = qnorm(0.975, mean = QMC.mean, sd =  QMC.sd),
      mean = QMC.mean,
      sd =  QMC.sd
    ),
    # Fraction of blood flow to the fat (2016 Lin)
    QFCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QFC.mean, sd =  QFC.sd),
      b = qnorm(0.975, mean = QFC.mean, sd =  QFC.sd),
      mean = QFC.mean,
      sd =  QFC.sd
    ),
    # Fraction of blood flow to the ovary 
    QovaryCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QovaryC.mean, sd =  QovaryC.sd),
      b = qnorm(0.975, mean = QovaryC.mean, sd =  QovaryC.sd),
      mean = QovaryC.mean,
      sd =  QovaryC.sd
    ),
    
    # Fraction of blood flow to the rest of body (total sum equals to 1)
    QRestCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = QrestC.mean, sd = QrestC.sd),
      b = qnorm(0.975, mean = QrestC.mean, sd = QrestC.sd),
      mean = QrestC.mean,
      sd = QrestC.sd
    ),
    BW = rtruncnorm(n = N, 
                    a = qnorm(0.025, mean = BW.mean, sd = BW.sd), 
                    b = qnorm(0.975, mean = BW.mean, sd = BW.sd), 
                    mean = BW.mean, sd = BW.sd), 
    # Fractional organ tissue volumes (unitless)
    # Fractional liver tissue 
    VLCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VLC.mean, sd =  VLC.sd),
      b = qnorm(0.975, mean = VLC.mean, sd =  VLC.sd) ,
      mean = VLC.mean,
      sd =  VLC.sd
    ),
    # Fractional kidney tissue 
    VKCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VKC.mean, sd =  VKC.sd),
      b = qnorm(0.975, mean = VKC.mean, sd =  VKC.sd),
      mean = VKC.mean,
      sd =  VKC.sd
    ),
    # Fractional fat tissue 
    VFCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VFC.mean, sd =  VFC.sd),
      b = qnorm(0.975, mean = VFC.mean, sd =  VFC.sd),
      mean = VFC.mean,
      sd =  VFC.sd
    ),
    # Fractional muscle tissue 
    VMCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VMC.mean, sd =  VMC.sd),
      b = qnorm(0.975, mean = VMC.mean, sd =  VMC.sd),
      mean = VMC.mean,
      sd =  VMC.sd
    ),
    # Fractional ovary tissue 
    VovaryCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VovaryC.mean, sd =  VovaryC.sd),
      b = qnorm(0.975, mean = VovaryC.mean, sd =  VovaryC.sd),
      mean = VovaryC.mean,
      sd =  VovaryC.sd
    ),
    # Fractional rest of body (total sum equals to 1)
    VRestCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = VRestC.mean, sd =  VRestC.sd),
      b = qnorm(0.975, mean = VRestC.mean, sd =  VRestC.sd),
      mean = VRestC.mean,
      sd =  VRestC.sd
    ),
    PL = rlnormTrunc(
      N,
      meanlog = m.log.PL,
      sdlog = sd.log.PL,
      min = qlnorm(0.025, meanlog = m.log.PL, sdlog = sd.log.PL),
      max = qlnorm(0.975, meanlog = m.log.PL, sdlog = sd.log.PL)
    ),
    PK = rlnormTrunc(
      N,
      meanlog = m.log.PK,
      sdlog = sd.log.PK,
      min = qlnorm(0.025, meanlog = m.log.PK, sdlog = sd.log.PK),
      max = qlnorm(0.975, meanlog = m.log.PK, sdlog = sd.log.PK)
    ),
    PM = rlnormTrunc(
      N,
      meanlog = m.log.PM,
      sdlog = sd.log.PM,
      min = qlnorm(0.025, meanlog = m.log.PM, sdlog = sd.log.PM),
      max = qlnorm(0.975, meanlog = m.log.PM, sdlog = sd.log.PM)
    ),
    PF = rlnormTrunc(
      N,
      meanlog = m.log.PF,
      sdlog = sd.log.PF,
      min = qlnorm(0.025, meanlog = m.log.PF, sdlog = sd.log.PF),
      max = qlnorm(0.975, meanlog = m.log.PF, sdlog = sd.log.PF)
    ),
    Povary = rlnormTrunc(
      N,
      meanlog = m.log.PF,
      sdlog = sd.log.PF,
      min = qlnorm(0.025, meanlog = m.log.PF, sdlog = sd.log.PF),
      max = qlnorm(0.975, meanlog = m.log.PF, sdlog = sd.log.PF)
    ),
    
    PRest = rlnormTrunc(
      N,
      meanlog = m.log.Prest,
      sdlog = sd.log.Prest,
      min = qlnorm(0.025, meanlog = m.log.Prest, sdlog = sd.log.Prest),
      max = qlnorm(0.975, meanlog = m.log.Prest, sdlog = sd.log.Prest)
    ),
    #{Kinetic Constants}
    ## gastric emptying rate
    Kst = rlnormTrunc(
      N,
      meanlog = m.log.Kst,
      sdlog = sd.log.Kst,
      min = qlnorm(0.025, meanlog = m.log.Kst, sdlog = sd.log.Kst),
      max = qlnorm(0.975, meanlog = m.log.Kst, sdlog = sd.log.Kst)
    ),
    ## intestinal transit rate
    Kint = rlnormTrunc(
      N,
      meanlog = m.log.Kint,
      sdlog = sd.log.Kint,
      min = qlnorm(0.025, meanlog = m.log.Kint, sdlog = sd.log.Kint),
      max = qlnorm(0.975, meanlog = m.log.Kint, sdlog = sd.log.Kint)
    ),
    ## intestinal absorption rate
    Ka = rlnormTrunc(
      N,
      meanlog = m.log.Ka,
      sdlog = sd.log.Ka,
      min = qlnorm(0.025, meanlog = m.log.Ka, sdlog = sd.log.Ka),
      max = qlnorm(0.975, meanlog = m.log.Ka, sdlog = sd.log.Ka)
    ),
    # Percentage Plasma Protein Binding unitless
    #   PB = rlnormTrunc(
    # N,
    #   meanlog = m.log.PB,
    #   sdlog = sd.log.PB,
    #    min = qlnorm(0.025, meanlog = m.log.PB, sdlog = sd.log.PB),
    #     max = 1
    #    ),
    #{Metabolic Rate Constant}
    KmC = rlnormTrunc(
      N,
      meanlog = m.log.KmC,
      sdlog = sd.log.KmC,
      min = qlnorm(0.025, meanlog = m.log.KmC, sdlog = sd.log.KmC),
      max = qlnorm(0.975, meanlog = m.log.KmC, sdlog = sd.log.KmC)
    ),
    # Urinary Elimination Rate Constants
    KurineC = rlnormTrunc(
      N,
      meanlog = m.log.KurineC,
      sdlog = sd.log.KurineC,
      min = qlnorm(0.025, meanlog = m.log.KurineC, sdlog = sd.log.KurineC),
      max = qlnorm(0.975, meanlog = m.log.KurineC, sdlog = sd.log.KurineC)
    )
  )


## Combine data and run the simulation


pred<- function (Dose, tlay1=12){
  
  tinterval = 12
  TDOSE     = 20
  BW        = 2.2
  DOSE      = Dose*BW
  route     = 'oral'
  

  if (route == "iv") {
    
    ev_1 <- ev (amt  = DOSE, ii = tinterval, tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "Aplas_free", replicate = FALSE)
    
    ev_2 <- ev (amt  = DOSE, ii = tinterval, tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
    
    ex <- ev_1+ev_2
    
    
  }
  
  
  
  if (route == "oral") {
    
    ev_1 <- ev (amt  = DOSE, ii = tinterval,
                
                addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
    
    ev_2 <- ev (amt  = DOSE, ii = tinterval,
                
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    
    ex <- ev_1+ev_2
    
  }
  
  
  
  tsamp  = tgrid(0, tinterval*(TDOSE-1)+85*24, 0.1)
  
  
  
  ## Simulation
  
  out <- mod %>% param(c(tlay1=tlay1))%>%
         data_set(ex) %>%
         update(atol= 1E-13, rtol=1E-6, maxsteps=50000) %>%
         idata_set(idata1) %>%
         mrgsim(obsonly=TRUE, tgrid=tsamp)
  
  
  output <- cbind.data.frame(Time = out$time,
                             Plasma  = out$Plasma,
                             Liver  = out$Liver,
                             Kidney  = out$Kidney,
                             Muscle  = out$Muscle,
                             Fat  = out$Fat,
                             Ovary = out$Ovary,
                             Yolk = out$Yolk, 
                             White = out$White, 
                             Bal = out$Bal)
  
  return(output)
  
}

output1 <- list()
tout2 <- list()
tout3 <- list()

t1 <- c((1:25)*24)
#t1 <- c(229, 232,234, 236, 240, 244, 252, 264, 276)






for (i in 1:N) {

idata1 <- idata[i,]%>%select(-ID)
#ex1 <- ex %>% filter(ID==i)

output1[[i]] <- pred(Dose =1)
 

for (j in 1:length(t1)) {
  output2 <- pred(Dose =1, tlay1 = t1[j]/24)
  tout2[[j]] <-output2 %>% filter(Time == t1[j])
  }
tout3[[i]] <- do.call(rbind, tout2)

}

out1 <- do.call(rbind, output1)  #plasma monte-carlo

tout4<- do.call(rbind, tout3)   # egg yolk and white montecarlo



tinterval = 12
TDOSE     = 20
out.sum <- out1 %>% as_data_frame %>%
  mutate(Time1 = Time/24 - (TDOSE-1)*tinterval/24) %>%
  dplyr::select(Time1, Plasma:White) %>%
  arrange(Time1)  %>%
  group_by(Time1) %>%
  dplyr::summarise(CV1 = quantile(Plasma, 0.01), CV50 = median(Plasma), CV99 = quantile(Plasma, 0.99),
                   CL1 = quantile(Liver, 0.01), CL50 = median(Liver), CL99 = quantile(Liver, 0.99),
                   CK1 = quantile(Kidney, 0.01), CK50 = median(Kidney), CK99 = quantile(Kidney, 0.99),
                   CM1 = quantile(Muscle, 0.01), CM50 = median(Muscle), CM99 = quantile(Muscle, 0.99),
                   CF1 = quantile(Fat, 0.01), CF50 = median(Fat), CF99 = quantile(Fat, 0.99),
                   Covary1 = quantile(Ovary, 0.01), Covary50 = median(Ovary), Covary99 = quantile(Ovary, 0.99))



output_plot=as.data.frame(out.sum)





out.sum1 <- tout4 %>% as_data_frame %>%
  mutate(Time1 = Time/24 - (TDOSE-1)*tinterval/24) %>%
  dplyr::select(Time1, Yolk:White) %>%
  arrange(Time1)  %>%
  group_by(Time1) %>%
  dplyr::summarise(Cyolk1 = quantile(Yolk, 0.01), Cyolk50 = median(Yolk), Cyolk99 = quantile(Yolk, 0.99),
                   Cwhite1 = quantile(White, 0.01), Cwhite50 = median(White), Cwhite99 = quantile(White, 0.99))

output_plot1=as.data.frame(out.sum1)


#########################################################################################################
#######draw plots for sfigures


tinterval = 12
TDOSE     = 20
out.sum <- out1 %>% as_data_frame %>%
  mutate(Time1 = Time/24 - (TDOSE-1)*tinterval/24) %>%
  dplyr::select(Time1, Plasma:White) %>%
  arrange(Time1)  %>%
  group_by(Time1) %>%
  dplyr::summarise(CV1 = quantile(Plasma, 0.01), CV50 = median(Plasma), CV99 = quantile(Plasma, 0.99),
                   CL1 = quantile(Liver, 0.01), CL50 = median(Liver), CL99 = quantile(Liver, 0.99),
                   CK1 = quantile(Kidney, 0.01), CK50 = median(Kidney), CK99 = quantile(Kidney, 0.99),
                   CM1 = quantile(Muscle, 0.01), CM50 = median(Muscle), CM99 = quantile(Muscle, 0.99),
                   CF1 = quantile(Fat, 0.01), CF50 = median(Fat), CF99 = quantile(Fat, 0.99),
                   Covary1 = quantile(Ovary, 0.01), Covary50 = median(Ovary), Covary99 = quantile(Ovary, 0.99))



output_plot=as.data.frame(out.sum)





out.sum1 <- tout4 %>% as_data_frame %>%
  mutate(Time1 = Time) %>%
  dplyr::select(Time1, Yolk:White) %>%
  arrange(Time1)  %>%
  group_by(Time1) %>%
  dplyr::summarise(Cyolk1 = quantile(Yolk, 0.01), Cyolk50 = median(Yolk), Cyolk99 = quantile(Yolk, 0.99),
                   Cwhite1 = quantile(White, 0.01), Cwhite50 = median(White), Cwhite99 = quantile(White, 0.99))

output_plot1=as.data.frame(out.sum1)

write.csv(out.sum1,'Monte Carlo Laying Hen1Sfigure 20 12.CSV', row.names = FALSE)

####################################################################################################################





















## Draw ggplot for Monte Carlo concentration
A = ggplot(output_plot, aes(Time1)) + 
  geom_ribbon(aes(ymin = CV1,
                  ymax = CV99
  ), fill = 'grey',
  show.legend = T, 
  size = 0.2,
  alpha = 0.5) +  # alpha is transparency parameter. 
  geom_line(aes(y = CV50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = CV99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = CV1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)')) + 
  ylab(expression(paste('Concentration (mg/kg)'))) + 
  #geom_line(aes(y = 2),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  theme_bw() + theme(axis.text=element_text(size=16), legend.title = element_blank())+
  #geom_vline(aes(xintercept = min(subset(output_plot, CV99<= 2 & Time1 > 0) %>% dplyr::select(Time1))), 
  #          size = 1, color = "red", linetype = 2,
  #          show.legend = F)+ 
  ggtitle("(A) Plasma, 1 mg/kg Oral dose") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +theme_bw() +theme(plot.title = element_text(size=16, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 14,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 14, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=14, face="bold", color="black", family="A"),
                                                                        strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                        legend.position = "none") +scale_y_log10()







B = ggplot(output_plot, aes(Time1)) + 
  geom_ribbon(aes(ymin = CL1,
                  ymax = CL99
  ), fill = 'grey',
  show.legend = T, 
  size = 0.2,
  alpha = 0.5) +  # alpha is transparency parameter. 
  geom_line(aes(y = CL50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = CL99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = CL1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)')) + 
  ylab(expression(paste('Concentration (mg/kg)'))) + 
  #geom_line(aes(y = 2),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  theme_bw() + theme(axis.text=element_text(size=16), legend.title = element_blank())+
  #geom_vline(aes(xintercept = min(subset(output_plot, CV99<= 2 & Time1 > 0) %>% dplyr::select(Time1))), 
  #          size = 1, color = "red", linetype = 2,
  #          show.legend = F)+ 
  ggtitle("(B) Liver, 1 mg/kg Oral dose") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +theme_bw() +theme(plot.title = element_text(size=16, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 14,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 14, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=14, face="bold", color="black", family="A"),
                                                                        strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                        legend.position = "none") +scale_y_log10()



C = ggplot(output_plot, aes(Time1)) + 
  geom_ribbon(aes(ymin = CK1,
                  ymax = CK99
  ), fill = 'grey',
  show.legend = T, 
  size = 0.2,
  alpha = 0.5) +  # alpha is transparency parameter. 
  geom_line(aes(y = CK50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = CK99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = CK1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)')) + 
  ylab(expression(paste('Concentration (mg/kg)'))) + 
  #geom_line(aes(y = 2),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  theme_bw() + theme(axis.text=element_text(size=16), legend.title = element_blank())+
  #geom_vline(aes(xintercept = min(subset(output_plot, CV99<= 2 & Time1 > 0) %>% dplyr::select(Time1))), 
  #          size = 1, color = "red", linetype = 2,
  #          show.legend = F)+ 
  ggtitle("(C) Kidney, 1 mg/kg Oral dose") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +theme_bw() +theme(plot.title = element_text(size=16, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 14,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 14, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=14, face="bold", color="black", family="A"),
                                                                        strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                        legend.position = "none") +scale_y_log10()



D = ggplot(output_plot, aes(Time1)) + 
  geom_ribbon(aes(ymin = CM1,
                  ymax = CM99
  ), fill = 'grey',
  show.legend = T, 
  size = 0.2,
  alpha = 0.5) +  # alpha is transparency parameter. 
  geom_line(aes(y = CM50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = CM99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = CM1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)')) + 
  ylab(expression(paste('Concentration (mg/kg)'))) + 
  #geom_line(aes(y = 2),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  theme_bw() + theme(axis.text=element_text(size=16), legend.title = element_blank())+
  #geom_vline(aes(xintercept = min(subset(output_plot, CV99<= 2 & Time1 > 0) %>% dplyr::select(Time1))), 
  #          size = 1, color = "red", linetype = 2,
  #          show.legend = F)+ 
  ggtitle("(D) Muscle, 1 mg/kg Oral dose") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +theme_bw() +theme(plot.title = element_text(size=16, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 14,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 14, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=14, face="bold", color="black", family="A"),
                                                                        strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                        legend.position = "none") +scale_y_log10()




E = ggplot(output_plot, aes(Time1)) + 
  geom_ribbon(aes(ymin = CF1,
                  ymax = CF99
  ), fill = 'grey',
  show.legend = T, 
  size = 0.2,
  alpha = 0.5) +  # alpha is transparency parameter. 
  geom_line(aes(y = CF50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = CF99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = CF1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)')) + 
  ylab(expression(paste('Concentration (mg/kg)'))) + 
  #geom_line(aes(y = 2),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  theme_bw() + theme(axis.text=element_text(size=16), legend.title = element_blank())+
  #geom_vline(aes(xintercept = min(subset(output_plot, CV99<= 2 & Time1 > 0) %>% dplyr::select(Time1))), 
  #          size = 1, color = "red", linetype = 2,
  #          show.legend = F)+ 
  ggtitle("(E) Fat, 1 mg/kg Oral dose") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +theme_bw() +theme(plot.title = element_text(size=16, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 14,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 14, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=14, face="bold", color="black", family="A"),
                                                                        strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                        legend.position = "none") +scale_y_log10()







FF = ggplot(output_plot, aes(Time1)) + 
  geom_ribbon(aes(ymin = Covary1,
                  ymax = Covary99
  ), fill = 'grey',
  show.legend = T, 
  size = 0.2,
  alpha = 0.5) +  # alpha is transparency parameter. 
  geom_line(aes(y = Covary50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = Covary99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = Covary1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)')) + 
  ylab(expression(paste('Concentration (mg/kg)'))) + 
  #geom_line(aes(y = 2),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  theme_bw() + theme(axis.text=element_text(size=16), legend.title = element_blank())+
  #geom_vline(aes(xintercept = min(subset(output_plot, CV99<= 2 & Time1 > 0) %>% dplyr::select(Time1))), 
  #          size = 1, color = "red", linetype = 2,
  #          show.legend = F)+ 
  ggtitle("(F) Ovary, 1 mg/kg Oral dose") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +theme_bw() +theme(plot.title = element_text(size=16, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 14,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 14, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=14, face="bold", color="black", family="A"),
                                                                        strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                        legend.position = "none") +scale_y_log10()










G = ggplot(output_plot1, aes(Time1)) + 
  geom_ribbon(aes(ymin = Cyolk1,
                  ymax = Cyolk99
  ), fill = 'grey',
  show.legend = T, 
  size = 0.2,
  alpha = 0.5) +  # alpha is transparency parameter. 
  geom_line(aes(y = Cyolk50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = Cyolk99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = Cyolk1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)')) + 
  ylab(expression(paste('Concentration (mg/kg)'))) + 
  #geom_line(aes(y = 2),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  theme_bw() + theme(axis.text=element_text(size=16), legend.title = element_blank())+
  #geom_vline(aes(xintercept = min(subset(output_plot, CV99<= 2 & Time1 > 0) %>% dplyr::select(Time1))), 
  #          size = 1, color = "red", linetype = 2,
  #          show.legend = F)+ 
  ggtitle("(G) Egg yolk, 1 mg/kg Oral dose") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +theme_bw() +theme(plot.title = element_text(size=16, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 14,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 14, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=14, face="bold", color="black", family="A"),
                                                                        strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                        legend.position = "none") +scale_y_log10()








H = ggplot(output_plot1, aes(Time1)) + 
  geom_ribbon(aes(ymin = Cwhite1,
                  ymax = Cwhite99
  ), fill = 'grey',
  show.legend = T, 
  size = 0.2,
  alpha = 0.5) +  # alpha is transparency parameter. 
  geom_line(aes(y = Cwhite50,
                color = '50 Percentile'), 
            size = 1, 
            show.legend = T) + # draw the mean value to the chart.
  geom_line(aes(y = Cwhite99,
                color = '99 Percentile'), 
            size = 1, 
            show.legend = T) +  
  geom_line(aes(y = Cwhite1,
                color = '1 Percentile'), 
            size = 1, 
            show.legend = T) + 
  scale_x_continuous(name = c('Time (Day)')) + 
  ylab(expression(paste('Concentration (mg/kg)'))) + 
  #geom_line(aes(y = 2),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F) + # add tolerance line
  theme_bw() + theme(axis.text=element_text(size=16), legend.title = element_blank())+
  #geom_vline(aes(xintercept = min(subset(output_plot, CV99<= 2 & Time1 > 0) %>% dplyr::select(Time1))), 
  #          size = 1, color = "red", linetype = 2,
  #          show.legend = F)+ 
  ggtitle("(H) Egg white, 1 mg/kg Oral dose") +
  
  ylab("Concentration (mg/kg)")+ xlab("Time (days)") +theme_bw() +theme(plot.title = element_text(size=16, face="bold", family="A"), axis.text    
                                                                        = element_text(size = 14,  color = "black", face="bold", family="A"),
                                                                        axis.title   = element_text(size = 14, face = "bold", color = "black", family="A"),
                                                                        strip.text.x = element_text(size=14, face="bold", color="black", family="A"),
                                                                        strip.text.y = element_text(size=11, face="bold", color="black", family="A"),
                                                                        legend.position = "none") +scale_y_log10()



## Combine all of the ggplots in one figure
windows()
figure <- ggarrange( A,B,C,D,E,FF,G,H)
figure



## Create a chart for median, 1 percentile, 99 percentile concentration for plasma, tissue and egg via given dose time and intervals.
write.csv(out.sum,'Monte Carlo Laying Hen 20 12.CSV', row.names = FALSE)
write.csv(out.sum1,'Monte Carlo Laying Hen1 20 12.CSV', row.names = FALSE)



































































