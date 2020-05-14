# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Hybrid Model plus:
# - Lung compartment attached to plasma
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Set PATH for mrgsolve
# Get version of R
  Rversion <- paste("R-", getRversion(), sep="")
# Set path for Rtools
  Sys.setenv(PATH = paste0(glue::glue(
    "c:/program files/R/{Rversion}/bin/x64/;", "c:/RTools/bin/;",
    "c:/RTools/mingw_64/bin/;", Sys.getenv("PATH")
  )))  # Sys.setenv
    
# Define model code
  code <- '
$INIT
// Define compartment names
  DEPOT = 0,  // oral dose depot
  PLAS = 0,  // unbound plasma
  PER1 = 0,  // fast-distribution peripheral
  PER2 = 0,  // slow-distribution peripheral
  MUSC = 0,  // muscle interstitial fluid
  SUBC = 0,  // subcutis intersitial fluid
  WBCC = 0,  // white blood cell cytosol
  MUSCB = 0,  // muscle interstitial fluid (bound drug)
  SUBCB = 0,  // subcutis intersitial fluid (bound drug)
  WBCCB = 0,  // white blood cell cytosol (bound drug)
  LUNG = 0,  // lung interstitial fluid
  LUNGB = 0,  // lung interstitial fluid (bound drug)
  ALMA = 0,  // alveolar macrophage cytosol
  ALMAB = 0,  // alveolar macrophage cytosol (bound drug)
  PLASAUC = 0,  // plasma AUC compartment
  MUSCAUC = 0,  // muscle AUC compartment
  SUBCAUC = 0,  // subcutis AUC compartment
  LUNGAUC = 0,  // lung AUC compartment
  WBCCAUC = 0,  // white blood cell cytosol AUC compartment 
  ALMAAUC = 0,  // alveolar macrophage cytosol AUC compartment

$PARAM
// Fixed-Effects Parameters (Three Compartment Model)
  POPKA = 0.259,  // first-order absorption constant, hr^-1
  POPCLF = 100,  // apparent clearance, L/hr
  POPVCF = 186,  // apparent central volume, L
  POPVP1F = 2290,  // apparent peripheral volume (fast), L
  POPQP1F = 180,  // apparent intercompartment clearance (fast), L/hr
  POPVP2F = 2610,  // apparent peripheral volume (slow), L
  POPQP2F = 10.6,  // apparent intercompartment clearance (slow), L/hr
  POPRUVPRO = 0.406,  // proportional residual unexplained variability
  
// Fixed-Effects Parameters (Tissue/Cell Kinetic Model)
  POPKIN = 0.16,  // uptake into tissue (rate constant), hr^-1
  POPKOUT = 0.15,  // distribution out of tissue (rate constant), hr^-1
  POPKON = 0.56,  // non-specific tissue binding (on rate constant), hr^-1
  POPKOFF = 0.05,  // non-specific tissue binding (off rate constant), hr^-1
  POPDFMUSC = 0.55,  // distribution factor (muscle), dimensionless
  POPDFSUBC = 0.25,  // distribution factor (subcutis), dimensionless
  POPDFWBCC = 77,  // distribution factor (WBC cytosol), dimensionless
  POPDFLUNG = 140,  // distribution factor (lung), dimensionless
  POPDFALMA = 730,  // distribution factor (AM cytosol), dimensionless
  
// Literature/Derived Values
  FUP = 0.007581,  // fraction unionised in plasma, dimensionless
  FUWBCC = 0.0012,  // fraction unionised in WBC cytosol, dimensionless
  RATWBCC = 0.351,  // ratio of WBC cytosol to total WBC conc., dimensionless
  RATWBCL = 13.33,  // ratio of WBC lysosome to total WBC conc., dimensionless
  
// Default covariate values
  WT = 79,
    
// Inter-individual Variability
// ETA values from input dataset are used by default
  ETA1 = 0,  // Random effect for CLF
  ETA2 = 0,  // Random effect for VCF
  ETA3 = 0,  // Random effect for KIN
  ETA4 = 0,  // Random effect for DFMUSC
  ETA5 = 0,  // Random effect for DFSUBC
  ETA6 = 0  // Random effect for DFWBCC

$OMEGA  
// Intra Individual Variability (full covariance-variance matrix)
  @annotated @block
  IIVCLF   : 0.097969 : intra-individual variability for CLF
  IIVVCF   : 0.28 1.2769 : intra-individual variability for VCF
  IIVKIN   : 0 0 0.00000484 : intra-individual variability for KIN
  IIVDFMUSC: 0 0 0 0.072361 : intra-individual variability for DFMUSC
  IIVDFSUBC: 0 0 0 0 0.099225 : intra-individual variability for DFSUBC
  IIVDFWBCC: 0 0 0 0 0 0.00000484 : intra-individual variability for DFWBCC

$SIGMA  
// Residual Unexplained Variability	(diagonal matrix)
  @annotated
  RUVPLASPRO: 1 : proportional RUV for plasma concentrations
  RUVMUSCPRO: 0.14 : proportional RUV for plasma concentrations
  RUVMUSCADD: 0.51 : proportional RUV for plasma concentrations
  RUVSUBCPRO: 0.34 : proportional RUV for plasma concentrations
  RUVSUBCADD: 0.000001 : proportional RUV for plasma concentrations
  RUVWBCCPRO: 0.23 : proportional RUV for plasma concentrations
  RUVWBCCADD: 0.000001 : proportional RUV for plasma concentrations

$MAIN
// Pharmacokinetic Model Parameters
  double KA = POPKA;
  double CLF = POPCLF*pow(WT/70, 0.75)*exp(ETA1);
  double VCF = POPVCF*pow(WT/70, 1)*exp(ETA2);
  double VP1F = POPVP1F*pow(WT/70, 1);
  double QP1F = POPQP1F*pow(WT/70, 0.75);
  double VP2F = POPVP2F;
  double QP2F = POPQP2F;
  
// Micro-rate Constants
  double K20 = CLF/VCF;
  double K23 = QP1F/VCF;
  double K32 = QP1F/VP1F;
  double K24 = QP2F/VCF;
  double K42 = QP2F/VP2F;
  
// Tissue Kinetic Parameters
  double KIN = POPKIN*exp(ETA3);
  double KOUT = POPKOUT;
  double KON = POPKON;
  double KOFF = POPKOFF;
  
// Distribution Factors
  double DFMUSC = POPDFMUSC*exp(ETA4);
  double DFSUBC = POPDFSUBC*exp(ETA5);
  double DFWBCC = POPDFWBCC*exp(ETA6);
  
  double DFLUNG = POPDFLUNG;
  double DFALMA = POPDFALMA;
  
// Scaling Factor
  double SF = 1000;  // dose in mg, concentration in ng/mL

$ODE
// Plasma concentration and saturable fraction unbound
  double CPT = SF*PLAS/VCF;
  double FU = 0.4984+0.5339*CPT/(230.9+CPT);

// Three-Compartment Pharmacokinetic Differential Equations
  dxdt_DEPOT = -KA*DEPOT;
  dxdt_PLAS = KA*DEPOT - K23*PLAS + K32*PER1 - K24*PLAS + K42*PER2 - K20*PLAS +
    -KIN*FU*PLAS + KOUT*MUSC - KIN*FU*PLAS + KOUT*SUBC - KIN*FUP*FU*PLAS + KOUT*WBCC -
    KIN*FU*PLAS + KOUT*LUNG;
  dxdt_PER1 = K23*PLAS + -K32*PER1;
  dxdt_PER2 = K24*PLAS + -K42*PER2;
  
// Tissue Pharmacokinetic Differential Equations
  dxdt_MUSC = KIN*FU*PLAS - KOUT*MUSC + KOFF*MUSCB - KON*MUSC;
  dxdt_SUBC = KIN*FU*PLAS - KOUT*SUBC + KOFF*SUBCB - KON*SUBC;
  dxdt_WBCC = KIN*FUP*FU*PLAS - KOUT*WBCC + KOFF*WBCCB - KON*WBCC;
  dxdt_MUSCB = KON*MUSC - KOFF*MUSCB;
  dxdt_SUBCB = KON*SUBC - KOFF*SUBCB;
  dxdt_WBCCB = KON*WBCC - KOFF*WBCCB;

  dxdt_LUNG = KIN*FU*PLAS - KOUT*LUNG - KON*LUNG + KOFF*LUNGB -KIN*FUP*LUNG + KOUT*ALMA;
  dxdt_LUNGB = KON*LUNG - KOFF*LUNGB;
  dxdt_ALMA = KIN*FUP*LUNG - KOUT*ALMA - KON*ALMA + KOFF*ALMAB;
  dxdt_ALMAB = KON*ALMA - KOFF*ALMAB;

// Tissue Concentration
  double CMUS = SF*MUSC*KOUT*DFMUSC/(KIN*VCF);
  double CSUB = SF*SUBC*KOUT*DFSUBC/(KIN*VCF);
  double CLUN = SF*LUNG*KOUT*DFLUNG/(KIN*VCF);
  double CWBC = ((SF*WBCC*KOUT*DFWBCC/(KIN*VCF))/FUWBCC)/RATWBCC;
  double CAMC = ((SF*ALMA*KOUT*DFALMA/(KIN*VCF))/FUWBCC)/RATWBCC;
  
// Calculate Plasma and Tissue AUC
  dxdt_PLASAUC = CPT;
  dxdt_MUSCAUC = CMUS;
  dxdt_SUBCAUC = CSUB;
  dxdt_LUNGAUC = CLUN;
  dxdt_WBCCAUC = CWBC;
  dxdt_ALMAAUC = CAMC;

// Calculate Plasma and Tissue Cmax
if (SOLVERTIME == 0) {
  double PLASCMAX = 0;
  double MUSCCMAX = 0;
  double SUBCCMAX = 0;
  double LUNGCMAX = 0;
  double WBCCCMAX = 0;
  double ALMACMAX = 0;
}
if (CPT > PLASCMAX) PLASCMAX = CPT;
if (CMUS > MUSCCMAX) MUSCCMAX = CMUS;
if (CSUB > SUBCCMAX) SUBCCMAX = CSUB;
if (CLUN > LUNGCMAX) LUNGCMAX = CLUN;
if (CWBC > WBCCCMAX) WBCCCMAX = CWBC;
if (CAMC > ALMACMAX) ALMACMAX = CAMC;
  
$TABLE  
// Individual Predictions
// Plasma Concentrations
  double CPLAST = SF*PLAS/VCF;  // plasma (total)
  double CPLAS = CPLAST*FU;  // plasma (unbound)
  
// Tissue Concentrations
  double CLUNG = SF*LUNG*KOUT*DFLUNG/(KIN*VCF);  // lung (unbound)
  double CMUSC = SF*MUSC*KOUT*DFMUSC/(KIN*VCF);  // muscle (unbound)
  double CSUBC = SF*SUBC*KOUT*DFSUBC/(KIN*VCF);  // subcutis (unbound)
  
// WBC Concentrations
  double CWBCC = SF*WBCC*KOUT*DFWBCC/(KIN*VCF);  // WBC cytosol (unionised)
  double CWBCCT = CWBCC/FUWBCC;  // WBC cytosol (total)
  double CWBCT = CWBCCT/RATWBCC;  // WBC (total)
  double CWBCLT = RATWBCL*CWBCT;  // WBC lysosome (total)
 
// AM Concentrations
  double CALMAC = SF*ALMA*KOUT*DFALMA/(KIN*VCF);  // AM cytosol (unionised)
  double CALMACT = CALMAC/FUWBCC;  // AM cytosol (total)
  double CALMAT = CALMACT/RATWBCC;  // AM (total)
  double CALMALT = RATWBCL*CALMAT;  // AM lysosome (total)
  
// Dependent Variable
// Scaling factor must be applied to additive component of error model
// Additive EPS (ug/mL), IPRED (ng/mL)
//  double DVPLAS = CPLAS*(1 + POPRUVPRO*EPS(1));
//  double DVMUSC = CPLAS*(1 + EPS(2)) + SF*EPS(3);
//  double DVSUBC = CPLAS*(1 + EPS(4)) + SF*EPS(5);
//  double DVWBCC = CPLAS*(1 + EPS(6)) + SF*EPS(7);
  
$CAPTURE  
// Capture output
//  DVPLAS DVMUSC DVSUBC DVWBCC // Dependent Variables
  CPLAS CPLAST CMUSC CSUBC CWBCC CWBCCT CWBCT CWBCLT // Individual Predictions
  CLUNG CALMAC CALMACT CALMAT CALMALT // Individual Predictions
  FU CLF VCF KIN DFMUSC DFSUBC DFWBCC  // Individual Parameters
  ETA1 ETA2 ETA3 ETA4 ETA5 ETA6  // Inter-individual Variability  
  PLASAUC MUSCAUC SUBCAUC LUNGAUC WBCCAUC ALMAAUC  // Area Under the Curve
  PLASCMAX MUSCCMAX SUBCCMAX LUNGCMAX WBCCCMAX ALMACMAX  // Max. concentration
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile model 
  mod <- mrgsolve::mcode("DevMod_Azithromycin", code)
