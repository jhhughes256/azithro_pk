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
  PMLC = 0,  // polymorphonuclear leukocyte cytosol
  MUSCB = 0,  // muscle interstitial fluid (bound drug)
  SUBCB = 0,  // subcutis intersitial fluid (bound drug)
  PMLCB = 0,  // polymorphonuclear leukocyte cytosol (bound drug)
  LUNG = 0,  // lung interstitial fluid
  LUNGB = 0,  // lung interstitial fluid (bound drug)
  ALMA = 0,  // alveolar macrophage cytosol
  ALMAB = 0,  // alveolar macrophage cytosol (bound drug)

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
  POPDFPMLC = 52,  // distribution factor (PML cytosol), dimensionless
  POPDFLUNG = 140,  // distribution factor (lung), dimensionless
  POPDFALMA = 250,  // distribution factor (AM cytosol), dimensionless
  
// Literature/Derived Values
  FUP = 0.007581,  // fraction unionised in plasma, dimensionless
  FUPMLC = 0.0012,  // fraction unionised in PML cytosol, dimensionless
  RATPMLC = 0.351,  // ratio of PML cytosol to total PML conc., dimensionless
  RATPMLL = 13.33,  // ratio of PML lysosome to total PML conc., dimensionless
  
// Default covariate values
  WT = 79,
    
// Inter-individual Variability
// ETA values from input dataset are used by default
  ETA1 = 0,  // Random effect for CLF
  ETA2 = 0,  // Random effect for VCF
  ETA3 = 0,  // Random effect for KIN
  ETA4 = 0,  // Random effect for DFMUSC
  ETA5 = 0,  // Random effect for DFSUBC
  ETA6 = 0  // Random effect for DFPMLC

$OMEGA  
// Intra Individual Variability (full covariance-variance matrix)
  @annotated @block
  IIVCLF   : 0.097969 : intra-individual variability for CLF
  IIVVCF   : 0.28 1.2769 : intra-individual variability for VCF
  IIVKIN   : 0 0 0.00000484 : intra-individual variability for KIN
  IIVDFMUSC: 0 0 0 0.072361 : intra-individual variability for DFMUSC
  IIVDFSUBC: 0 0 0 0 0.099225 : intra-individual variability for DFSUBC
  IIVDFPMLC: 0 0 0 0 0 0.00000484 : intra-individual variability for DFPMLC

$SIGMA  
// Residual Unexplained Variability	(diagonal matrix)
  @annotated
  RUVPLASPRO: 1 : proportional RUV for plasma concentrations
  RUVMUSCPRO: 0.14 : proportional RUV for plasma concentrations
  RUVMUSCADD: 0.51 : proportional RUV for plasma concentrations
  RUVSUBCPRO: 0.34 : proportional RUV for plasma concentrations
  RUVSUBCADD: 0.000001 : proportional RUV for plasma concentrations
  RUVPMLCPRO: 0.23 : proportional RUV for plasma concentrations
  RUVPMLCADD: 0.000001 : proportional RUV for plasma concentrations

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
  double DFPMLC = POPDFPMLC*exp(ETA6);
  
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
    -KIN*FU*PLAS + KOUT*MUSC - KIN*FU*PLAS + KOUT*SUBC - KIN*FUP*FU*PLAS + KOUT*PMLC -
    KIN*FU*PLAS + KOUT*LUNG;
  dxdt_PER1 = K23*PLAS + -K32*PER1;
  dxdt_PER2 = K24*PLAS + -K42*PER2;
  
// Tissue Pharmacokinetic Differential Equations
  dxdt_MUSC = KIN*FU*PLAS - KOUT*MUSC + KOFF*MUSCB - KON*MUSC;
  dxdt_SUBC = KIN*FU*PLAS - KOUT*SUBC + KOFF*SUBCB - KON*SUBC;
  dxdt_PMLC = KIN*FUP*FU*PLAS - KOUT*PMLC + KOFF*PMLCB - KON*PMLC;
  dxdt_MUSCB = KON*MUSC - KOFF*MUSCB;
  dxdt_SUBCB = KON*SUBC - KOFF*SUBCB;
  dxdt_PMLCB = KON*PMLC - KOFF*PMLCB;

  dxdt_LUNG = KIN*FU*PLAS - KOUT*LUNG - KON*LUNG + KOFF*LUNGB -KIN*FUP*LUNG + KOUT*ALMA;
  dxdt_LUNGB = KON*LUNG - KOFF*LUNGB;
  dxdt_ALMA = KIN*FUP*LUNG - KOUT*ALMA - KON*ALMA + KOFF*ALMAB;
  dxdt_ALMAB = KON*ALMA - KOFF*ALMAB;
  
$TABLE  
// Individual Predictions
// Plasma Concentrations
  double CPLAST = SF*PLAS/VCF;  // plasma (total)
  double CPLAS = CPLAST*FU;  // plasma (unbound)
  
// Tissue Concentrations
  double CLUNG = SF*LUNG*KOUT*DFLUNG/(KIN*VCF);  // lung (unbound)
  double CMUSC = SF*MUSC*KOUT*DFMUSC/(KIN*VCF);  // muscle (unbound)
  double CSUBC = SF*SUBC*KOUT*DFSUBC/(KIN*VCF);  // subcutis (unbound)
  
// PML Concentrations
  double CPMLC = SF*PMLC*KOUT*DFPMLC/(KIN*VCF);  // PML cytosol (unionised)
  double CPMLCT = CPMLC/FUPMLC;  // PML cytosol (total)
  double CPMLT = CPMLCT/RATPMLC;  // PML (total)
  double CPMLLT = RATPMLL*CPMLT;  // PML lysosome (total)
 
// AM Concentrations
  double CALMAC = SF*ALMA*KOUT*DFALMA/(KIN*VCF);  // AM cytosol (unionised)
  double CALMACT = CALMAC/FUPMLC;  // AM cytosol (total)
  double CALMAT = CALMACT/RATPMLC;  // AM (total)
  double CALMALT = RATPMLL*CALMAT;  // AM lysosome (total)
  
// Dependent Variable
// Scaling factor must be applied to additive component of error model
// Additive EPS (ug/mL), IPRED (ng/mL)
  double DVPLAS = CPLAS*(1 + POPRUVPRO*EPS(1));
  double DVMUSC = CPLAS*(1 + EPS(2)) + SF*EPS(3);
  double DVSUBC = CPLAS*(1 + EPS(4)) + SF*EPS(5);
  double DVPMLC = CPLAS*(1 + EPS(6)) + SF*EPS(7);
  
$CAPTURE  
// Capture output
  DVPLAS DVMUSC DVSUBC DVPMLC // Dependent Variables
  CPLAS CPLAST CMUSC CSUBC CPMLC CPMLCT CPMLT CPMLLT // Individual Predictions
  CLUNG CALMAC CALMACT CALMAT CALMALT // Individual Predictions
  FU CLF VCF KIN DFMUSC DFSUBC DFPMLC  // Individual Parameters
  ETA1 ETA2 ETA3 ETA4 ETA5 ETA6  // Inter-individual Variability  
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile model 
  mod <- mrgsolve::mcode("DevMod_Azithromycin", code)

