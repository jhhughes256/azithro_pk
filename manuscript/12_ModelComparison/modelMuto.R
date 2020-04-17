# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Model sourced from:
#   Muto, C., et al. (2010). "Pharmacokinetic-pharmacodynamic analysis of 
#   azithromycin extended release in Japanese patients with common respiratory 
#   tract infectious disease." Journal of Antimicrobial 
#   Chemotherapy 66(1): 165-174.
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
  PERI = 0,  // fast-distribution peripheral

$PARAM
// Fixed-Effects Parameters
  POPKA = 0.725,  // first-order absorption constant, hr^-1
  POPALAG = 0.457,  // absorption lag, hr
  POPCLF = 103,  // apparent clearance, L/hr
  POPVCF = 1830,  // apparent central volume, L
  POPVPF = 4340,  // apparent peripheral volume (fast), L
  POPQPF = 138,  // apparent intercompartment clearance (fast), L/hr
  
// Covariate-Effect Parameters
  COVWTCLF = 0.971,  // Effect of body weight on CLF
  COVWTVCF = 1.03,  // Effect of body weight on VCF
  COVWTVPF = 1,  // Effect of body weight on VPF
  COVWTQPF = 0.75,  // Effect of body weight on QPF
  COVAGECLF = -0.166,  // Effect of age on CLF
  COVAGEVCF = -0.256,  // Effect of age on VCF
  
// Default covariate values
  WT = 60,  // median body weight from Muto et al.
  AGE = 36,  // median age from Muto et al.
    
// Inter-individual Variability
// ETA values from input dataset are used by default
  ETA1 = 0,  // Random effect for KA
  ETA2 = 0,  // Random effect for CLF
  ETA3 = 0,  // Random effect for VCF

$OMEGA  
// Intra Individual Variability (diagonal matrix)
  @annotated @block
  IIVKA  : 1.17 : intra-individual variability for KA
  IIVCLF : 0.102 0.115 : intra-individual variability for CLF
  IIVVCF : 0.246 0.0943 0.185 : intra-individual variability for VCF

$SIGMA  
// Residual Unexplained Variability	(diagonal matrix)
  @annotated
  RUVPRO: 0.0876 : proportional RUV for plasma concentrations
  
$MAIN
// Pharmacokinetic Model Parameters
  double KA = POPKA*exp(ETA1);
  double CLF = POPCLF*pow(WT/70, COVWTCLF)*pow(AGE/45, COVAGECLF)*exp(ETA2);
  double VCF = POPVCF*pow(WT/70, COVWTVCF)*pow(AGE/45, COVAGEVCF)*exp(ETA3);
  double VPF = POPVPF*pow(WT/70, COVWTVPF);
  double QPF = POPQPF*pow(WT/70, COVWTQPF);
  
// Micro-rate Constants
  double K20 = CLF/VCF;
  double K23 = QPF/VCF;
  double K32 = QPF/VPF;
  
// Scaling Factor
  double SF = 1000;  // dose in mg, concentration in ng/mL
  
// mrgsolve Absorption Parameters
  ALAG_DEPOT = POPALAG;

$ODE
// Three-Compartment Pharmacokinetic Differential Equations
  dxdt_DEPOT = -KA*DEPOT;
  dxdt_PLAS = KA*DEPOT - K23*PLAS + K32*PERI - K20*PLAS;
  dxdt_PERI = K23*PLAS + -K32*PERI;

$TABLE  
// Individual Predictions
  double CPLAST = SF*PLAS/VCF;  // plasma
  
// Dependent Variable
  double DVPLAST = CPLAST*(1 + EPS(1));
  
$CAPTURE  
// Capture output
  DVPLAST // Dependent Variables
  CPLAST // Individual Predictions
  KA CLF VCF // Individual Parameters
  ETA1 ETA2 ETA3 // Inter-individual Variability  
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile model 
  modMuto <- mrgsolve::mcode("Muto2010_Azithromycin", code)

