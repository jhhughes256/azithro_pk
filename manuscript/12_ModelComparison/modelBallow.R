# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Model sourced from:
#   Ballow, C. H., et al. (1998). "Pharmacokinetics of Oral Azithromycin in 
#   Serum, Urine, Polymorphonuclear Leucocytes and Inflammatory vs 
#   Non-Inflammatory Skin Blisters in Healthy Volunteers." Clinical Drug 
#   Investigation 15(2): 159-167.
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

$PARAM
// Fixed-Effects Parameters
  POPALAG = 1.11,  // absorption lag, hr
  POPKA = 0.619,  // first-order absorption constant, hr^-1
  POPCLF = 143,  // apparent clearance, L/hr
  POPCLR = 10.6,  // renal clearance, L/hr
  POPVCF = 597,  // apparent central volume, L
  POPVP1F = 4305,  // apparent peripheral volume (fast), L
  POPQP1F = 253,  // apparent intercompartment clearance (fast), L/hr
  POPVP2F = 3709,  // apparent peripheral volume (slow), L
  POPQP2F = 61.3,  // apparent intercompartment clearance (slow), L/hr
  
$MAIN
// Pharmacokinetic Model Parameters
  double KA = POPKA;
  double CLF = POPCLF;
  double VCF = POPVCF;
  double VP1F = POPVP1F;
  double QP1F = POPQP1F;
  double VP2F = POPVP2F;
  double QP2F = POPQP2F;
  
// Micro-rate Constants
  double K20 = CLF/VCF;
  double K23 = QP1F/VCF;
  double K32 = QP1F/VP1F;
  double K24 = QP2F/VCF;
  double K42 = QP2F/VP2F;
  
// Scaling Factor
  double SF = 1000;  // dose in mg, concentration in ng/mL
  
// mrgsolve Absorption Parameters
  ALAG_DEPOT = POPALAG;

$ODE
// Three-Compartment Pharmacokinetic Differential Equations
  dxdt_DEPOT = -KA*DEPOT;
  dxdt_PLAS = KA*DEPOT - K23*PLAS + K32*PER1 - K24*PLAS + K42*PER2 - K20*PLAS;
  dxdt_PER1 = K23*PLAS + -K32*PER1;
  dxdt_PER2 = K24*PLAS + -K42*PER2;

$TABLE  
// Individual Predictions
  double CPLAST = SF*PLAS/VCF;  // plasma
  
// Dependent Variable
  double DVPLAST = CPLAST;
  
$CAPTURE  
// Capture output
  DVPLAST // Dependent Variables
  CPLAST // Individual Predictions
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile model 
  modBall <- mrgsolve::mcode("Ballow1998_Azithromycin", code)

