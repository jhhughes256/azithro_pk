# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Model sourced from:
#   Zhao, Q., et al. (2014). "Population pharmacokinetics of azithromycin and 
#   chloroquine in healthy adults and paediatric malaria subjects following oral 
#   administration of fixed-dose azithromycin and chloroquine combination 
#   tablets." Malaria journal 13: 36-36.
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
  POPKA = 0.259,  // first-order absorption constant, hr^-1
  POPCLF = 100,  // apparent clearance, L/hr
  POPVCF = 186,  // apparent central volume, L
  POPVP1F = 2890,  // apparent peripheral volume (fast), L
  POPQP1F = 180,  // apparent intercompartment clearance (fast), L/hr
  POPVP2F = 2610,  // apparent peripheral volume (slow), L
  POPQP2F = 10.6,  // apparent intercompartment clearance (slow), L/hr
  POPRUVPRO = 0.406,  // proportional residual unexplained variability
  
// Default covariate values
  WT = 79,
    
// Inter-individual Variability
// ETA values from input dataset are used by default
  ETA1 = 0,  // Random effect for CLF
  ETA2 = 0,  // Random effect for VCF

$OMEGA  
// Intra Individual Variability (diagonal matrix)
  @annotated @block
  IIVCLF   : 0.097969 : intra-individual variability for CLF
  IIVVCF   : 0.28 1.2769 : intra-individual variability for VCF

$SIGMA  
// Residual Unexplained Variability	(diagonal matrix)
  @annotated
  RUVPRO: 1 : proportional RUV for plasma concentrations
  
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
  
// Scaling Factor
  double SF = 1000;  // dose in mg, concentration in ng/mL

$ODE
// Three-Compartment Pharmacokinetic Differential Equations
  dxdt_DEPOT = -KA*DEPOT;
  dxdt_PLAS = KA*DEPOT - K23*PLAS + K32*PER1 - K24*PLAS + K42*PER2 - K20*PLAS;
  dxdt_PER1 = K23*PLAS + -K32*PER1;
  dxdt_PER2 = K24*PLAS + -K42*PER2;

$TABLE  
// Individual Predictions
  double CPLAST = SF*PLAS/VCF;  // plasma (unbound)
  
// Dependent Variable
  double DVPLAST = CPLAST*exp(POPRUVPRO*EPS(1));
  
$CAPTURE  
// Capture output
  DVPLAST // Dependent Variables
  CPLAST // Individual Predictions
  CLF VCF // Individual Parameters
  ETA1 ETA2 // Inter-individual Variability  
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile model 
  modZhao <- mrgsolve::mcode("Zhao2014_Azithromycin", code)

