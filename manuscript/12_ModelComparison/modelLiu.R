# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Model sourced from:
#   Liu, P., et al. (2007). "Comparative pharmacokinetics of azithromycin in 
#   serum and white blood cells of healthy subjects receiving a single-dose 
#   extended-release regimen versus a 3-day immediate-release regimen." 
#   Antimicrob Agents Chemother 51(1): 103-109.
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
  PLAS = 0,  // unbound plasma
  PERI = 0,  // fast-distribution peripheral

$PARAM
// Fixed-Effects Parameters
  POPR1 = 184,  // zero-order absorption constant, mg/hr
  POPCLF = 122,  // apparent clearance, L/hr
  POPVCF = 518,  // apparent central volume, L
  POPVPF = 3400,  // apparent peripheral volume (fast), L
  POPQPF = 102,  // apparent intercompartment clearance (fast), L/hr
    
// Inter-individual Variability
// ETA values from input dataset are used by default
  ETA1 = 0,  // Random effect for CLF
  ETA2 = 0,  // Random effect for VCF
  ETA3 = 0,  // Random effect for VPF
  ETA4 = 0,  // Random effect for QPF
  ETA5 = 0,  // Random effect for R1

$OMEGA  
// Intra Individual Variability (diagonal matrix)
  @annotated
  IIVCLF : 0.109 : intra-individual variability for CLF
  IIVVCF : 0.992 : intra-individual variability for VCF
  IIVVPF : 0.00382 : intra-individual variability for VPF
  IIVQPF : 0.233 : intra-individual variability for QPF
  IIVR1  : 0.0108 : intra-individual variability for R1

$SIGMA  
// Residual Unexplained Variability	(diagonal matrix)
  @annotated
  RUVPRO: 0.0881 : proportional RUV for plasma concentrations
  
$MAIN
// Pharmacokinetic Model Parameters
  double R1 = POPR1*exp(ETA5);
  double CLF = POPCLF*exp(ETA1);
  double VCF = POPVCF*exp(ETA2);
  double VPF = POPVPF*exp(ETA3);
  double QPF = POPQPF*exp(ETA4);
  
// Micro-rate Constants
  double K20 = CLF/VCF;
  double K23 = QPF/VCF;
  double K32 = QPF/VPF;
  
// Scaling Factor
  double SF = 1000;  // dose in mg, concentration in ng/mL
  
// mrgsolve Absorption Parameters
  R_PLAS = R1;

$ODE
// Three-Compartment Pharmacokinetic Differential Equations
  dxdt_PLAS = -K23*PLAS + K32*PERI - K20*PLAS;
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
  CLF VCF VPF QPF R1 // Individual Parameters
  ETA1 ETA2 ETA3 ETA4 ETA5 // Inter-individual Variability  
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile model 
  modLiu <- mrgsolve::mcode("Liu2007_Azithromycin", code)

