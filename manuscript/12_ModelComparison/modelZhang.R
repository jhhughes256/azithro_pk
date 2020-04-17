# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Model sourced from:
#   Zhang, X. Y., et al. (2010). "Population pharmacokinetics study of 
#   azithromycin oral formulations using NONMEM." Int J Clin Pharmacol 
#   Ther 48(10): 662-669.
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
  POPKA = 1.09,  // first-order absorption constant, hr^-1
  POPCLF = 121,  // apparent clearance, L/hr
  POPVCF = 1939,  // apparent central volume, L
  POPVPF = 5650,  // apparent peripheral volume (fast), L
  POPQPF = 282,  // apparent intercompartment clearance (fast), L/hr
  POPRUVPRO = 0.3296,  // proportional RUV, SD
  POPRUVADD = 0.01,  // additive RUV, SD (ng/mL)
  
// Covariate-Effect Parameters
  COVWTCLF = 0.0379,  // Effect of body weight on CLF
  COVAGEVCF = 0.161,  // Effect of age on VCF
  
// Default covariate values
  WT = 60.43,  // mean body weight from Zhang et al.
  AGE = 19.03,  // mean age from Zhang et al.
    
// Inter-individual Variability
// ETA values from input dataset are used by default
  ETA1 = 0,  // Random effect for CLF
  ETA2 = 0,  // Random effect for QPF
  ETA3 = 0,  // Random effect for VCF
  ETA4 = 0,  // Random effect for VPF
  ETA5 = 0,  // Random effect for KA

$OMEGA  
// Intra Individual Variability (diagonal matrix)
  @annotated 
  IIVCLF : 0.02128681 : intra-individual variability for CLF
  IIVQPF : 0.08838729 : intra-individual variability for QPF
  IIVVCF : 0.00000001 : intra-individual variability for VCF
  IIVVPF : 0.19027044 : intra-individual variability for VPF
  IIVKA  : 0.69288976 : intra-individual variability for KA

$SIGMA  
// Residual Unexplained Variability	(diagonal matrix)
  @annotated
  RUVPRO: 0.10863616 : proportional RUV for plasma concentrations
  RUVADD: 0.0001 : additive RUV for plasma concentrations
  
$MAIN
// Pharmacokinetic Model Parameters
  double KA = POPKA*exp(ETA5);
  double CLF = (POPCLF - COVWTCLF*(60.43 - WT))*exp(ETA1);
  double VCF = (POPVCF - COVAGEVCF*(AGE - 19.03))*exp(ETA3);
  double VPF = POPVPF*exp(ETA4);
  double QPF = POPQPF*exp(ETA2);
  
// Micro-rate Constants
  double K20 = CLF/VCF;
  double K23 = QPF/VCF;
  double K32 = QPF/VPF;
  
// Scaling Factor
  double SF = 1000;  // dose in mg, concentration in ng/mL

$ODE
// Three-Compartment Pharmacokinetic Differential Equations
  dxdt_DEPOT = -KA*DEPOT;
  dxdt_PLAS = KA*DEPOT - K23*PLAS + K32*PERI - K20*PLAS;
  dxdt_PERI = K23*PLAS + -K32*PERI;

$TABLE  
// Individual Predictions
  double CPLAST = SF*PLAS/VCF;  // plasma (unbound)
  
// Dependent Variable
  double DVPLAST = CPLAST*(1 + EPS(1)) + EPS(2);
  
$CAPTURE  
// Capture output
  DVPLAST // Dependent Variables
  CPLAST // Individual Predictions
  CLF QPF VCF VPF KA // Individual Parameters
  ETA1 ETA2 ETA3 ETA4 ETA5 // Inter-individual Variability  
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile model 
  modZhang <- mrgsolve::mcode("Zhang2010_Azithromycin", code)

