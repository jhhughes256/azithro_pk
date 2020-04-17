# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Model sourced from:
#   Salman, S., et al. (2010). "Pharmacokinetic Properties of Azithromycin in 
#   Pregnancy." Antimicrobial Agents and Chemotherapy 54(1): 360.
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
  POPDUR = 1.55,  // zero-order absorption into DEPOT, hr
  POPKA = 0.525,  // first-order absorption constant, hr^-1
  POPCLF = 158,  // apparent clearance, L/hr
  POPVCF = 384,  // apparent central volume, L
  POPVP1F = 4080,  // apparent peripheral volume (fast), L
  POPQP1F = 325,  // apparent intercompartment clearance (fast), L/hr
  POPVP2F = 5040,  // apparent peripheral volume (slow), L
  POPQP2F = 66.4,  // apparent intercompartment clearance (slow), L/hr
  POPRUVPRO = 0.3296,  // proportional RUV, SD
  POPRUVADD = 0.01,  // additive RUV, SD (ng/mL)
  
// Covariate-Effect Parameters
  COVPREGVCF = 330,  // Effect of pregnancy on VCF
  
// Default covariate values
  WT = 51.7,  // mean body weight from Salman et al.
  PREG = 0,  // pregnancy status
    
// Inter-individual Variability
// ETA values from input dataset are used by default
  ETA1 = 0,  // Random effect for VCF
  ETA2 = 0,  // Random effect for CLF
  ETA3 = 0,  // Random effect for VPF
  ETA4 = 0,  // Random effect for DUR

$OMEGA  
// Intra Individual Variability (diagonal matrix)
  @annotated 
  IIVVCF : 0.992016 : intra-individual variability for VCF
  IIVCLF : 0.080089 : intra-individual variability for CLF
  IIVVPF : 0.126736 : intra-individual variability for VPF
  IIVDUR : 0.591361 : intra-individual variability for DUR

$SIGMA  
// Residual Unexplained Variability	(diagonal matrix)
  @annotated
  RUVPRO: 0.097344 : proportional RUV for plasma concentrations
  
$MAIN
// Pharmacokinetic Model Parameters
  double DUR = POPDUR*exp(ETA4);
  double KA = POPKA;
  double CLF = POPCLF*pow(WT/70, 0.75)*exp(ETA2);
  double VCF = (POPVCF+COVPREGVCF)*pow(WT/70, 1.00)*exp(ETA1);
  double VP1F = POPVP1F*pow(WT/70, 1.00)*exp(ETA3);
  double QP1F = POPQP1F*pow(WT/70, 0.75);
  double VP2F = POPVP2F*pow(WT/70, 1.00);
  double QP2F = POPQP2F*pow(WT/70, 0.75);
  
// Micro-rate Constants
  double K20 = CLF/VCF;
  double K23 = QP1F/VCF;
  double K32 = QP1F/VP1F;
  double K24 = QP2F/VCF;
  double K42 = QP2F/VP2F;
  
// Scaling Factor
  double SF = 1000;  // dose in mg, concentration in ng/mL
  
// mrgsolve Absorption Parameters
  D_DEPOT = DUR;

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
  double DVPLAST = CPLAST*(1 + EPS(1));
  
$CAPTURE  
// Capture output
  DVPLAST // Dependent Variables
  CPLAST // Individual Predictions
  VCF CLF VP1F DUR // Individual Parameters
  ETA1 ETA2 ETA3 ETA4 // Inter-individual Variability  
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile model 
  modSalm <- mrgsolve::mcode("Salman2010_Azithromycin", code)

