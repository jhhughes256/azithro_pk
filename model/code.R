# Azithromycin Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Model sourced from:
#   Songmao Zheng, Peter Matzneller, Markus Zeitlinger, Stephan Schmidt, 
#   Development of a Population Pharmacokinetic Model Characterizing the Tissue
#   Distribution of Azithromycin in Healthy Subjects, Antimicrobial Agents and 
#   Chemotherapy Oct 2014, 58 (11) 6675-6684; DOI: 10.1128/AAC.02904-14
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

$PARAM
// Fixed-Effects Parameters
  POPTLAG = 1.45,  // absorption lag, hr
  POPKA = 0.88,  // first-order absorption constant, hr^-1
  POPCLF = 258,  // apparent clearance, L/hr
  POPVCF = 160,  // apparent central volume, L
  POPVP1F = 1190,  // apparent peripheral volume (fast), L
  POPQP1F = 207,  // apparent intercompartment clearance (fast), L/hr
  POPVP2F = 9721,  // apparent peripheral volume (slow), L
  POPQP2F = 101,  // apparent intercompartment clearance (slow), L/hr
  POPKIN = 0.16,  // uptake into tissue (rate constant), hr^-1
  POPKOUT = 0.15,  // distribution out of tissue (rate constant), hr^-1
  POPKON = 0.56,  // non-specific tissue binding (on rate constant), hr^-1
  POPKOFF = 0.05,  // non-specific tissue binding (off rate constant), hr^-1
  POPDFMUSC = 0.55,  // distribution factor (muscle), dimensionless
  POPDFSUBC = 0.25,  // distribution factor (subcutis), dimensionless
  POPDFPMLC = 52,  // distribution factor (PML cytosol), dimensionless
  
// Literature/Derived Values
  FUNPLAS = 0.0076,  // fraction unionised in plasma, dimensionless
  FUNPMLC = 0.0012,  // fraction unionised in PML cytosol, dimensionless
  RATPMLC = 0.351,  // ratio of PML cytosol to total PML conc., dimensionless
  RATPMLL = 13.33,  // ratio of PML lysosome to total PML conc., dimensionless
    
// Inter-individual Variability
// ETA values from input dataset are used by default
  ETA1 = 0,  // Random effect for TLAG
  ETA2 = 0,  // Random effect for CLF
  ETA3 = 0,  // Random effect for VCF
  ETA4 = 0,  // Random effect for KIN
  ETA5 = 0,  // Random effect for DFMUSC
  ETA6 = 0,  // Random effect for DFSUBC
  ETA7 = 0  // Random effect for DFPMLC

$OMEGA  
// Intra Individual Variability (diagonal matrix)
  @annotated 
  IIVTLAG  : 0.030976 : intra-individual variability for TLAG
  IIVCLF   : 0.085849 : intra-individual variability for CLF
  IIVVCF   : 2.832489 : intra-individual variability for VCF
  IIVKIN   : 0.00000484 : intra-individual variability for KIN
  IIVDFMUSC: 0.072361 : intra-individual variability for DFMUSC
  IIVDFSUBC: 0.099225 : intra-individual variability for DFSUBC
  IIVDFPMLC: 0.00000484 : intra-individual variability for DFPMLC

$SIGMA  
// Residual Unexplained Variability	(diagonal matrix)
  @annotated
  RUVPLASPRO: 0.14 : proportional RUV for plasma concentrations
  RUVPLASADD: 35.2 : proportional RUV for plasma concentrations
  RUVMUSCPRO: 0.14 : proportional RUV for plasma concentrations
  RUVMUSCADD: 0.51 : proportional RUV for plasma concentrations
  RUVSUBCPRO: 0.34 : proportional RUV for plasma concentrations
  RUVSUBCADD: 0.000001 : proportional RUV for plasma concentrations
  RUVPMLCPRO: 0.23 : proportional RUV for plasma concentrations
  RUVPMLCADD: 0.000001 : proportional RUV for plasma concentrations

$MAIN
// Pharmacokinetic Model Parameters
  double TLAG = POPTLAG*exp(ETA1);
  double KA = POPKA;
  double CLF = POPCLF*exp(ETA2);
  double VCF = POPVCF*exp(ETA3);
  double VP1F = POPVP1F;
  double QP1F = POPQP1F;
  double VP2F = POPVP2F;
  double QP2F = POPQP2F;
  
// Tissue Kinetic Parameters
  double KIN = POPKIN*exp(ETA4);
  double KOUT = POPKOUT;
  double KON = POPKON;
  double KOFF = POPKOFF;
  
// Distribution Factors
  double DFMUSC = POPDFMUSC*exp(ETA5);
  double DFSUBC = POPDFSUBC*exp(ETA6);
  double DFPMLC = POPDFPMLC*exp(ETA7);
  
// Scaling Factor
  double SF = 1000;  // dose in mg, concentration in ng/mL
  
// mrgsolve Absorption Parameters
  ALAG_DEPOT = TLAG;

$ODE
// Unionised drug in plasma
  double PLASU = PLAS*FUNPLAS;

// Three-Compartment Pharmacokinetic Differential Equations
  dxdt_DEPOT = -KA*DEPOT;
  dxdt_PLAS = KA*DEPOT + QP1F*PER1/VP1F - QP1F*PLAS/VCF + QP2F*PER2/VP2F -  
    QP2F*PLAS/VCF - KIN*PLAS + KOUT*MUSC - KIN*PLAS + KOUT*SUBC - 
    KIN*PLASU + KOUT*PMLC - CLF*PLAS/VCF;
  dxdt_PER1 = QP1F*PLAS/VCF + -QP1F*PER1/VP1F;
  dxdt_PER2 = QP2F*PLAS/VCF + -QP2F*PER2/VP2F;
  
// Tissue Pharmacokinetic Differential Equations
  dxdt_MUSC = KIN*PLAS - KOUT*MUSC + KOFF*MUSCB - KON*MUSC;
  dxdt_SUBC = KIN*PLAS - KOUT*SUBC + KOFF*SUBCB - KON*SUBC;
  dxdt_PMLC = KIN*PLASU - KOUT*PMLC + KOFF*PMLCB - KON*PMLC;
  dxdt_MUSCB = KON*MUSC - KOFF*MUSCB;
  dxdt_SUBCB = KON*SUBC - KOFF*SUBCB;
  dxdt_PMLCB = KON*PMLC - KOFF*PMLCB;

$TABLE  
// Plasma Concentrations
  double CPLAS = SF*PLAS/VCF;  // plasma (unbound)
// fraction unbound is saturable for azithromycin
// determined by FU = 0.4984+0.5339*CPLAST/(230.9+CPLAST)
// substitute CPLAST = CPLAS/FU then solve for FU:
// Method: clear denominators, transpose to one side, quadratic formula
  double FU = (-(CPLAS - 0.4984*230.9) + sqrt(pow((CPLAS - 0.4984*230.9), 2) - 
    4*230.9*-(0.4984+0.5339)*CPLAS))/(2*230.9);  // fraction unbound
  double CPLAST = CPLAS/FU;  // plasma (total)
  
// Tissue Concentrations
  double CMUSC = SF*MUSC*KOUT*DFMUSC/(KIN*VCF);  // muscle (unbound)
  double CSUBC = SF*SUBC*KOUT*DFSUBC/(KIN*VCF);  // subcutis (unbound)
  
// PML Concentrations
  double CPMLC = SF*PMLC*KOUT*DFPMLC/(KIN*VCF);  // PML cytosol (unionised)
  double CPMLCT = CPMLC/FUNPMLC;  // PML cytosol (total)
  double CPMLT = CPMLCT/RATPMLC;  // PML (total)
  double CPMLLT = RATPMLL*CPMLT;  // PML lysosome (total)
  
$CAPTURE  
// Capture output
  FU CPLAS CPLAST CMUSC CSUBC CPMLC CPMLCT CPMLT CPMLLT // Outputs
//  TLAG CLF VCF KIN DFMUSC DFSUBC DFPMLC  // Individual Parameters
//  ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7  // Inter-individual Variability  
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile model 
  mod <- mrgsolve::mcode("Zheng2014_Azithromycin", code)

