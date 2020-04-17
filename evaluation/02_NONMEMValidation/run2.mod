; Zheng et al. Azithromycin Pharmacokinetic model
; ------------------------------------------------------------------------------
; Model sourced from:
;   Songmao Zheng, Peter Matzneller, Markus Zeitlinger, Stephan Schmidt, 
;   Development of a Population Pharmacokinetic Model Characterizing the Tissue
;   Distribution of Azithromycin in Healthy Subjects, Antimicrobial Agents and 
;   Chemotherapy Oct 2014, 58 (11) 6675-6684; DOI: 10.1128/AAC.02904-14
; ------------------------------------------------------------------------------
; PK concentrations (ug/mL)
; AMT in mg
; TIME in hours
$PROBLEM run2.mod
; Base model (unbound plasma concentrations)
;   3-compartment model
;   First-order absorption and absorption lag
;   First-order elimination
;   Random effects on ALAG1, CL and V2
;   Combined error model
; Expanded model (unbound/unionised tissue concentrations)
;   Three tissue compartments + three deep-tissue compartments
;   Distribution into and out of tissue (KIN-KOUT) 
;   Tissue binding within each tissue (KON-KOFF)
;   Distribution factor estimated for each tissue (DF...)
;   Unionised plasma available for distribution to PML tissue (FUP)
;   Random effects on KIN and three DF... parameters
;   Combined error model for each tissue
; Dependent Variables
;   CMT2: Unbound plasma concentrations
;   CMT3: Unbound muscle interstitial fluid concentrations
;   CMT4: Unbound subcutis interstitial fluid concentrations
;   CMT9: Unionised polymorphonuclear leukocyte cytosol concentrations

$INPUT C,ID,TIME,AMT,CMT,EVID,MDV,DV,DVID
       
$DATA AZPop_30MAR2020.csv IGNORE = C

$SUBROUTINE ADVAN13 TOL=9

$MODEL
COMP=(1_DEPOT)  
COMP=(2_PLASMA)
COMP=(3_MUSCLE)
COMP=(4_SUBCU)
COMP=(5_OTHERF)
COMP=(6_OTHERS)
COMP=(7_MUS1)
COMP=(8_SUB1)
COMP=(9_LEUKO)
COMP=(10_LEU1)
;COMP=(9_DEPOT)  
;COMP=(10_PLASMA)
;COMP=(11_MUSCLE)
;COMP=(12_SUBCU)
;COMP=(13_OTHERF)
;COMP=(14_OTHERS)
;COMP=(15_MUS1)
;COMP=(16_SUB1)
;COMP=(17_LEUKO)
;COMP=(18_LEU1)

$PK ; Define basic PK relationships

KA = THETA(1)*EXP(ETA(1))
CL = THETA(2)*EXP(ETA(2))
V2 = THETA(3)*EXP(ETA(3))

KIN = THETA(4)*EXP(ETA(4)) ; MUSCLE and subcutis
KOUT = THETA(5)*EXP(ETA(5)) ; MUSCLE and subcutis

TVLAG=THETA(6)*EXP(ETA(6))
ALAG1=TVLAG

KPMU = THETA(7)*EXP(ETA(7)); MUSCLE
KPSU = THETA(8)*EXP(ETA(8)); SUBCUTIS

V5 = THETA(9)*EXP(ETA(9)) ; volume for other organs FAST
Q5 = THETA(10)*EXP(ETA(10)) ;plasma flow for other organs FAST

V6 = THETA(11)*EXP(ETA(11)); volume for other organs SLOW
Q6 = THETA(12)*EXP(ETA(12)); plasma flow for other organs SLOW

KON= THETA(13)*EXP(ETA(13)); BINDING/OTHER PROCESS IN TISSUE THAT DECREASE UNBOUND CONC. IN EXCHANGE WITH PLASMA
KOFF = THETA(14)*EXP(ETA(14)); OPPOSITE PROCESS RELATED TO KON

KPLE=THETA(15)*EXP(ETA(15))

K25  = Q5/V2
K52  = Q5/V5

K26  = Q6/V2
K62  = Q6/V6

K20  = CL/V2

FUP = 0.007581

;S2=V2
;S3= KIN*V2/(KOUT*KPMU); KPMU  PARTITION COEFFIENT FOR MUSLE
;S4= KIN*V2/(KOUT*KPSU); KPSU PARTITION COEFFIENT FOR SUBCUTIS
;S9= KIN*V2/(KOUT*KPLE); KPLE PARTITION COEFFIENT FOR LEUKOCYTE

$DES ; Define differential equations
;Dosing Compartment
DADT(1) = -KA*A(1)

;Plasma (UNBOUND)
DADT(2) = KA*A(1)-(K20+KIN+KIN+KIN*FUP+K26+K25)*A(2)+KOUT*A(3)+KOUT*A(4)+KOUT*A(9)+K52*A(5)+K62*A(6)

;muscle unbound
DADT(3)=KIN*A(2)+KOFF*A(7)-KOUT*A(3)-KON*A(3)

;subcutis unbound
DADT(4)=KIN*A(2)+KOFF*A(8)-KOUT*A(4)-KON*A(4)

;other FAST
DADT(5)=K25*A(2)-K52*A(5)

;other slow
DADT(6)=K26*A(2)-K62*A(6)

DADT(7)=KON*A(3)-KOFF*A(7)

DADT(8)=KON*A(4)-KOFF*A(8)

;LEUKO UNIONIZED
DADT(9)=KIN*FUP*A(2)+KOFF*A(10)-KOUT*A(9)-KON*A(9)

DADT(10)=KON*A(9)-KOFF*A(10)

$ERROR

; Not original code - Added by JIM for debug
APK = A(2)
AMU = A(3)
ASU = A(4)
ALE = A(9)
;

IPRED=0
I1=0
I2=0
I3=0
I4=0

IF(CMT==2) THEN
I1=1
IPRED=A(2)/V2
ENDIF

IF(CMT==3) THEN
I2=1
IPRED=A(3)*KOUT*KPMU/(KIN*V2)
ENDIF

IF(CMT==4) THEN
I3=1
IPRED=A(4)*KOUT*KPSU/(KIN*V2)
ENDIF

IF(CMT==9) THEN
I4=1
IPRED=A(9)*KOUT*KPLE/(KIN*V2)
ENDIF

Y= IPRED*(1+EPS(1)*I1+EPS(3)*I2+EPS(5)*I3+EPS(7)*I4)+(EPS(2)*I1+EPS(4)*I2+EPS(6)*I3+EPS(8)*I4)
IRES=IPRED-Y

$THETA
  0.88 FIX ; POPKA
  258 FIX ; POPCL
  160 ; POPV2
  0.16 ; POPKIN
  0.15 ; POPKOUT
  1.45 FIX ; POPTLAG
  0.55 ; POPKPMU
  0.25 ; POPKPSU
  1190 ; POPV5
  207 ; POPQ5
  9721 ; POPV6
  101 ; POPQ6
  0.56 ; POPKON
  0.05 ; POPKOFF
  52 ; POPKPLE

$OMEGA
  0 FIX ; PPVKA
  0.085849 ; PPVCL
  2.832489 ; PPVV2
  0.00000484 ; PPVKIN
  0 FIX ; PPVKOUT
  0.030976 ; PPVTLAG
  0.072361 ; PPVKPMU
  0.099225 ; PPVKPSU
  0 FIX ; PPVV5
  0 FIX ; PPVQ5
  0 FIX ; PPVV6
  0 FIX ; PPVQ6
  0 FIX ; PPVKON
  0 FIX ; PPVKOFF
  0.00000484 ; PPVKPLE

$SIGMA
  0.14 ; ERRPKPRO
  35.2 ; ERRPKADD
  0.14 ; ERRMUPRO
  0.51 ; ERRMUADD
  0.34 ; ERRSCPRO
  0.000001 ; ERRSCADD
  0.23 ; ERRLEPRO
  0.000001 ; ERRLEADD

$SIMULATION(555) ONLYSIM SUBPROB=1

$TABLE ID TIME AMT CMT EVID MDV CL V2 KIN TVLAG KPMU KPSU KPLE
  ETAS(1:LAST) IPRED PRED IRES APK AMU ASU ALE
  NOPRINT ONEHEADER FILE = run2.fit
; ------------------------------------------------------------------------------

