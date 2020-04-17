# Azithromycin Concentrations in Cells, Epithelial Lining Fluid and Serum
# -----------------------------------------------------------------------------
# Data digitised from:
#   Olsen, K. M., et al. (1996). "Intrapulmonary pharmacokinetics of 
#   azithromycin in healthy volunteers given five oral doses." Antimicrobial 
#   Agents and Chemotherapy 40(11): 2582.
# -----------------------------------------------------------------------------

olsendf <- tibble::tribble(
  ~Sample, ~MN_4, ~SD_4, ~MN_28, ~SD_28, ~MN_76, ~SD_76, ~MN_124, ~SD_124, ~MN_172, ~SD_172, ~MN_244, ~SD_244, ~MN_340, ~SD_340, ~MN_508, ~SD_508,
     "AM",    74,    26,    180,    105,    217,     64,     464,      65,     304,      71,     112,      20,     110,      40,      82,      88,
    "PBM",    58,    13,    120,     28,    110,     31,     107,      18,      77,      22,      47,      12,      21,      14,       0,       0,
    "ELF",  0.45,  0.15,   1.53,   0.31,   2.67,   0.85,    3.12,    0.93,    0.61,    0.23,       0,       0,       0,       0,       0,       0, 
  "Serum", 0.178,  0.05,  0.122,  0.055,  0.093,  0.036,   0.054,   0.008,   0.031,   0.055,   0.015,   0.005,       0,       0,       0,       0
)

# -----------------------------------------------------------------------------
# AM - alveolar macrophage
# PBM - peripheral blood monocytes
# ELF - Epithelial lining fluid (lung)
# Number of subjects was five at 28, 76, 124 and 508 hours
# All other time points had six subjects
# Drug concentrations are in ug/mL for ELF and Serum.
# Drug concentrations are in ug/mL of monocyte/macrophage volume for AM and PBM.
