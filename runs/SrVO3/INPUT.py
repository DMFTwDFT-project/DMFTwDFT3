
################   Input parameters for DFT+DMFT calculations   ##################

###### Main loop parameters ###########

p = {"Niter":     5,               # Number of DFT+DMFT iterations
     "Nit":       1,               # Number of DMFT iterations
     "n_tot":     19,            # Number of total electrons
     "nspin":     1,            # Number of total electrons
     "atomnames": ['V','O'],       # The name of atoms"],
     "orbs"     : ['d','p'],       # The name  of orbitals
     "L_rot"     : [0,0],       # Whether rotate local axis or not
     "cor_at":    [['V1']],      # Correlated atoms, put degenerate atoms in the same list
     "U":         [5.0],        # Intra-U for each cor_at
     "J":         [1.0],        # Hund's coupling"],
     "Uprime":    [4.8],        # Double counting U'"],
     "cor_orb":   [[['d_z2','d_x2y2'],['d_xz','d_yz','d_xy']]], # DMFT orbitals, other orbitals are treated by HF"],
     "q":         [20,20,20],       # [Nq_x,Nq_y,Nq_z]"],
     "ewin":      [-7,7],           # Energy Window
     "nom":       6000,             # Number of Matsubara frequencies"],
     "noms":      1200,            # Number of Matsubara frequencies"],
     "nomlog":    30,             # Number of Matsubara frequencies"],
     "Nforce":    0,             # How many time force will be computed"],
     "nfine":     1,             # fine mesh"],
     "T_high":    0.01,             # fine mesh"],
     "dc_type":   1,             # dc type"],
     "mu_iter":   200,                # The chemical potential convergence step"],
     "mix_sig":   0.3,              # Mixing parameter for Sigma"],
     "Nd_qmc":    1,             # DMFT Nd values are obtained from QMC sampling"],
     "print_at":  ['V1'],     # The local Green functions are printed"],
     "Nrelax" : 0,
     "path_bin":  "/users/ukh0001/projects/DFTDMFT/bin/",     # Path to bin files
     }

##### VASP parameters ########
pV = {"System=":     ["SrVO3",       "# The name of system"],
      "ENCUT=":      [600.0,                 "# Energy cutoff"], 
      "ISPIN=":      [1,               "#ISPIN"],
      "NBANDS=":     [46,                     "# LMAX"],
      "LMAXMIX=":    [4,             "#LMAX"],
      "NCORE=":    [1,             "#LMAX"],
      "IALGO=":    [48,             "#LMAX"],
      #"NELM=":    [20,             "#LMAX"],
      "ISMEAR=":      [-4,                    "# ISMEAR"],
      "EDIFF=":      [1e-4,                    "# ISMEAR"],
      "GGA=":    ["PS",             "#GGA"],
      "ADDGRID=":    [".TRUE.",             "#GGA"],
      "LWANNIER90=":      [".TRUE.",                    "# ISMEAR"],
      #"LDMFT=":      [".TRUE.",                    "#  DMFT"],
     } 

##### CTQMC parameters #########
pC = {"exe":         ["ctqmc",    "# Path to executable"],
      "Delta":       ["Delta.inp",                 " # Input bath function hybridization"],
      "cix":         ["impurity.cix",              " # Input file with atomic state"],
      "mu":          [0,                           " # Chemical potential"],
      "beta":        [100.0,                       " # Inverse temperature"],
      "M" :          [5000000,                    " # Number of Monte Carlo steps"],
      "nom":         [60,                          " # number of Matsubara frequency points to sample"],
      "Nmax":        [1400,                        " # maximal number of propagators"],
      "Ntau":        [1000,                        " # Ntau"],
      "SampleGtau":  [1000,                        " # Sample Gtau"],
      "sderiv":      [0.005,                       " # maximal discripancy"],
      "CleanUpdate": [50000,                       " # clean update after QMC steps"],
      "aom":         [8,                           " # number of frequency points to determin high frequency tail"],
      "warmup":      [250000,                      "  # Warmup"],
      "GlobalFlip":  [1000000,                      "  # Global flip"],
      "PChangeOrder":[0.9,                         "  # Probability to add/remove interval"],
      "TwoKinks":    [0.,                          "  # Two kinks"],
      "Ncout":       [500000,                      " # Ncout"],
      "Naver":       [80000000,                    "  # Naver"],
      }
############## CIX parameters ###########
pD = {"para="     : 1,
      "l="        : 2,
      "n="        : [6, 7, 8],
      "OCA_G="    : False,
      "CoulombF=" : "Ising",
      }

