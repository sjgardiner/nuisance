
#--------------------
# Fitter Requirements
#-------------------- 

# ROOT Dependency
export ROOTSYS="/home/stowell/software/root/"

# dependency on NIWG_DATA will eventually be removed
export NIWG_DATA="/home/stowell/t2krep/NIWG/external_data/" #path to external data    

# -------------------------
# Generator/RW requirements
# -------------------------

#GENIE
export GENIE="" #/path/to/genie/

# NUWRO
export NUWRO="/home/stowell/software/nuwro/" # /path/to/nuwro

# Neut
export NEUT="/home/stowell/t2krep/NIWG/neut/branches/neut_5.3.4_v1r25/" #/path/to/neut

# NIWG
export NIWGREWEIGHT="/home/stowell/t2krep/GlobalAnalysisTools/NIWGReWeight/head/" #/path/to/niwgreweight

# T2K RW Stuff
export T2KREWEIGHT="/home/stowell/t2krep/GlobalAnalysisTools/T2KReWeight/head/" # path to t2krew

#-------------------
# Other Dependencies
#-------------------

# Cern LIB  (Required by NEUT and T2KREW)
export CERN="/home/stowell/software/CERNLIB"
export CERN_LEVEL=2005

# PYTHIA (Required by GENIE and NuWro)
export PYTHIA6=""

# LIBXML2 (Required By GENIE)
export LIBXML2_INC=""
export LIBXML2_LIB=""

#LHAPDFINC (Required By GENIE)
export LHAPDFINC=""
export LHAPDFLIB=""

#LOG4CPP (Required by GENIE)
export LOG4CPP_INC=""
export LOG4CPP_LIB=""

# LHAPATH (Required by GENIE)
export LHAPATH=""
