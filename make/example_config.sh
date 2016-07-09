
# Example config script. Default is all engines disabled. Before 
# enabling engines, make sure environment is set correctly in 
# example_setup.sh
#
#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/setup.sh

# List of configure options
# -- -- sub level denotes extra flags to be used with the first level flag

# --en/disable-t2kreweight

# --en/disable-neut
# -- --with-cern=$CERN_ROOT \

# --en/disable-niwg \    

# --en/disable-genie \
# -- --with-pythia6-lib=$PYTHIA6 \
# -- --with-libxml2-inc=$LIBXML2_INC \
# -- --with-libxml2-lib=$LIBXML2_LIB \
# -- --with-lhapdf-inc=$LHAPDFINC  \
# -- --with-lhapdf-lib=$LHAPDFLIB \
# -- --with-log4cpp-inc=$LOG4CPP_INC \
# -- --with-log4cpp-lib=$LOG4CPP_LIB \

# --en/disable-nuwro \
# -- --with-pythia6-lib=$PYTHIA6 \
# -- --en/disable-nuwro-reweight \

# For help with configure flags run
# ./configure --help

# Default builds just against NEUT + NIWG + T2KRW
$EXT_FIT/configure \
    --enable-t2kreweight \
    --enable-neut \
    --with-cern=$CERN_ROOT \
    --enable-niwg \
    --disable-genie \
    --enable-nuwro \
    --enable-nuwro-reweight \
