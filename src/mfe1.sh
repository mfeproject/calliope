#!/bin/bash

#get the real location of this script
THIS_SCRIPT=$(readlink -f $(type -p $0))

MFE_DIR=${THIS_SCRIPT%bin/*}
MFE_LIBDIR=${MFE_DIR}/lib
MFE_EXE=${MFE_LIBDIR}/mfe1

exec ${MFE_EXE} -L ${MFE_LIBDIR} $@
