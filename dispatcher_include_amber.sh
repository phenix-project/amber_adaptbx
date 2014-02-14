#
# Environment additions for AMBER interface
#

# include before command
if [ ! -z "$AMBERHOME" ]; then
  if [ "$LIBTBX_OS_NAME" = "Darwin" ]; then
    if [ -z "$DYLD_LIBRARY_PATH" ]; then
      export DYLD_LIBRARY_PATH=${AMBERHOME}/lib
    else
      export DYLD_LIBRARY_PATH=${AMBERHOME}/lib:${DYLD_LIBRARY_PATH}
    fi
  else
    if [ -z "$LD_LIBRARY_PATH" ]; then
      export LD_LIBRARY_PATH=${AMBERHOME}/lib
    else
      export LD_LIBRARY_PATH=${AMBERHOME}/lib:${LD_LIBRARY_PATH}
    fi
  fi
fi
