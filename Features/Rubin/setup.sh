#! /bin/bash
#
# file: setup.sh
#
# This bash script performs any setup necessary in order to test your
# entry.  It is run only once, before running any other code belonging
# to your entry.

set -e
set -o pipefail

# Remove (or set it to 0) if you are not using Matlab
NEED_MATLAB=1

# Example: compile a C module for Matlab
cd matlab
mex viterbi_Springer.c
cd ..

cd python_speech_features
python setup.py install --user
cd ..

# Install matlab-engine for python
MATLAB_BIN=$(dirname $(which matlab))
cd $MATLAB_BIN/../extern/engines/python && python setup.py build --build-base=/tmp/mlpy install --user

