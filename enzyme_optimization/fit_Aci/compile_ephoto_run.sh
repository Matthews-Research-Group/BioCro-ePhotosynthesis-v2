#!/bin/bash

# Enable aliases in the script
alias myclangcpp='clang++ -std=c++20'
# Get the directory of the script
SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname "$SCRIPT")
# Now go back twice
MAIN_DIR=$(dirname "$(dirname "$SCRIPTPATH")")
ephoto_path="$MAIN_DIR/models/ePhotosynthesis_C"
# Get the directory of the conda env named "ephoto"
conda_path=$(conda info --envs | grep ephoto | awk '{print $NF}')

myclangcpp -o myephoto.exe -I$conda_path/include -I$ephoto_path/include run_ephoto_EPS.cpp -L$ephoto_path/build -lEPhotosynthesis -Wl,-rpath,$ephoto_path/build 
