#!/bin/bash

boost_path=/home/n-z/yufenghe/.conda/envs/ephoto/include
sundial_path=/home/n-z/yufenghe/.conda/envs/ephoto/
ephoto_path=/home/n-z/yufenghe/GitHub/BioCro-ePhotosynthesis/models/ePhotosynthesis_C/
g++ -o enzyme_opt.exe run_ephoto_opt.cpp -O3 -g -fopenmp -Wall -ggdb3 -I$boost_path -I$ephoto_path/include  -I$sundial_path/include/ -lm -L/home/n-z/yufenghe/.conda/envs/ephoto/lib -lnlopt -L$ephoto_path/build -lEPhotosynthesis -Wl,-rpath=$ephoto_path/build

