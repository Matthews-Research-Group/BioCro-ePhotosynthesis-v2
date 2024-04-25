#!/bin/bash

#no longer need to specify boost path. On biocluster,it's a loaded module
#boost_path=/home/n-z/yufenghe/.conda/envs/ephoto/include
sundial_path=/home/n-z/yufenghe/.conda/envs/ephoto/
ephoto_path=/home/n-z/yufenghe/GitHub/BioCro-ePhotosynthesis-v2//models/ePhotosynthesis_C/

g++ -o enzyme_opt_global.exe run_ephoto_opt_Global.cpp -O3 -g -fopenmp -Wall -ggdb3 -I$ephoto_path/include  -I$sundial_path/include/ -lm -L/home/n-z/yufenghe/.conda/envs/ephoto/lib -lnlopt -lboost_regex -L$ephoto_path/build -lEPhotosynthesis -Wl,-rpath=$ephoto_path/build

#g++ -o enzyme_opt_LagMLSL_LDS.exe run_ephoto_opt_AUGLAG.cpp -O3 -g -fopenmp -Wall -ggdb3 -I$ephoto_path/include  -I$sundial_path/include/ -lm -L/home/n-z/yufenghe/.conda/envs/ephoto/lib -lnlopt -L$ephoto_path/build -lEPhotosynthesis -Wl,-rpath=$ephoto_path/build

#g++ -o enzyme_opt_local_BOBYQA.exe run_ephoto_opt_Local.cpp -O3 -g -fopenmp -Wall -ggdb3 -I$ephoto_path/include  -I$sundial_path/include/ -lm -L/home/n-z/yufenghe/.conda/envs/ephoto/lib -lnlopt -L$ephoto_path/build -lEPhotosynthesis -Wl,-rpath=$ephoto_path/build

#g++ -o test1.exe nlopt_ex1.cpp -O3 -g -fopenmp -Wall -ggdb3 -I$ephoto_path/include  -I$sundial_path/include/ -lm -L/home/n-z/yufenghe/.conda/envs/ephoto/lib -lnlopt -L$ephoto_path/build -lEPhotosynthesis -Wl,-rpath=$ephoto_path/build

#g++ -o test2.exe nlopt_ex1_AUGLAG.cpp -O3 -g -fopenmp -Wall -ggdb3 -I$ephoto_path/include  -I$sundial_path/include/ -lm -L/home/n-z/yufenghe/.conda/envs/ephoto/lib -lnlopt -L$ephoto_path/build -lEPhotosynthesis -Wl,-rpath=$ephoto_path/build

#g++ test.cpp -O3 -g -fopenmp -Wall -ggdb3 -I$ephoto_path/include  -I$sundial_path/include/ -lm -L/home/n-z/yufenghe/.conda/envs/ephoto/lib -lnlopt -L$ephoto_path/build -lEPhotosynthesis -Wl,-rpath=$ephoto_path/build
