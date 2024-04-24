# This folder contains scripts for running enzyme optimizations using [nlopt](https://nlopt.readthedocs.io/en/latest/)
### Important!!!
Both **Einput7.txt** and **ProteinContentCal.txt** are needed to run the optimization. Currently, they both have the Vmax values of some same enzymes. Please double check them to make sure their Vmax values are identical!

### Prerequisites: 
- miniconda/anaconda - (optional but recommended)
- cmake
- ePhotosynthesis(C++)
- nlopt
- any compilers that support c++11 and above
#### How to check the default and set the c++ standard?
For example, `g++ -x c++  -E -dM -< /dev/null | grep __cplusplus` gives the default standard. 
To specify a version, use the `-std=c++xx` flag.
### File description
- ***.txt**: input data required by the ePhotosynthesis or the optimizations
- **compile_test.sh**: the compile command I use on BioCluster
- **run_ephoto_opt_AUGLAG.cpp**: optimization with AUGLAG methods (since we need a constraint)
- **run_ephoto_opt_Global.cpp**: optimization with global methods that directly support constraints
