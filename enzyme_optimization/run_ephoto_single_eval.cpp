#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <stdlib.h> 

#include <nlopt.hpp>

#include <boost/algorithm/string_regex.hpp>
#include <boost/regex.hpp>

#include "modules/trDynaPS.hpp"
#include "modules/CM.hpp"
#include "modules/EPS.hpp"
#include "modules/PR.hpp"
#include "drivers/drivers.hpp"
#include "Variables.hpp"
#include "globals.hpp"

#define number_of_parameters 25
//remove PR ones, we then have 18 to optimize
//#define number_of_parameters 18 
using namespace ePhotosynthesis;
using namespace ePhotosynthesis::drivers;
using namespace ePhotosynthesis::modules;

const boost::regex token("\\s+");

struct EPS_inputs
{
    double stoptime=5000.0;
    double begintime=0.0;
    double stepsize=0.5;
    double abstol=1e-6;
    double reltol=1e-6;
    double PAR= 1500.0;
    double Tp = 25.0;
    double Ci = 300; 
    int maxSubSteps=2500;
    double totalE=0.0;
    std::vector<std::string> ename{}; 
    std::vector<double> vmax{}; 
    std::vector<double> kcat{}; 
    std::vector<double> mweight{}; 
};

void print_map(const std::string& comment, const std::map<std::string, int>& m)
{
  std::cout << comment;

// C++11 alternative:
  for (const auto& n : m)
  {
      std::cout << n.first << " = " << n.second << "; ";
  }
  std::cout << '\n';
}

void print(std::vector <double> const &a) {
   for(int i=0; i < a.size(); i++)
   std::cout << a.at(i) << ' ';
}

double myvconstraint(const std::vector<double> &parameters, std::vector<double> &grad, void *data)
{
    EPS_inputs *d = reinterpret_cast<EPS_inputs*>(data);
    std::vector<double> vmax = d->vmax;
    std::vector<double> kcat = d->kcat;
    std::vector<double> mweight = d->mweight;
    if (!grad.empty()) {
       std::cout<<"grad should be empty!YH";
       exit(3);
    }
    double constrain_E = 0.0;
    for (int it = 0; it < parameters.size(); ++it) {
      constrain_E += parameters[it]*vmax[it]/kcat[it]*mweight[it]/1000.0;
    }
    //take the difference between this totalE and the default total 
    //std::cout<<"constrain_E "<<constrain_E;
    return abs(constrain_E - (d->totalE)); 
} 

double EPS_run(const double& begintime, 
             const double& stoptime, 
             const double& stepsize, 
             const double& abstol, 
             const double& reltol, 
             const double& Tp,
             const double& PAR,
             const double& Ci, 
             const int&    maxSubSteps,
             const double& totalE, //not used in this function
             const std::vector<std::string> &ename,
             const std::vector<double> &vmax, //not used in this function
             const std::vector<double> &kcat, //not used in this function
             const std::vector<double> &mweight, //not used in this function
             const std::vector<double> &parameter_sf) //parameters to be optimized
{
//since we check if metabolites reach steady-state, we have to use this to calculate the metabolites.
       bool record = true;
       bool saveMetabolite = true;//this outputs all metabolites' time series
       std::string evn="InputEvn.txt";
       std::string atpcost="InputATPCost.txt";
       std::string enzymeFile="Einput7.txt";
       std::map<std::string, std::string> inputs;

       readFile(evn, inputs);
       readFile(atpcost, inputs);
       Variables *theVars = new Variables();
       //YH:readFile actually changes units by dividing 30!
       readFile(enzymeFile, theVars->EnzymeAct); 
       int k=1;
       for(auto& x:ename)
       {
          theVars->EnzymeAct.at(x) *= parameter_sf[k-1];
//        std::cout<<x<<","<<theVars->EnzymeAct.at(x)<<std::endl;
          k++;
       }
       //remove Ca and Light inputs from the text file. Instead, they are function inputs
       //theVars->TestCa = static_cast<double>(stof(inputs.at("CO2"), nullptr));
       //theVars->TestLi = static_cast<double>(stof(inputs.at("PAR"), nullptr));
       theVars->CO2_in = Ci; 
       theVars->TestLi = PAR; 
       if (stoi(inputs.at("SucPath"), nullptr) > 0)  CM::setTestSucPath(true);
       theVars->TestATPCost = stoi(inputs.at("ATPCost"), nullptr);
       theVars->record = record;
       theVars->saveMetabolite = saveMetabolite;
       theVars->useC3 = true;     //for EPSDriver
       theVars->RUBISCOMETHOD = 2;
       PR::setRUBISCOTOTAL(3);

       Driver *maindriver;
       maindriver = new EPSDriver(theVars, begintime, stepsize, stoptime, maxSubSteps, abstol, reltol, 1, 1, Tp);
       std::vector<double> ResultRate = maindriver->run();

       std::cout<<" "<<ResultRate[0]<<std::endl; 

       //write to file output.data, Append
//       std::ofstream outfile("output.data",std::ios_base::app);
//      outfile << ResultRate[0] << std::endl;
//     outfile.close();
        if (theVars != nullptr) {
            maindriver->inputVars= nullptr;
            delete theVars;
        }
       delete maindriver;

       return ResultRate[0];
}

std::vector<double> readValuesFromFile(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Error opening file" << std::endl;
        return {};
    }

    std::vector<double> values;
    std::string line;

    // Read the file line by line
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        double value;
        // Read each value separated by space and add to the vector
        while (iss >> value) {
            values.push_back(value);
        }
    }

    inputFile.close();
    return values;
}

int main(int argc, char* argv[])
{
//read in txt inputs
   EPS_inputs my_inputs;
   //passing in command line arguments
   //start with argv[1] since argv[0] is the program exe
   for (int i = 1; i < argc; i++) { /* We will iterate over argv[] to get the parameters stored inside.
                               * Note that we're starting on 1 because we don't need to know the 
                               * path of the program, which is stored in argv[0] */
// Check that we haven't finished parsing already
// And, the argument is not a number, which is at even number position
       if (i + 1 != argc && i%2 != 0) 
       { 
          if (strcmp(argv[i], "-r")==0) {
            my_inputs.PAR = atof(argv[i + 1]);
          } else if (strcmp(argv[i], "-t")==0) {
            my_inputs.Tp = atof(argv[i + 1]);
          } else if (strcmp(argv[i], "-c")==0) {
            my_inputs.Ci = atof(argv[i + 1]);
          } else {
            std::cout << "Not enough or invalid arguments, please try again.\n";
            exit(0);
          }
      }
   }
   std::cout<< "PAR Tp Ci are "<<my_inputs.PAR<<" "<<my_inputs.Tp<<" "<<my_inputs.Ci<<std::endl;
   std::ifstream infile;
   infile.open("ProteinContentCal.txt");// 
   std::string line;
   while (std::getline(infile, line))
   {
       std::stringstream ss(line);
       std::string name;
       double vmax, kcat, mweight;
       if (ss >> name>> vmax >> kcat >> mweight)
       {
          if(name=="V16") continue; // remove V16: ATP synthase
//          if(name=="V112"||name=="V113"||name=="V121"||name=="V122"||name=="V123"||name=="V124"||name=="V131") continue; // remove photorespiration from optimization 
          my_inputs.ename.push_back(name);
          my_inputs.vmax.push_back(vmax);
          my_inputs.kcat.push_back(kcat);
          my_inputs.mweight.push_back(mweight);
          my_inputs.totalE += vmax/kcat * mweight/1000.0; 
       }
   }
   infile.close();// 
   std::cout<<"default TotalE is "<<my_inputs.totalE<<"\n";

    EPS_inputs *f_data_ptr =  &my_inputs; 
//define scaling factors
    const std::string filename = "optimized_enzyme_scaling_factors.txt";
    std::vector<double> parameter_sf = readValuesFromFile(filename);
//run EPS Driver 
    EPS_run(f_data_ptr->begintime,
                   f_data_ptr->stoptime, 
                   f_data_ptr->stepsize, 
                   f_data_ptr->abstol,
                   f_data_ptr->reltol, 
                   f_data_ptr->Tp,
                   f_data_ptr->PAR,
                   f_data_ptr->Ci,
                   f_data_ptr->maxSubSteps,
                   f_data_ptr->totalE,
                   f_data_ptr->ename,
                   f_data_ptr->vmax,
                   f_data_ptr->kcat,
                   f_data_ptr->mweight,
                   parameter_sf);
    return (0);
}
