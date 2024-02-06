#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

#include <nlopt.hpp>

#include <boost/algorithm/string_regex.hpp>
#include <boost/regex.hpp>

#include "modules/trDynaPS.hpp"
#include "modules/CM.hpp"
#include "modules/EPS.hpp"
#include "modules/PR.hpp"
#include "drivers/drivers.hpp"
#include "Variables.hpp"

#define number_of_parameters 26 
using namespace ePhotosynthesis;
using namespace ePhotosynthesis::drivers;
using namespace ePhotosynthesis::modules;

const boost::regex token("\\s+");

typedef struct {
    double a, b;
} my_constraint_data;

struct EPS_inputs
{
    double stoptime=5000.0;
    double begintime=0.0;
    double stepsize=0.5;
    double abstol=9.9e-6;
    double reltol=1e-5;
    double PAR= 1500.0;
    double Tp = 25.0;
    double Ci = 300; 
    int maxSubSteps=2500;
};

void print_map(const std::string& comment, const std::map<std::string, int>& m)
{
  std::cout << comment;
//    // Iterate using C++17 facilities
//    for (const auto& [key, value] : m)
//        std::cout << '[' << key << "] = " << value << "; ";

// C++11 alternative:
  for (const auto& n : m)
  {
      std::cout << n.first << " = " << n.second << "; ";
  }

// C++98 alternative:
//  for (std::map<std::string, int>::const_iterator it = m.begin(); it != m.end(); ++it)
//      std::cout << it->first << " = " << it->second << "; ";

    std::cout << '\n';
}

void print(std::vector <double> const &a) {
   for(int i=0; i < a.size(); i++)
   std::cout << a.at(i) << ' ';
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
             const std::vector<double> &parameter_sf)
{
       bool record = false;
       std::string evn="InputEvn.txt";
       std::string atpcost="InputATPCost.txt";
       std::string enzymeFile="Einput7.txt";
       std::map<std::string, std::string> inputs;

       readFile(evn, inputs);
       readFile(atpcost, inputs);
       Variables *theVars = new Variables();
       readFile(enzymeFile, theVars->EnzymeAct);

//       const int number_of_parameters = 26;//these enzymes will be scaled
       int k=1;
       for(auto& x:theVars->EnzymeAct)
       {
         if (k<=number_of_parameters)
         {
           x.second *= parameter_sf[k-1];
         }
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
       theVars->useC3 = true;     //for EPSDriver
       theVars->RUBISCOMETHOD = 2;
       PR::setRUBISCOTOTAL(3);

       Driver *maindriver;
       maindriver = new EPSDriver(theVars, begintime, stepsize, stoptime, maxSubSteps, abstol, reltol, 1, 1, Tp);
       std::vector<double> ResultRate = maindriver->run();

       std::cout<<"with An = "<<ResultRate[0]<<std::endl; 

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

int count=0;
double myvfunc(const std:: vector<double> &parameter_sf, std::vector<double> &grad, void *f_data)
{
    EPS_inputs *f_data_ptr = (EPS_inputs *) f_data; 
    ++count;
    std::cout << "this is the evalution # "<<count<<std::endl;  
    std::cout<< "the scaling factors are ";
    print(parameter_sf);
    return EPS_run(f_data_ptr->begintime,
                   f_data_ptr->stoptime, 
                   f_data_ptr->stepsize, 
                   f_data_ptr->abstol,
                   f_data_ptr->reltol, 
                   f_data_ptr->Tp,
                   f_data_ptr->PAR,
                   f_data_ptr->Ci,
                   f_data_ptr->maxSubSteps,
                   parameter_sf);
}

int main()
{
    EPS_inputs my_inputs;
    
    nlopt::opt opt(nlopt::GN_DIRECT, number_of_parameters);

//lower and upper bounds
    std::vector<double> lb(number_of_parameters,0.5);
    std::vector<double> ub(number_of_parameters,2.0);

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_max_objective(myvfunc, &my_inputs);

    opt.set_xtol_rel(1e-4);

//initial guess
    std::vector<double> x(number_of_parameters,1.0);
    double maxf;
    
    try{
        nlopt::result result = opt.optimize(x, maxf);
        std::cout << "found max at f("; 
        print(x);
        std::cout << "="; 
        std::cout << std::setprecision(10) << maxf << std::endl;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
  
    return (0);
}
