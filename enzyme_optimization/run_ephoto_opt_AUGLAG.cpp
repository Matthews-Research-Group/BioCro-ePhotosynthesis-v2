#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>

#include <nlopt.hpp>

#include <boost/algorithm/string_regex.hpp>
#include <boost/regex.hpp>

#include "modules/trDynaPS.hpp"
#include "modules/CM.hpp"
#include "modules/EPS.hpp"
#include "modules/PR.hpp"
#include "drivers/drivers.hpp"
#include "Variables.hpp"

#define number_of_parameters 25
using namespace ePhotosynthesis;
using namespace ePhotosynthesis::drivers;
using namespace ePhotosynthesis::modules;

const boost::regex token("\\s+");

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
       std::cout<<"WARNING: grad is not empty!YH";
    }
    double constrain_E = 0.0;
    for (int it = 0; it < parameters.size(); ++it) {
      constrain_E += parameters[it]*vmax[it]/kcat[it]*mweight[it]/1000.0;
    }
    //take the difference between this totalE and the default total 
    //std::cout<<"constrain_E "<<constrain_E;
    return constrain_E - (d->totalE); 
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
             const double& totalE,
             const std::vector<std::string> &ename,
             const std::vector<double> &vmax,
             const std::vector<double> &kcat,
             const std::vector<double> &mweight,
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

int count=0;
double myvfunc(const std:: vector<double> &parameter_sf, std::vector<double> &grad, void *f_data)
{
    EPS_inputs *f_data_ptr = (EPS_inputs *) f_data; 
    ++count;
//    std::cout << "this is the evalution # "<<count<<std::endl;  
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
                   f_data_ptr->totalE,
                   f_data_ptr->ename,
                   f_data_ptr->vmax,
                   f_data_ptr->kcat,
                   f_data_ptr->mweight,
                   parameter_sf);
}

int main()
{
//read in txt inputs
   EPS_inputs my_inputs;
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
          my_inputs.ename.push_back(name);
          my_inputs.vmax.push_back(vmax);
          my_inputs.kcat.push_back(kcat);
          my_inputs.mweight.push_back(mweight);
          my_inputs.totalE += vmax/kcat * mweight/1000.0; 
       }
   }
   infile.close();// 
   std::cout<<"default TotalE is"<<my_inputs.totalE<<"\n";
//define optimizer    
    nlopt::opt opt(nlopt::AUGLAG, number_of_parameters);
    nlopt::opt local_opt(nlopt::GN_DIRECT_L, number_of_parameters);
    local_opt.set_xtol_rel(1e-3);
    opt.set_local_optimizer(local_opt);
//lower and upper bounds
    std::vector<double> lb(number_of_parameters,0.1);
    std::vector<double> ub(number_of_parameters,10.0);
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
//objective function
    opt.set_max_objective(myvfunc, &my_inputs);
//constraints
    opt.add_inequality_constraint(myvconstraint, &my_inputs, 1e-8);
//tolerance for stopping
    opt.set_xtol_rel(1e-3);
//initial guess
    std::vector<double> x(number_of_parameters,1.0);
    double maxf;
//run optimization    
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