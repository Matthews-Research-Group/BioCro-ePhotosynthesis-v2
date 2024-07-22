#include <iostream>
#include <boost/algorithm/string_regex.hpp>
#include <boost/regex.hpp>
#include <vector>

#include "modules/trDynaPS.hpp"
#include "modules/CM.hpp"
#include "modules/EPS.hpp"
#include "modules/PR.hpp"
#include "drivers/drivers.hpp"
#include "Variables.hpp"
using namespace ePhotosynthesis;
using namespace ePhotosynthesis::drivers;
using namespace ePhotosynthesis::modules;

const boost::regex token("\\s+");
struct EPS_inputs
{
    double alpha1=1.0;
    double alpha2=1.0;
};

void EPS_run(double begintime, double stoptime, double stepsize, double abstol, double reltol, double Tp,double PAR, double Ci, int maxSubSteps,const EPS_inputs& EPSinput)
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

       theVars->EnzymeAct.at("V1") *= EPSinput.alpha1;
       theVars->EnzymeAct.at("V2") *= EPSinput.alpha2;
       theVars->EnzymeAct.at("V3") *= EPSinput.alpha2;
       theVars->EnzymeAct.at("V5") *= EPSinput.alpha2;
       theVars->EnzymeAct.at("V6") *= EPSinput.alpha2;
       theVars->EnzymeAct.at("V7") *= EPSinput.alpha2;
       theVars->EnzymeAct.at("V8") *= EPSinput.alpha2;
       theVars->EnzymeAct.at("V9") *= EPSinput.alpha2;
       theVars->EnzymeAct.at("V10") *= EPSinput.alpha2;
       theVars->EnzymeAct.at("V13") *= EPSinput.alpha2;
       theVars->EnzymeAct.at("V23") *= EPSinput.alpha2;

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
       theVars->sensitivity_sf = 1.0;

       PR::setRUBISCOTOTAL(3);

       Driver *maindriver;
       maindriver = new EPSDriver(theVars, begintime, stepsize, stoptime, maxSubSteps, abstol, reltol, 1, 1, Tp,true);
       std::vector<double> ResultRate = maindriver->run();
//       std::cout<<"assim,carboxy,pr are "<<ResultRate[0]<<","<<ResultRate[1]<<","<<ResultRate[2]<<std::endl; 
//       std::cout<<"assim is "<<ResultRate[0]<<std::endl; 

       //write to file output.data, Append
       std::ofstream outfile("output.data", std::ios::out | std::ios::app);
       outfile << PAR<<","<<Tp<<","<<Ci<<","<<ResultRate[0] << std::endl;
       outfile.close();
        if (theVars != nullptr) {
            maindriver->inputVars= nullptr;
            delete theVars;
        }
       delete maindriver;
}

int main(int argc, char* argv[])
{
   EPS_inputs my_inputs;
   double stoptime=5000.0, begintime=0.0, stepsize=0.5;
//   double abstol=1e-5, reltol=1e-4;
   double abstol=1e-6, reltol=1e-6;
   int maxSubSteps=2500;
   int i,curve_option; 
   double Ci;
   //passing in command line arguments
   //start with argv[1] since argv[0] is the program exe
   for (int i = 1; i < argc; i++) { /* We will iterate over argv[] to get the parameters stored inside.
                               * Note that we're starting on 1 because we don't need to know the 
                               * path of the program, which is stored in argv[0] */
     if (i==1) {
       my_inputs.alpha1 = atof(argv[i]);
     } else if(i==2) {
       my_inputs.alpha2 = atof(argv[i]);
     } else{
       curve_option = atoi(argv[i]);
     }
   }
   //std::cout<<my_inputs.alpha1<<","<<my_inputs.alpha2<<std::endl;
   std::remove("output.data");  //remove the old output file.
   if(curve_option==1){//A-Ci
     double PAR= 2000.0;
     double Tp = 25.0;
     std::vector<double> Cis = {100, 150, 200, 250, 300, 400, 500, 600, 800, 1200};
     for (i=0;i < Cis.size();i++) {
      double Ci = Cis[i];
      EPS_run(begintime, stoptime, stepsize, abstol, reltol, Tp, PAR, Ci, maxSubSteps,my_inputs);
     }
   }else if(curve_option==2){//A-Q
     std::vector<double> PARs= {100,300,500,700,900,1100,1300,1500};
     double Tp = 25.0;
     double Ci = 400.0; 
     for (i=0;i < PARs.size();i++) {
      double PAR = PARs[i];
      EPS_run(begintime, stoptime, stepsize, abstol, reltol, Tp, PAR, Ci, maxSubSteps,my_inputs);
     }
   }else{//A-T
     double PAR= 400.0;
     std::vector<double> Tps = {0,5,10,15,20,25,30,35,40};
     double Ci = 400.0; 
     for (i=0;i < Tps.size();i++) {
      double Tp = Tps[i];
      EPS_run(begintime, stoptime, stepsize, abstol, reltol, Tp, PAR, Ci, maxSubSteps,my_inputs);
     }
   }
   return (0);
}
