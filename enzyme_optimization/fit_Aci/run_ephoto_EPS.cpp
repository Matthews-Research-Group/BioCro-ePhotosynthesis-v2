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

void EPS_run(double begintime, double stoptime, double stepsize, double abstol, double reltol, double Tp,double PAR, double Ci, int maxSubSteps)
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

//       std::cout << "EnzymeAct size = " << theVars->EnzymeAct.size() << '\n';

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
       std::cout<<"assim is "<<ResultRate[0]<<std::endl; 

       //write to file output.data, Append
       std::ofstream outfile("output.data",std::ios_base::app);
       outfile << ResultRate[0] << std::endl;
       outfile.close();
        if (theVars != nullptr) {
            maindriver->inputVars= nullptr;
            delete theVars;
        }
       delete maindriver;
}

int main()
{
    double stoptime=5000.0, begintime=0.0, stepsize=0.5;
//    double abstol=1e-5, reltol=1e-4;
    double abstol=9.9e-6, reltol=1e-5;
    double PAR= 1500.0;
    double Tp = 25.0;
//    double Cis[20] = {10,60,110,160,210,260,310,360,410,460,510,560,610,660,710,760,810,860,910,960};
//    std::vector<double> Cis = {110,210,260,310,360,410,460,510,560,610,660};
    std::vector<double> Cis = {460,510};
    int maxSubSteps=2500;
    int i; 
    double Ci;
   for (i=0;i < Cis.size();i++) {
    std::cout<<"i is "<<i<<std::endl; 
    Ci = Cis[i];
    EPS_run(begintime, stoptime, stepsize, abstol, reltol, Tp, PAR, Ci, maxSubSteps);
   }

    return (0);
}
