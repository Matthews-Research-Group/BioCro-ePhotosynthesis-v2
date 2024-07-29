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
    double Q10_1  = 1.93;
    double Q10_2  = 2.0;
    double Q10_3  = 2.0;
    double Q10_5  = 2.0;
    double Q10_6  = 2.0;
    double Q10_7  = 2.0;
    double Q10_8  = 2.0;
    double Q10_9  = 2.0;
    double Q10_10 = 2.0;
    double Q10_13 = 2.0;
    double Q10_23 = 2.0;
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

       theVars->Q10_1  = EPSinput.Q10_1; 
       theVars->Q10_2  = EPSinput.Q10_2; 
       theVars->Q10_3  = EPSinput.Q10_3; 
       theVars->Q10_5  = EPSinput.Q10_5; 
       theVars->Q10_6  = EPSinput.Q10_6; 
       theVars->Q10_7  = EPSinput.Q10_7; 
       theVars->Q10_8  = EPSinput.Q10_8; 
       theVars->Q10_9  = EPSinput.Q10_9; 
       theVars->Q10_10 = EPSinput.Q10_10;
       theVars->Q10_13 = EPSinput.Q10_13;
       theVars->Q10_23 = EPSinput.Q10_23;
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

       std::cout << ResultRate[0] << std::endl;
       //write to file output.data, Append

       //std::ofstream outfile("output.data");
       //outfile << ResultRate[0] << std::endl;
       //outfile.close();

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
   int i; 
   double PAR,Tp,Ci;
   // Create an array of pointers to the struct members
    std::vector<double*> Q10_pointers = {
        &my_inputs.Q10_1,  &my_inputs.Q10_2, &my_inputs.Q10_3,
        &my_inputs.Q10_5,  &my_inputs.Q10_6, &my_inputs.Q10_7,
        &my_inputs.Q10_8,  &my_inputs.Q10_9, &my_inputs.Q10_10,
        &my_inputs.Q10_13, &my_inputs.Q10_23
    };
   //passing in command line arguments
   //start with argv[1] since argv[0] is the program exe
   if (argc==6){
     for (int i = 1; i < argc; i++) { /* We will iterate over argv[] to get the parameters stored inside.
                                 * Note that we're starting on 1 because we don't need to know the 
                                 * path of the program, which is stored in argv[0] */
       if (i==1) {
         my_inputs.alpha1 = atof(argv[i]);
       } else if(i==2) {
         my_inputs.alpha2 = atof(argv[i]);
       } else if(i==3) {
         PAR = atof(argv[i]);
       } else if(i==4) {
         Tp  = atof(argv[i]);
       } else if(i==5) {
         Ci  = atof(argv[i]);
       }
     }
   }else if (argc==17){
     for (int i = 1; i < argc; i++) { /* We will iterate over argv[] to get the parameters stored inside.
                                 * Note that we're starting on 1 because we don't need to know the 
                                 * path of the program, which is stored in argv[0] */
       if (i==1) {
         my_inputs.alpha1 = atof(argv[i]);
       } else if(i==2) {
         my_inputs.alpha2 = atof(argv[i]);
       } else if(i==3) {
         PAR = atof(argv[i]);
       } else if(i==4) {
         Tp  = atof(argv[i]);
       } else if(i==5) {
         Ci  = atof(argv[i]);
       }else{
        *(Q10_pointers[i]) = atof(argv[i]);
       }
     }
   }else{
        // If an error occurs
     std::cerr << "Error: incorrect number of input arguments" << std::endl;
     std::exit(EXIT_FAILURE); // or std::exit(1);
   }
   EPS_run(begintime, stoptime, stepsize, abstol, reltol, Tp, PAR, Ci, maxSubSteps,my_inputs);
   return (0);
}
