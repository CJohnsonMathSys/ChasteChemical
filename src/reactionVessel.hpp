# pragma once

#include <string>
#include <tuple>
#include <vector>

class reactionVessel
{
private:
  // file names
  std::string inputReactionFilename;
  std::string outputReactionFilename;
  std::string inputMixtureFilename;
  std::string ouputMixtureFilename;

  // chemistry
  std::vector<std::string> reactingSpeciesNames;
  std::vector<std::vector<int>> reactingStoichMatrix;
  std::vector<std::vector<double>> kineticValue;
  std::vector<std::vector<std::string>> reactionInfo;

  // composition
  std::vector<std::string> mixtureSpeciesNames;
  std::vector<double> mixtureInitialConcentrations;
  std::vector<double> reactionRates;
  std::vector<double> mixtureConcentrations;
public:
  //constructor
  reactionVessel();
  reactionVessel(std::string,std::string,std::string,std::string); // pass file names ot intitialise the vessel
  //destructor
  ~reactionVessel();

  // file names
  // setter methods
  void setReactionFile(std::string); // only input file
  void setReactionFile(std::string, std::string); // both input and output reacion file
  void setMixtureFile(std::string); // only input mixture file
  void setMixtureFile(std::string, std::string); // both input and output mixture files
  // getter methods
  std::string getReactionFilename(int); // return reaction file names; int=1 for input, int=0 for output, default input
  std::string getMixtureFilename(int); // return mixture file names; int=1 for input, int=0 for output, default input

  // chemistry
  // setter methods
  void setReactingSpeciesNames(std::vector<std::string>); // method to set the species name array
  void setReactingStoichiometry(std::vector<std::vector<int>>); // method to set the stoichiometric matrix
  void setReactionKinetics(std::vector<std::vector<double>>); // set kineitc values vector
  void setKineticLabels(std::vector<std::vector<std::string>>); // set labels for the kineitc vector
  // getter methods
  std::vector<std::string> getReactingSpeciesNames();
  std::vector<std::vector<int>> getReactingStoichiometry();
  std::vector<std::vector<double>> getReactionKinetics(); // return the reaciton kinetic values as a vector
  std::vector<std::vector<std::string>> getKineticLabels(); // return the label vector associated with the reaciton kinetics

  // chemistry functions
  // function to add a reaction to the system
  void addReactionFromFile(std::string);
  void addReactionFromString(std::string);
  // removing reactions non trivial, how to know whether reaction already exists?


  // composition
  // setter
  void setMixtureSpeciesNames(std::vector<std::string>);
  void setMixtureInitialConcentrations(std::vector<double>);
  void setReactionRates(std::vector<double>);
  void setMixtureConcentrations(std::vector<double>);
  // getter
  std::vector<std::string> getMixtureSpeciesNames();
  std::vector<double> getMixtureInitialConcentrations();
  std::vector<double> getReactionRates();
  std::vector<double> getMixtureConcentrations();

  std::vector<double> calculateReactionQuotient(std::vector<bool>);
  std::vector<double> calculateReactionRates(double);
  std::vector<double> calculateConcentrationChange();
  // write methods

  void printReactionVesselToFile(std::string);
  void printReactionSystemToFile(std::string);
  void printVesselMixtureToFile(std::string);

};
