#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>

#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "PetscException.hpp"

#include "reactionVessel.hpp"
#include "reactionInputProcessing.hpp"

int main(int argc, char *argv[])
{
    ExecutableSupport::StandardStartup(&argc, &argv);
    int exit_code = ExecutableSupport::EXIT_OK;
    try
    {


      std::string mixtureFilename = readMixtureComposition();
      std::cout<<"Mixture filename: "<<mixtureFilename<<std::endl;

      //std::cout<<"Reaction filename: "<<reactionFilename<<std::endl;

      //reactionVessel cell(reactionFilename);
      reactionVessel cell;
      std::cout<<"Print the stoichiometric matrix for the system: "<<std::endl;
      printVec(cell.getReactingSpeciesNames());
      printVecVec(cell.getReactingStoichiometry());
      std::cout<<"Print the mixture composition for the system: "<<std::endl;
      printVec(cell.getMixtureSpeciesNames());
      printVec(cell.getmixtureInitialConcentrations());
      //std::string newReactionString = "4C + 2Y <-> ATP | kf = 3.4 kr = 2.4";
      std::string newReactionString = "A <-> 2Z + 4C | kf = 1.4";
      std::cout<<"Add new reaction string: "+newReactionString<<" "<<std::endl;

      cell.addReactionFromString(newReactionString);
      std::cout<<"Print the stoichiometric matrix for the system: "<<std::endl;
      printVec(cell.getReactingSpeciesNames());
      printVecVec(cell.getReactingStoichiometry());
      std::cout<<"Print the mixture composition for the system: "<<std::endl;
      printVec(cell.getMixtureSpeciesNames());
      printVec(cell.getmixtureInitialConcentrations());

      std::string newReactionFile = "/home/chaste/projects/reactionDiffusion/apps/src/Data/reversible.txt";
      std::cout<<"Add reaction from file: reversible.txt"<<" "<<std::endl;
      cell.addReactionFromFile(newReactionFile);
      std::cout<<"Print the stoichiometric matrix for the system: "<<std::endl;
      printVec(cell.getReactingSpeciesNames());
      printVecVec(cell.getReactingStoichiometry());
      std::cout<<"Print the mixture composition for the system: "<<std::endl;
      printVec(cell.getMixtureSpeciesNames());
      printVec(cell.getmixtureInitialConcentrations());
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }
    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
