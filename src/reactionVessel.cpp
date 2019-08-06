#include "reactionInputProcessing.hpp"
#include "reactionVessel.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <tuple>
#include <cmath>
#include <stdlib.h>
#include <algorithm>

//constructor
reactionVessel::reactionVessel(){

  // file names
  inputReactionFilename = "/home/chaste/projects/reactionDiffusion/apps/src/Data/brusselator.txt";
  outputReactionFilename = "/home/chaste/projects/reactionDiffusion/apps/src/Data/defaultreactionOut.txt";
  inputMixtureFilename = "/home/chaste/projects/reactionDiffusion/apps/src/Data/mixture.csv";
  ouputMixtureFilename = "/home/chaste/projects/reactionDiffusion/apps/src/Data/defaultMixutreOut.csv";

  std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> >> system = parseReactionsFromFile(inputReactionFilename);
  std::tuple<std::vector<std::string>,std::vector<std::vector<int>>> stoichTuple = composeStoichMatrix(system);
  std::tuple<std::vector<std::string>, std::vector<double>> compositionTuple;
  compositionTuple = composeCompositionVectorFromFile(inputMixtureFilename,std::get<0>(stoichTuple));

  // chemistry
  reactingSpeciesNames = std::get<0>(stoichTuple);
  reactingStoichMatrix = std::get<1>(stoichTuple);

  std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<double>>> kineticTuple = composeKinetics(system);
  reactionInfo = std::get<0>(kineticTuple);
  kineticValue = std::get<1>(kineticTuple);

  mixtureSpeciesNames = std::get<0>(compositionTuple);
  mixtureInitialConcentrations = std::get<1>(compositionTuple);
  // set concentratioons to initial values
  mixtureConcentrations = mixtureInitialConcentrations;
  std::vector<double> nullRates(reactionInfo.size(),0.0);
  reactionRates = nullRates;

};

reactionVessel::reactionVessel(std::string iRF="brusselator.txt",std::string oRF="defoRF",std::string iMF="mixture.csv",std::string oMF="defoMF"){

  // file names
  inputReactionFilename = iRF;
  outputReactionFilename = oRF;
  inputMixtureFilename = iMF;
  ouputMixtureFilename = oMF;

  std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> >> system = parseReactionsFromFile(iRF);
  std::tuple<std::vector<std::string>,std::vector<std::vector<int>>> stoichTuple = composeStoichMatrix(system);
  std::tuple<std::vector<std::string>, std::vector<double>> compositionTuple;
  compositionTuple = composeCompositionVectorFromFile(iMF,std::get<0>(stoichTuple));

  // chemistry
  reactingSpeciesNames = std::get<0>(stoichTuple); // vector of species names as strings
  reactingStoichMatrix = std::get<1>(stoichTuple); // vecvec [reaction][species] = stoichvalue


  std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<double>>> kineticTuple = composeKinetics(system);
  reactionInfo = std::get<0>(kineticTuple);
  kineticValue = std::get<1>(kineticTuple);

  mixtureSpeciesNames = std::get<0>(compositionTuple);
  mixtureInitialConcentrations = std::get<1>(compositionTuple);
  // set concentratioons to initial values
  mixtureConcentrations = mixtureInitialConcentrations;
  std::vector<double> nullRates(reactionInfo.size(),0.0);
  reactionRates = nullRates;
};
//destructor
reactionVessel::~reactionVessel(){};

void reactionVessel::setReactionFile(std::string file){
  inputReactionFilename = file;
  return;
};

void reactionVessel::setReactionFile(std::string inFile, std::string outFile){
  inputReactionFilename = inFile;
  outputReactionFilename = outFile;
  return;
};

void reactionVessel::setMixtureFile(std::string file){
  inputMixtureFilename = file;
  return;
};

void reactionVessel::setMixtureFile(std::string inFile, std::string outFile){
  inputMixtureFilename = inFile;
  ouputMixtureFilename = outFile;
  return;
};

std::string reactionVessel::getReactionFilename(int select = 1){
  // return reaction file names; int=1 for input, int=0 for output, default input
  if(select == 0){
    // return output reaction filename
    return outputReactionFilename;
  }else{
    return inputReactionFilename;
  }
};

std::string reactionVessel::getMixtureFilename(int select = 1){
  // return mixture file names; int=1 for input, int=0 for output, default input
  if(select == 0){
    return ouputMixtureFilename;
  }else{
    return inputMixtureFilename;
  }
};

// chemistry
// setter methods
void reactionVessel::setReactingSpeciesNames(std::vector<std::string> names){
  reactingSpeciesNames = names;
  return;
};

void reactionVessel::setReactingStoichiometry(std::vector<std::vector<int>> stoich){
  reactingStoichMatrix = stoich;
  return;
};

void reactionVessel::setReactionKinetics(std::vector<std::vector<double>> values){
  kineticValue = values;
  return;
};

void reactionVessel::setKineticLabels(std::vector<std::vector<std::string>> labels){
  reactionInfo = labels;
  return;
};

// getter methods
std::vector<std::string> reactionVessel::getReactingSpeciesNames(){
  return reactingSpeciesNames;
};

std::vector<std::vector<int>> reactionVessel::getReactingStoichiometry(){
  return reactingStoichMatrix;
};

std::vector<std::vector<double>> reactionVessel::getReactionKinetics(){
  return kineticValue;
};

std::vector<std::vector<std::string>> reactionVessel::getKineticLabels(){
  return reactionInfo;
};

// chemistry functions
void reactionVessel::addReactionFromFile(std::string addFilename){
  // check filename is in correct format
  std::string fileType = ".txt";
  if (addFilename.find(fileType) != std::string::npos){
    // filename string in correct form
  }else{
    // not a .txt file, append
    std::cout<< "Warning: Given file is not a .txt file, appending."<<std::endl;
    addFilename.append(fileType);
  }
  // read and parse new reaction file
  std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> >> newSystem = parseReactionsFromFile(addFilename);


  // get current reactions
  std::vector<std::string> currentRSN = getReactingSpeciesNames();
  std::vector<std::vector<int>> currentRSM = getReactingStoichiometry();
  std::vector<std::vector<double>> currentKV = getReactionKinetics();
  std::vector<std::vector<std::string>> currentRI = getKineticLabels();

  std::vector<std::vector<int>> newRSM;

  std::vector<std::string> tempRSN;
  std::vector<std::vector<int>> tempRSM;
  std::vector<std::vector<double>> tempKV;
  std::vector<std::vector<std::string>> tempRI;

  // reformat the parsed reaction
  std::tuple<std::vector<std::string>,std::vector<std::vector<int>>> stoichTuple = composeStoichMatrix(newSystem);
  // returns speciesArray, stoichMatrix


  tempRSN = std::get<0>(stoichTuple);
  tempRSM = std::get<1>(stoichTuple);

  // combine the species arrays
  //auto iter = newRSN.insert(newRSN.begin(), tempRSN);
  //newRSN.insert(iter, currentRSN);

  std::vector<std::string> newRSN(tempRSN.size() + currentRSN.size());
  std::copy(currentRSN.begin(), currentRSN.end(), newRSN.begin());
  std::copy(tempRSN.begin(), tempRSN.end(), newRSN.begin() + currentRSN.size());
  // remove duplicated species name entries
  sort(newRSN.begin(), newRSN.end());
  auto last = unique(newRSN.begin(), newRSN.end());
  newRSN.erase(last, newRSN.end());


  //unsigned int numberOfReactions = newRSM.size();
  // check if new reaction introduces new species, if so update species vectors and initialise a zero concentration
  if(newRSN.size() == currentRSN.size()){
    // no new species have appeared, append the reaction stoichiometry
    // species ordering should be the same
    // combine the stoich vectors to add new reaction(s)
    newRSM.resize(currentRSM.size() + tempRSM.size());
    std::copy(currentRSM.begin(), currentRSM.end(), newRSM.begin());
    std::copy(tempRSM.begin(), tempRSM.end(), newRSM.begin() + currentRSM.size());


  }else{
    // new species require reforming of the stochiometric matrix
    unsigned int numberOfSpecies = newRSN.size();
    unsigned int numberOfReactions = currentRSM.size() + tempRSM.size();
    //std::vector<std::vector<int>> newRSM(numberOfReactions, std::vector<int>(numberOfSpecies,0));
    newRSM.resize(numberOfReactions);
    std::fill(newRSM.begin(), newRSM.end(), std::vector<int>(numberOfSpecies,0));
    // iterate through the new species vector

    std::vector<std::string>::iterator iterStoich;
    std::vector<std::string>::iterator iterReactants;

    for(unsigned int r=0; r<currentRSM.size(); r++){
      unsigned int stoichCounter =0;
      for(iterStoich = newRSN.begin(); iterStoich != newRSN.end(); iterStoich++, stoichCounter++){
        // for each element of the new expanded species set
        unsigned int reactCounter =0;
        for(iterReactants = currentRSN.begin(); iterReactants != currentRSN.end(); iterReactants++, reactCounter++){
          // for each element of the current smaller species set
          if(*iterStoich == *iterReactants){
            // if the species names match
            // issue when 2x + Y -> 3X coded as Y -> X, are these equivalent?
            newRSM[r][stoichCounter] = currentRSM[r][reactCounter];
          }

        }
      }
    }
    // repeat for the new reactions
    for(unsigned int r=0; r<tempRSM.size(); r++){
      unsigned int stoichCounter =0;
      for(iterStoich = newRSN.begin(); iterStoich != newRSN.end(); iterStoich++, stoichCounter++){
        // for each element of the new expanded species set
        unsigned int reactCounter =0;
        for(iterReactants = tempRSN.begin(); iterReactants != tempRSN.end(); iterReactants++, reactCounter++){
          // for each element of the current smaller species set
          if(*iterStoich == *iterReactants){
            // if the species names match
            // issue when 2x + Y -> 3X coded as Y -> X, are these equivalent?
            newRSM[r+currentRSM.size()][stoichCounter] = tempRSM[r][reactCounter];
          }
        }
      }
    }
  }
  // now have complete stoich matrix and species name vector

  // update kinetics and compositions

  // determine reaction kinetics for new reactions, append to current reactions
  std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<double>>> kineticTuple = composeKinetics(newSystem);
  tempRI = std::get<0>(kineticTuple);
  tempKV = std::get<1>(kineticTuple);
  // append the new reaction details to the current reaction details
  std::vector<std::vector<std::string>> newRI(tempRI.size() + currentRI.size());
  std::move(currentRI.begin(), currentRI.end(), newRI.begin());
  std::move(tempRI.begin(), tempRI.end(), newRI.begin() + currentRI.size());

  std::vector<std::vector<double>> newKV(tempKV.size() + currentKV.size());
  std::move(currentKV.begin(), currentKV.end(), newKV.begin());
  std::move(tempKV.begin(), tempKV.end(), newKV.begin() + currentKV.size());

  // mixture composition
  std::vector<std::string> currentMSN = getMixtureSpeciesNames();
  std::vector<double> currentMC = getMixtureInitialConcentrations();
  std::vector<std::string> newMSN = currentMSN;
  std::vector<double> newMC = currentMC;
  // check the new reaction species names for any new species not in mixture
  std::vector<std::string>::iterator iterC;
  std::vector<std::string>::iterator iterN;
  bool found = false;
  for(iterN = newRSN.begin(); iterN != newRSN.end(); iterN++){
    // for each species in the new reaciton system, search previous mixture
    found = false;
    for(iterC = currentMSN.begin(); iterC != currentMSN.end(); iterC++){
      if(*iterC == *iterN){
        // if species found in mixture do nothing
        found =true;
        break;
      }
    }
    if(!found){
      // if species is not found in the mixture, add to mixture vocab at 0 concentration
      newMSN.push_back(*iterN);
      newMC.push_back(0.0);
    }
  }
  std::vector<std::string> sortedMSN = newMSN;
  std::vector<double> sortedMC = newMC;

  sort(sortedMSN.begin(), sortedMSN.end());
  auto lastM = unique(sortedMSN.begin(), sortedMSN.end());
  sortedMSN.erase(lastM, sortedMSN.end());

  std::vector<std::string>::iterator iterUnsorted;
  std::vector<std::string>::iterator iterSorted;
  unsigned int unsortedValue = 0;

  for(iterUnsorted = newMSN.begin(); iterUnsorted != newMSN.end(); iterUnsorted++, unsortedValue++){
    unsigned int sortedValue=0;
    for(iterSorted = sortedMSN.begin(); iterSorted != sortedMSN.end(); iterSorted++, sortedValue++){
      if(*iterUnsorted == *iterSorted){
        // if match on the labels
        sortedMC[sortedValue] = newMC[unsortedValue];
        break;
      }
    }
  }
  /*
  std::cout<<"--------------------------"<<std::endl;
  printVec(newMSN);
  printVec(newMC);
  std::cout<<"Sorted: "<<std::endl;
  printVec(sortedMSN);
  printVec(sortedMC);
  std::cout<<"--------------------------"<<std::endl;
  */
  newMSN = sortedMSN;
  newMC = sortedMC;

  // set new reactions system
  setReactingSpeciesNames(newRSN);
  setReactingStoichiometry(newRSM);
  setReactionKinetics(newKV);
  setKineticLabels(newRI);
  setMixtureSpeciesNames(newMSN);
  setMixtureInitialConcentrations(newMC);

  return;
};

void reactionVessel::addReactionFromString(std::string addReactionString){
  // get current reactions
  std::vector<std::string> currentRSN = getReactingSpeciesNames();
  std::vector<std::vector<int>> currentRSM = getReactingStoichiometry();
  std::vector<std::vector<double>> currentKV = getReactionKinetics();
  std::vector<std::vector<std::string>> currentRI = getKineticLabels();

  std::vector<std::vector<int>> newRSM;

  std::vector<std::string> tempRSN;
  std::vector<std::vector<int>> tempRSM;
  std::vector<std::vector<double>> tempKV;
  std::vector<std::vector<std::string>> tempRI;

  // read new reaction from string
  std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string>> reactionTuple = parseReactionString(addReactionString);
  // vectorise the read string for passing to composeStoichMatrix
  std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string>>> reactionTupleVectorised = {reactionTuple};
  // reformat the parsed reaction
  std::tuple<std::vector<std::string>,std::vector<std::vector<int>>> stoichTuple = composeStoichMatrix(reactionTupleVectorised);
  // returns speciesArray, stoichMatrix

  tempRSN = std::get<0>(stoichTuple);
  tempRSM = std::get<1>(stoichTuple);

  // combine the species arrays
  std::vector<std::string> newRSN(tempRSN.size() + currentRSN.size());
  std::copy(currentRSN.begin(), currentRSN.end(), newRSN.begin());
  std::copy(tempRSN.begin(), tempRSN.end(), newRSN.begin() + currentRSN.size());
  // remove duplicated species name entries
  sort(newRSN.begin(), newRSN.end());
  auto last = unique(newRSN.begin(), newRSN.end());
  newRSN.erase(last, newRSN.end());


  //unsigned int numberOfReactions = newRSM.size();
  // check if new reaction intorduces new species, if so update species vectors and initialise a zero concentration
  if(newRSN.size() == currentRSN.size()){

    // no new species have appeared, append the reaction stoichiometry
    // species ordering should be the same
    // combine the stoich vectors to add new reaction(s)
    newRSM.resize(currentRSM.size() + tempRSM.size());
    std::copy(currentRSM.begin(), currentRSM.end(), newRSM.begin());
    std::copy(tempRSM.begin(), tempRSM.end(), newRSM.begin() + currentRSM.size());

  }else{
    // new species require reforming of the stochiometric matrix
    unsigned int numberOfSpecies = newRSN.size();
    unsigned int numberOfReactions = currentRSM.size() + tempRSM.size();
    //std::vector<std::vector<int>> newRSM(numberOfReactions, std::vector<int>(numberOfSpecies,0));

    newRSM.resize(numberOfReactions);
    std::fill(newRSM.begin(), newRSM.end(), std::vector<int>(numberOfSpecies,0));

    // iterate through the new species vector

    std::vector<std::string>::iterator iterStoich;
    std::vector<std::string>::iterator iterReactants;

    for(unsigned int r=0; r<currentRSM.size(); r++){
      unsigned int stoichCounter =0;
      for(iterStoich = newRSN.begin(); iterStoich != newRSN.end(); iterStoich++, stoichCounter++){
        // for each element of the new expanded species set
        unsigned int reactCounter =0;
        for(iterReactants = currentRSN.begin(); iterReactants != currentRSN.end(); iterReactants++, reactCounter++){
          // for each element of the current smaller species set

          if(*iterStoich == *iterReactants){
            // if the species names match
            // issue when 2x + Y -> 3X coded as Y -> X, are these equivalent?
            newRSM[r][stoichCounter] = currentRSM[r][reactCounter];
          }

        }
      }
    }

    // repeat for the new reactions
    for(unsigned int r=0; r<tempRSM.size(); r++){
      //std::cout<"Reaction: "<<r<<std::endl;
      unsigned int stoichCounter =0;
      for(iterStoich = newRSN.begin(); iterStoich != newRSN.end(); iterStoich++, stoichCounter++){
        // for each element of the new expanded species set
        unsigned int reactCounter =0;
        for(iterReactants = tempRSN.begin(); iterReactants != tempRSN.end(); iterReactants++, reactCounter++){
          // for each element of the current smaller species set
          if(*iterStoich == *iterReactants){
            // if the species names match
            // issue when 2x + Y -> 3X coded as Y -> X, are these equivalent?
            newRSM[r+currentRSM.size()][stoichCounter] = tempRSM[r][reactCounter];
          }
        }
      }
    }
  }
  // now have complete stoich matrix and species name vector
  // update kinetics and compositions

  // determin reaciton kinetics for new reactions, append to current reactions

  // vectorise the read string for passing to composeKinetics
  //std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string>>> reactionTupleVectorised = {reactionTuple};

  std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<double>>> kineticTuple = composeKinetics(reactionTupleVectorised);
  tempRI = std::get<0>(kineticTuple);
  tempKV = std::get<1>(kineticTuple);
  // append the new reaction details to the current reaction details
  std::vector<std::vector<std::string>> newRI(tempRI.size() + currentRI.size());
  std::move(currentRI.begin(), currentRI.end(), newRI.begin());
  std::move(tempRI.begin(), tempRI.end(), newRI.begin() + currentRI.size());

  std::vector<std::vector<double>> newKV(tempKV.size() + currentKV.size());
  std::move(currentKV.begin(), currentKV.end(), newKV.begin());
  std::move(tempKV.begin(), tempKV.end(), newKV.begin() + currentKV.size());

  // mixture composition
  std::vector<std::string> currentMSN = getMixtureSpeciesNames();
  std::vector<double> currentMC = getMixtureInitialConcentrations();
  std::vector<std::string> newMSN = currentMSN;
  std::vector<double> newMC = currentMC;
  // check the new reaction species names for any new species not in mixture
  std::vector<std::string>::iterator iterC;
  std::vector<std::string>::iterator iterN;
  bool found = false;
  for(iterN = newRSN.begin(); iterN != newRSN.end(); iterN++){
    // for each species in the new reaciton system, search previous mixture
    found = false;
    for(iterC = currentMSN.begin(); iterC != currentMSN.end(); iterC++){
      if(*iterC == *iterN){
        // if species found in mixture do nothing
        found =true;
        break;
      }
    }
    if(!found){
      // if species is not found in the mixture, add to mixture vocab at 0 concentration
      newMSN.push_back(*iterN);
      newMC.push_back(0.0);
    }
  }
  std::vector<std::string> sortedMSN = newMSN;
  std::vector<double> sortedMC = newMC;

  sort(sortedMSN.begin(), sortedMSN.end());
  auto lastM = unique(sortedMSN.begin(), sortedMSN.end());
  sortedMSN.erase(lastM, sortedMSN.end());

  std::vector<std::string>::iterator iterUnsorted;
  std::vector<std::string>::iterator iterSorted;
  unsigned int unsortedValue = 0;

  for(iterUnsorted = newMSN.begin(); iterUnsorted != newMSN.end(); iterUnsorted++, unsortedValue++){
    unsigned int sortedValue=0;
    for(iterSorted = sortedMSN.begin(); iterSorted != sortedMSN.end(); iterSorted++, sortedValue++){
      if(*iterUnsorted == *iterSorted){
        // if match on the labels
        sortedMC[sortedValue] = newMC[unsortedValue];
        break;
      }
    }
  }
  newMSN = sortedMSN;
  newMC = sortedMC;
  // set new reactions system
  setReactingSpeciesNames(newRSN);
  setReactingStoichiometry(newRSM);
  setReactionKinetics(newKV);
  setKineticLabels(newRI);
  setMixtureSpeciesNames(newMSN);
  setMixtureInitialConcentrations(newMC);

  return;
};


// composition
// settor
void reactionVessel::setMixtureSpeciesNames(std::vector<std::string> names){
  mixtureSpeciesNames = names;
  return;
};

void reactionVessel::setMixtureInitialConcentrations(std::vector<double> concentrations){
  mixtureInitialConcentrations = concentrations;

};

void reactionVessel::setMixtureConcentrations(std::vector<double> concentrations){
  mixtureConcentrations = concentrations;
};

void reactionVessel::setReactionRates(std::vector<double> rates){
  reactionRates = rates;
};

// getter
std::vector<std::string> reactionVessel::getMixtureSpeciesNames(){
  return mixtureSpeciesNames;
};

std::vector<double> reactionVessel::getMixtureInitialConcentrations(){
  return mixtureInitialConcentrations;
};

std::vector<double> reactionVessel::getReactionRates(){
  return reactionRates;
};

std::vector<double> reactionVessel::getMixtureConcentrations(){
  return mixtureConcentrations;
};

std::vector<double> reactionVessel::calculateReactionQuotient(std::vector<bool> reversible = {true}){
  // take in bool vector for reversibility
  std::vector<std::vector<int>> stoich = getReactingStoichiometry();
  std::vector<double> conc = getMixtureInitialConcentrations();

  unsigned int numberOfReactions = stoich.size();
  unsigned int numberOfSpecies = stoich[0].size();
  // initialise quotient vector
  std::vector<double> quotient(numberOfReactions, 1.0);
  double products = 1.0;
  double substrates = 1.0;
  if (reversible.size() != 1){
    for(unsigned int r=0; r<numberOfReactions; r++){
      // for each reation test if reverisble, if so calculate quotient
      if(reversible[r]){
        //std::cout<<"Reversible"<<std::endl;
        // take reaction row and depending on sign of stoich parse products and
        // substrates, calculate product of these
        products = 1.0;
        substrates = 1.0;
        for(unsigned int s =0; s<numberOfSpecies; s++){
          if(stoich[r][s] != 0){
            // species is present in the reaction
            if(stoich[r][s]>0){
              // product of the reaction
              products *= std::pow(conc[s], stoich[r][s]);
            //  std::cout<<"reaction: "<<r<<" Species "<<s<<" pro "<<products<<" contri "<<std::pow(conc[s], stoich[r][s])<<std::endl;
            //  std::cout<<"From conc "<<conc[s]<<" and stoich "<<stoich[r][s]<<std::endl;
            }else{
              // substrate of the reaction
              substrates *= std::pow(conc[s], abs(stoich[r][s]));
            //  std::cout<<"reaction: "<<r<<" Species "<<s<<" sub "<<substrates<<" contri "<<std::pow(conc[s], abs(stoich[r][s]))<<std::endl;
            //  std::cout<<"From conc "<<conc[s]<<" and stoich "<<stoich[r][s]<<std::endl;
            }
          }
        }
        if(std::abs(substrates-0.0)<1e-9){
          // substrates are not present
        //  std::cout<<"Nill substrates: "<<substrates<<std::endl;
          quotient[r] = 0.0;
        }else{
        //  std::cout<<"here: products "<<products<<" subs "<<substrates<<std::endl;
          quotient[r] = products/substrates;
        }

      }else{
        // else skip irreversible reacitons, quotinet is irrelevant
      //  std::cout<<"irreversible"<<std::endl;
      }
    }
  }else{
    // then all reactions to be considered reversible
    for(unsigned int r=0; r<numberOfReactions; r++){
      for(unsigned int s =0; s<numberOfSpecies; s++){
        if(stoich[r][s] != 0){
          // species is present in the reaction
          if(stoich[r][s]>0){
            // product of the reaction
            products *= std::pow(conc[s], stoich[r][s]);
          }else{
            // substrate of the reaction
            substrates *= std::pow(conc[s], stoich[r][s]);
          }
        }
      }
      if(std::abs(substrates-0.0)<1e-9){
        // substrates are not present
        quotient[r] = 0.0;
      }else{
        quotient[r] = products/substrates;
      }
    }
  };

  return quotient;
};

std::vector<double> reactionVessel::calculateReactionRates(double Temp = 310.0){
  // calculate the forward/ reverse reaction rates for each reaction
  // take in temperature at each reaciton step
  std::vector<std::vector<int>> stoich = getReactingStoichiometry();
  std::vector<double> conc = getMixtureInitialConcentrations();
  std::vector<std::vector<double>> kinetics = getReactionKinetics();
  std::vector<std::vector<std::string>> kineticLabels = getKineticLabels();

  unsigned int numberOfReactions = stoich.size();
  unsigned int numberOfSpecies = stoich[0].size();

  std::vector<bool> reversibleTest(numberOfReactions,false);

  // initialise vector for reaciton rates
  std::vector<double> reactionRates(numberOfReactions, 0.0);

  // initialise vector for each reaction containing current kf, kr, deltaG0
  std::vector<std::vector<double>> transformedKinetics(numberOfReactions, std::vector<double>(3,0.0));
  // can record these kinetics for plotting

  for(unsigned int r=0; r<numberOfReactions; r++){
    // test if reaction is reverisble
    if(kineticLabels[r].size() == 2){
      reversibleTest[r] = true;
      transformedKinetics[r][0] = kinetics[r][0];
      transformedKinetics[r][1] = kinetics[r][1];
      // get the deltaG value
      transformedKinetics[r][2] = KftodeltaG(kinetics[r][0], kinetics[r][1], Temp);

    }else if(kineticLabels[r][0] == "deltaG"){
      // reaction information in deltaG form
      reversibleTest[r] = true;
      transformedKinetics[r][2] = kinetics[r][0];
      // assume kr is 1.0
      transformedKinetics[r][1] = 1.0;
      // calculate kf
      transformedKinetics[r][0] = deltaGtoKf(kinetics[r][0], Temp, 1.0);

    }else{
      // reaction is irreversible
      // reversibleTest intilised to false, kr, deltaG are 0
      transformedKinetics[r][0] = kinetics[r][0];
    }
  }
  // now have Rx3 vector of kinetic values

  // calculate reaction quotient
  std::vector<double> quotient = calculateReactionQuotient(reversibleTest);
  //std::cout<<"Reaction quotient: "<<std::endl;
  //printVec(quotient);
  double products = 1.0;
  double substrates = 1.0;
  // calculate the reaciton gibbs free energy
  for(unsigned int r=0; r<numberOfReactions; r++){
    products = 1.0;
    substrates = 1.0;
    if(reversibleTest[r]){
      // if reaction is reversible
      for(unsigned int s =0; s<numberOfSpecies; s++){
        if(stoich[r][s] != 0){
          // species is present in the reaction
          if(stoich[r][s]>0){
            // product of the reaction
            //std::cout<<"Product: "<<std::pow(conc[s], stoich[r][s])<<std::endl;
            products *= std::pow(conc[s], stoich[r][s]);
          }else{
            // substrate of the reaction
            //std::cout<<"Substrate: "<<std::pow(conc[s], abs(stoich[r][s]))<<std::endl;
            substrates *= std::pow(conc[s], abs(stoich[r][s]));
          }
        }
      }
      //std::cout<<"Products: "<<products<<std::endl;
      //std::cout<<"Subst "<<substrates<<std::endl;
      //std::cout<<"here "<<deltaGtoKf(calculateGibbsFromQuotient(transformedKinetics[r][3], abs(quotient[r]), Temp), Temp, 1.0)*substrates<<std::endl;
      reactionRates[r] = deltaGtoKf(calculateGibbsFromQuotient(transformedKinetics[r][3], abs(quotient[r]), Temp), Temp, 1.0)*substrates - 1.0*products;
      //std::cout<<"reaction rate: "<<reactionRates[r]<<std::endl;
    }else{
      // reaction is irreversible
      for(unsigned int s =0; s<numberOfSpecies; s++){
        if(stoich[r][s] != 0){
          // species is present in the reaction
          if(stoich[r][s]<0){
            // product of the reaction
            substrates *= std::pow(conc[s], stoich[r][s]);
          }
        }
      }
      reactionRates[r] = transformedKinetics[r][0]*substrates;
    }
  }

  return reactionRates;
};

std::vector<double> reactionVessel::calculateConcentrationChange(){
  // calculate the change in individual species concentraiton from reaction rates
  std::vector<std::vector<int>> stoich = getReactingStoichiometry();
  std::vector<double> conc = getMixtureConcentrations();
  std::vector<double> reactionRates = calculateReactionRates();
  //std::cout<<"(conc change) Reaction rates: "<<std::endl;
  //printVec(reactionRates);
  unsigned int numberOfReactions = stoich.size();
  unsigned int numberOfSpecies = stoich[0].size();

  std::vector<double> concChange(numberOfSpecies,0.0);

  for(unsigned int s=0; s<numberOfSpecies; s++){
    for(unsigned int r=0; r<numberOfReactions; r++){
      concChange[s] += stoich[r][s]*reactionRates[r];
    }
  }
  return concChange;
};

// write methods

void reactionVessel::printReactionVesselToFile(std::string outputFilename = "reactionVessel.txt"){

  std::ofstream file;
  if(outputFilename != "reactionVessel.txt"){
    file.open(outputFilename, std::fstream::app | std::fstream::out);
  }else{
    // delete default contents file before writing.
    file.open(outputFilename, std::fstream::trunc | std::fstream::out);
  }

  if(file.is_open()){

    file<<"input reation filename: "<<inputReactionFilename<<"\n";
    file<<"output reation filename: "<<outputReactionFilename<<"\n";
    file<<"input mixture filename: "<<inputMixtureFilename<<"\n";
    file<<"output mixture filename: "<<ouputMixtureFilename<<"\n";
    file<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<"\n";
    printReactionSystemToFile(outputFilename);
    file<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<"\n";
    printVesselMixtureToFile(outputFilename);

    file.close();
  }else{
    std::cout<<"Unable to open file: "<<outputFilename<<std::endl;
  }


  return;
};

void reactionVessel::printReactionSystemToFile(std::string outputFilename="reactionSystem.txt"){

  std::ofstream file;
  // if file already exsts append to the file
  if(outputFilename != "reactionSystem.txt"){
    file.open(outputFilename, std::fstream::app | std::fstream::out);
  }else{
    // delete default contents file before writing.
    file.open(outputFilename, std::fstream::trunc | std::fstream::out);
  }

  if(file.is_open()){
    // retrieve system reactions
    std::vector<std::string> speciesNames = getReactingSpeciesNames();
    std::vector<std::vector<int>> reactionStoich = getReactingStoichiometry();
    std::vector<std::vector<double>> reactionKinetics = getReactionKinetics();
    std::vector<std::vector<std::string>> kineticLabels = getKineticLabels();

    std::string stringDelimiter = " ";
    std::string dataDelimiter = "|";
    std::string irreverDelimiter = "->";
    std::string reverDelimiter = "<->";

    unsigned int numberOfReactions = reactionKinetics.size();
    std::vector<std::string>::iterator iterR;
    std::vector<std::string>::iterator iterK;
    unsigned int stoichCounter =0;

    // print to file line by line, reaction by reaction basis
    for(unsigned int r = 0; r<numberOfReactions; r++){
      std::string leftHandString;
      std::string rightHandString;
      std::string kineticString;
      std::string stoichValue;
      std::string kineticValues;
      // pull out reactants and products
      stoichCounter=0;
      for(iterR = speciesNames.begin(); iterR != speciesNames.end(); iterR++, stoichCounter++){

        if(reactionStoich[r][stoichCounter] > 0){
          // add space for formatting
          rightHandString.append(" ");
          stoichValue = std::to_string(reactionStoich[r][stoichCounter]);
          rightHandString.append(stoichValue + *iterR);

        }else if(reactionStoich[r][stoichCounter] < 0){
          stoichValue = std::to_string(reactionStoich[r][stoichCounter]);
          leftHandString.append(stoichValue + *iterR);
          // add space for formatting
          leftHandString.append(" ");
        };
      }
      // strings now have reactants

      // determine whether the reaction was reversible
      if( kineticLabels[r].size() == 2 || kineticLabels[r][0] == "deltaG" ){
        // reaction is reversible
        leftHandString.append(reverDelimiter);
        leftHandString.append(rightHandString);
        unsigned int kineticCounter=0;
        for(iterK = kineticLabels[r].begin(); iterK != kineticLabels[r].end(); iterK++, kineticCounter++){

            kineticString.append(*iterK);
            kineticString.append(" = ");
            kineticValues = std::to_string(reactionKinetics[r][0]);
            kineticString.append(kineticValues + " ");
        }


      }else{
        // reaction is irreversible
        leftHandString.append(irreverDelimiter);
        leftHandString.append(rightHandString);
        kineticString.append("kf = ");
        kineticValues = std::to_string(reactionKinetics[r][0]);
        kineticString.append(kineticValues);
      };

      // check if last element is a space for formatting
      //if(std::strcmp(leftHandString.back(), std::string(1, ' '))){
      if(isspace(leftHandString.back())){
        // correct format
        leftHandString.append(dataDelimiter);
      }else{
        // add space
        leftHandString.append(" "+dataDelimiter);
      }
      // append kinetic string
      leftHandString.append(kineticString);

      file<<leftHandString<<"\n";
    };


    file.close();
  }else{
    std::cout<<"Unable to open file: "<<outputFilename<<std::endl;
    file.close();
  }


  return;
};

void reactionVessel::printVesselMixtureToFile(std::string outputFilename = "vesselMixture.csv"){

  std::ofstream file;
  // if file already exsts append to the file
  if(outputFilename != "vesselMixture.csv"){
    file.open(outputFilename, std::fstream::app | std::fstream::out);
  }else{
    // delete default contents file before writing.
    file.open(outputFilename, std::fstream::trunc | std::fstream::out);
  }

  if(file.is_open()){
    // retrive mixture species and concentrations
    std::vector<std::string> speciesNames = getMixtureSpeciesNames();
    std::vector<double> speciesConcentrations = getMixtureInitialConcentrations();
    unsigned int numberOfSpecies = speciesNames.size();
    // print name,cocnentration on line by line, species by species basis
    for(unsigned int line = 0; line<numberOfSpecies; line++){
      file << speciesNames[line] << "," << speciesConcentrations[line]<<'\n';
    }
    file.close();
  }else{
    std::cout<<"Unable to open file: "<<outputFilename<<std::endl;
  }
  return;
};
