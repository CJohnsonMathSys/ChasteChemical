#include "reactionInputProcessing.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <stdlib.h>
#include <algorithm>

std::string readReactionFilename(int argc, char** argv){
  // read in a text file with reaction information, return string
  std::string filename;
  std::string fileType = ".txt";
  //std::cout<<"argc: "<<argc<<std::endl;
  //std::cout<<"argv: "<<argv<<std::endl;

  if (argc < 2){
    // no reaction file presented, default to "brusselator.txt"
    filename = "brusselator.txt";
  }else{
    //std::cout<<"Test: "<<(std::string) argv[1]<<std::endl;
    filename = (std::string) argv[1];
    if (filename.find(fileType) != std::string::npos){
      // filename string in correct form
    }else{
      // not a .txt file, append
      std::cout<< "Warning: Given file is not a .txt file, appending."<<std::endl;
      filename.append(fileType);
    }
  }

  std::cout<<"Drawing reactions from: "<<filename<<std::endl;

return filename;
};

std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> >> parseReactionsFromFile(std::string filename){

  std::string line;
  //std::string relativePath = "../Data/";
  //std::string absolutePath = "/home/chaste/projects/reactionDiffusion/apps/src/Data/";
  //std::ifstream inputFile(relativePath+filename);
  std::ifstream inputFile(filename);
  std::vector<std::string> Reaction;
  std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> >> system;
  if(inputFile.is_open()){
    // open the reaction file
    while (getline(inputFile,line)){
      // for each non-empty reation file line, parse the reactions into
      // data structures on line by line basis
      if(!line.empty()){
        //std::cout<<"String length: "<<line.length()<<std::endl;
        //std::cout<<"First space: "<<line.find_first_of(" ");
        system.push_back(parseReactionString(line));
      }
    }
    inputFile.close();
  return system;
  }else{
    std::cout<<"Unable to open file: "<<filename<<std::endl;
    exit (EXIT_FAILURE);
  }
};

std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> > parseReactionString(std::string line){
  std::string stringDelimiter = "";
  std::string dataDelimiter = "|";
  std::string irreverDelimiter = "->";
  std::string reverDelimiter = "<->";
  std::string delim;
  std::string react;
  std::string info;
  std::string tempString;
  size_t pos= line.find(dataDelimiter);
  react = line.substr(0,pos-1);
  info = line.substr(pos+2,std::string::npos);

  // test for reversibility in reaction string
  if(line.find(reverDelimiter) != std::string::npos){
    // reaction is reversible
    delim = reverDelimiter;
  }else{
    // reaction is irreversible
    delim = irreverDelimiter;
  }

  std::vector<std::vector<std::string>> species;
  std::vector<std::vector<int>> stoich;
  std::vector<std::string> complexes;
  // parse into complexes
  size_t posC = line.find(delim);
  complexes.push_back(react.substr(0,posC-1));
  complexes.push_back(react.substr(posC+3,std::string::npos));

  for(unsigned int complexNumber=0; complexNumber<complexes.size(); complexNumber++){
    std::string str=complexes[complexNumber];
//          std::cout<<"find: "<<complexes[complexNumber]<<" : "<<str.find(" + ",0)<<std::endl;
    size_t posSnew = 0;
    size_t posSold = 0;
    std::vector<std::string> reactants;
    std::vector<int> stoichVector;

    // remove potential whitespace from zeroth character if there, case <->
    if(isspace(str.c_str()[0])){
      str.erase(str.begin());
    }

    while(posSnew != std::string::npos){

      posSnew = str.find(" + ",posSold);
    //  std::cout<<"posSnew: "<<posSnew<<std::endl;
    //  std::cout<<"posSold: "<<posSold<<std::endl;
      std::string strT = str.substr(posSold,posSnew-posSold);
  //    std::cout<<"Species: "<<strT<<std::endl;
    //  std::cout<<"Size: "<<strT.size()<<std::endl;
      posSold=posSnew+3;
      // determine the stoich ratio of the reactant

      int stoichValue=1;

      if(isdigit(strT.c_str()[0])){
        stoichValue=atoi(strT.c_str());

        int i=0;
        while(isdigit(strT.c_str()[i])){i++;}
          tempString=strT.substr(i,std::string::npos);
      }else{tempString=strT;}
    //  std::cout<<"Stoich value: "<<stoichValue<<std::endl;
      if(complexNumber==0){
        stoichValue = -1*stoichValue;
      }
      reactants.push_back(tempString);
      stoichVector.push_back(stoichValue);
    }
    species.push_back(reactants);
    stoich.push_back(stoichVector);
  }

  // parse the reaction information
//        std::cout<<"String info: "<<info<<std::endl;
  std::vector<double> kineticValue;
  std::vector<std::string> reactionInfo;
  double polarity;
  double value;
  if(delim==irreverDelimiter){
    // reactions are irreversible and assume kinetic constant given
    if(info.find("kr") != std::string::npos){
      // kinetic in reverse form, invert the value
      polarity = -1.0;
    }else{polarity = 1.0;}
    size_t pos= info.find(" = ");
    value = atof(info.substr(pos+3,std::string::npos).c_str());
//          std::cout<<"Value: "<<value<<std::endl;
    kineticValue.push_back(polarity*value);
    reactionInfo.push_back("kf");
  }else{
    // reaction is reversible
    // kinetic info may be in form of reaction rates or deltaG
    if(info.find("deltaG") != std::string::npos){
        // reaction info in thermodyanmic terms
        size_t pos= info.find(" = ");
        value = atof(info.substr(pos+3,std::string::npos).c_str());
//          std::cout<<"Value: "<<value<<std::endl;
        kineticValue.push_back(value);
        reactionInfo.push_back("deltaG");
    }else{

      // reaction info given in kinetic form
      // parse kf first
      size_t posF= info.find("kf = ");
      size_t posR= info.find("kr = ");
      if(posR > posF){
        // kf is written before kr in the string
        // report back with kf first in both symmetries
        value = atof(info.substr(posF+4,posR-1).c_str());
//        std::cout<<"KF|KR KF Value: "<<value<<std::endl;
        kineticValue.push_back(value);
        reactionInfo.push_back("kf");
        value = atof(info.substr(posR+4, std::string::npos).c_str());
  //      std::cout<<"KF|KR KR Value: "<<value<<std::endl;
        kineticValue.push_back(value);
        reactionInfo.push_back("kr");
      }else{
        // kr is written before kf
        value = atof(info.substr(posF+4, std::string::npos).c_str());
  //      std::cout<<"KR|KF KF Value: "<<value<<std::endl;
        kineticValue.push_back(value);
        reactionInfo.push_back("kf");
        value = atof(info.substr(posR+4,posF-1).c_str());
  //      std::cout<<"KR|KF KR Value: "<<value<<std::endl;
        kineticValue.push_back(value);
        reactionInfo.push_back("kr");
      }
    }
  }

  return std::make_tuple(species, stoich, kineticValue, reactionInfo);
}

void printVecVec(std::vector<std::vector<std::string>> VecVec){
  for (unsigned int i = 0; i < VecVec.size(); i++)
  {
    for (unsigned int j = 0; j < VecVec[i].size(); j++)
    {
      std::cout << VecVec[i][j]<< '\t';
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;

return;
}

void printVecVec(std::vector<std::vector<double>> VecVec){
  for (unsigned int i = 0; i < VecVec.size(); i++)
  {
    for (unsigned int j = 0; j < VecVec[i].size(); j++)
    {
      std::cout << VecVec[i][j]<< '\t';
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;

return;
}

void printVecVec(std::vector<std::vector<int>> VecVec){
  for (unsigned int i = 0; i < VecVec.size(); i++)
  {
    for (unsigned int j = 0; j < VecVec[i].size(); j++)
    {
      std::cout << VecVec[i][j]<< '\t';
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;

return;
}

void printVec(std::vector<std::string> Vec){
  for (unsigned int i = 0; i < Vec.size(); i++)
  {
    std::cout << Vec[i]<< '\t';
  }
  std::cout<<std::endl;
return;
}

void printVec(std::vector<double> Vec){
  for (unsigned int i = 0; i < Vec.size(); i++)
  {
    std::cout << Vec[i]<< '\t';
  }
  std::cout<<std::endl;
return;
}

void printVec(std::vector<int> Vec){
  for (unsigned int i = 0; i < Vec.size(); i++)
  {
    std::cout << Vec[i]<< '\t';
  }
  std::cout<<std::endl;
return;
}

std::tuple<std::vector<std::string>,std::vector<std::vector<int>>> composeStoichMatrix(std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> >> system){
  // produce a stoichiometric matrix and label vector for reaction system
  unsigned int numberOfReactions=system.size();

  // for each reaction read the species list for unique entries
  std::vector<std::string> speciesArray;
//  std::vector<std::string>::iterator iter;
  for(unsigned int r=0; r<numberOfReactions; r++){
  //  printVecVec(std::get<0>(system[r]));
    //std::cout<<"place: "<< std::get<0>(system[r]).size()<<std::endl;
    for(int c=0; c<2; c++){
      // for each (head, tail) complex
      unsigned int numberReactants= std::get<0>(system[r])[c].size();
      for(unsigned int s=0; s<numberReactants; s++){
        speciesArray.push_back(std::get<0>(system[r])[c][s]);
        }
    }
  }
  sort(speciesArray.begin(),speciesArray.end());
  auto last = unique(speciesArray.begin(), speciesArray.end());
  speciesArray.erase(last, speciesArray.end());
  unsigned int numberOfSpecies = speciesArray.size();
  std::vector<std::vector<int>> stoichMatrix(numberOfReactions, std::vector<int>(numberOfSpecies,0));

  // match the stoichiometric values to the label array
  std::vector<std::string>::iterator iterReactants;
  std::vector<std::string>::iterator iterStoich;
  for(unsigned int r=0; r<numberOfReactions; r++){
    // for each reaction row
    for(int c=0; c<2; c++){
      // for each of the complesx (head, tail) of the species vector
      unsigned int stoichCounter =0;
      unsigned int reactCounter =0;
      // for each unique species label, count through the stoihc metrix labels
      for(iterStoich = speciesArray.begin(); iterStoich != speciesArray.end(); iterStoich++, stoichCounter++){
        reactCounter =0;
        for(iterReactants = std::get<0>(system[r])[c].begin(); iterReactants != std::get<0>(system[r])[c].end(); iterReactants++){

          if(*iterStoich == *iterReactants){
            // issue when 2x + Y -> 3X coded as Y -> X, are these equivalent?
            stoichMatrix[r][stoichCounter] = stoichMatrix[r][stoichCounter] + std::get<1>(system[r])[c][reactCounter];

          }
          reactCounter++;
        }
      }

    }
  }

  return std::make_tuple(speciesArray, stoichMatrix);
}

/*
std::tuple<std::vector<std::string>,std::vector<std::vector<int>>> composeStoichMatrix(std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<double>>, std::vector<double>, std::vector<std::string> > system){
  // produce a stoichiometric matrix and label vector for reaction system
  unsigned int numberOfReactions=system.size();
  // for each reaction read the species list for unique entries
  std::vector<std::string> speciesArray;
//  std::vector<std::string>::iterator iter;

  //  printVecVec(std::get<0>(system[r]));
    //std::cout<<"place: "<< std::get<0>(system[r]).size()<<std::endl;
  for(int c=0; c<2; c++){
    // for each (head, tail) complex
    unsigned int numberReactants= std::get<0>(system)[c].size();
    for(unsigned int s=0; s<numberReactants; s++){
      speciesArray.push_back(std::get<0>(system)[c][s]);
      }
  }

  sort(speciesArray.begin(),speciesArray.end());
  auto last = unique(speciesArray.begin(), speciesArray.end());
  speciesArray.erase(last, speciesArray.end());
  //std::cout<<"Print unique reactant array: "<<std::endl;
  unsigned int numberOfSpecies = speciesArray.size();
  std::vector<std::vector<double>> stoichMatrix(numberOfReactions, std::vector<int>(numberOfSpecies,0));

  // match the stoichiometric values to the label array
  std::vector<std::string>::iterator iterReactants;
  std::vector<std::string>::iterator iterStoich;

  // for each reaction row
  for(int c=0; c<2; c++){
    // for each of the complesx (head, tail) of the species vector
    unsigned int stoichCounter =0;
    unsigned int reactCounter =0;
    for(iterStoich = speciesArray.begin(); iterStoich != speciesArray.end(); iterStoich++, stoichCounter++){
      for(iterReactants = std::get<0>(system)[c].begin(); iterReactants != std::get<0>(system)[c].end(); iterReactants++){
        if(*iterStoich == *iterReactants){
            // issue when 2x + Y -> 3X coded as Y -> X, are these equivalent?
            stoichMatrix[1][stoichCounter] = stoichMatrix[1][stoichCounter] + std::get<1>(system)[c][reactCounter];
            reactCounter++;
          }
        }
      }
    }

  return std::make_tuple(speciesArray, stoichMatrix);
}
*/
std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<double>>> composeKinetics(std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> >> system){

  std::vector<std::vector<std::string>> kineticLabels;
  std::vector<std::vector<double>> kineticValues;
  int numberOfReactions = system.size();

  for(int r=0; r<numberOfReactions; r++){
    std::vector<double> reactionKineticsLine;
    std::vector<std::string> reacitonLabelLine;

    for(unsigned int v=0; v<std::get<2>(system[r]).size();v++){
      reactionKineticsLine.push_back(std::get<2>(system[r])[v]);
      reacitonLabelLine.push_back(std::get<3>(system[r])[v]);
    }
    kineticValues.push_back(reactionKineticsLine);
    kineticLabels.push_back(reacitonLabelLine);
  }


  return std::make_tuple(kineticLabels, kineticValues);
};
/*
std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<double>>> composeKinetics(std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<double>>, std::vector<double>, std::vector<std::string> > system){

  std::vector<std::vector<std::string>> kineticLabels;
  std::vector<std::vector<double>> kineticValues;

  std::vector<double> reactionKineticsLine;
  std::vector<std::string> reacitonLabelLine;

  for(unsigned int v=0; v<std::get<2>(system).size();v++){
    reactionKineticsLine.push_back(std::get<2>(system)[v]);
    reacitonLabelLine.push_back(std::get<3>(system)[v]);
  }
  kineticValues.push_back(reactionKineticsLine);
  kineticLabels.push_back(reacitonLabelLine);

  return std::make_tuple(kineticLabels, kineticValues);
};
*/
std::string readMixtureComposition(std::string filename){
  // read in a csv file with mixture information, return string
  std::string fileType = ".csv";
  if (filename.find(fileType) != std::string::npos){
    // filename string in correct form
  }else{
    // not a .txt file, append
    std::cout<< "Warning: Given file is not a .txt file, appending."<<std::endl;
    filename.append(fileType);
  }


  std::cout<<"Drawing mixture composition from: "<<filename<<std::endl;

return filename;
}

std::string readMixtureComposition(){
  // read in a csv file with mixture information, return string
return "mixture.csv";
}

std::tuple<std::vector<std::string>, std::vector<double>> composeCompositionVectorFromFile(std::string filename, std::vector<std::string> vocab){
  // read in a composition file (filename) and parse the relevant info

  std::string line;
  std::vector<std::string> speciesVector;
  // initialise cocnentration vector to 0.0 for each vocab species
  std::vector<double> concentrationVector(vocab.size(),0.0);
  std::string individualSpecies;
  double individualConcentration;
  std::string delim =",";
  //std::string relativePath = "../Data/";
  //std::ifstream inputFile(relativePath+filename);
  std::ifstream inputFile(filename);
  std::vector<std::string>::iterator iter;
  int vocabCounter=0;
  bool found = false;
  // map to vocab to the mixture vector
  speciesVector = vocab;

  if(inputFile.is_open()){
    // open the reaction file
    while (getline(inputFile,line)){
      if(!line.empty()){

        size_t pos = line.find(delim);
        individualSpecies = line.substr(0,pos);
        individualConcentration = atof(line.substr(pos+1,std::string::npos).c_str());
  //      std::cout<<"Species: "<<individualSpecies<<std::endl;
  //      std::cout<<"Concentration: "<<individualConcentration<<std::endl;
        vocabCounter=0;
        found = false;
        for(iter = vocab.begin(); iter != vocab.end(); iter++, vocabCounter++){
          if(*iter == individualSpecies){
            concentrationVector[vocabCounter] = individualConcentration;
            found =true;
            break;
          }
        }
        if(!found){
          speciesVector.push_back(individualSpecies);
          concentrationVector.push_back(individualConcentration);
        }
      }
    }
    inputFile.close();
    return std::make_tuple(speciesVector, concentrationVector);
    }else{
    std::cout<<"Unable to open file: "<<filename<<std::endl;
    exit (EXIT_FAILURE);
    }
};
