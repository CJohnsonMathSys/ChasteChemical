#pragma once

#include <string>
#include <tuple>
#include <vector>


std::string readReactionFilename(int, char**);

std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> >> parseReactionsFromFile(std::string);

std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> > parseReactionString(std::string);

void printVecVec(std::vector<std::vector<double>> );

void printVecVec(std::vector<std::vector<int>> );

void printVecVec(std::vector<std::vector<bool>> );

void printVecVec(std::vector<std::vector<std::string>> );

void printVec(std::vector<double>);

void printVec(std::vector<int>);

void printVec(std::vector<bool>);

void printVec(std::vector<std::string>);

//std::vector<std::vector<double>> composeStoichMatrix(std::vector<std::vector<double>> stoich, std::vector<std::vector<std::string>> species);
std::tuple<std::vector<std::string>,std::vector<std::vector<int>>> composeStoichMatrix(std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> >>);

std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<double>>> composeKinetics(std::vector<std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<int>>, std::vector<double>, std::vector<std::string> >>);

//std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<double>>> composeKinetics(std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<double>>, std::vector<double>, std::vector<std::string> >);

std::string readMixtureComposition(std::string);

std::string readMixtureComposition();

std::tuple<std::vector<std::string>, std::vector<double>> composeCompositionVectorFromFile(std::string, std::vector<std::string>);

double deltaGtoKf(double, double, double);

double KftodeltaG(double, double, double);

double calculateGibbsFromQuotient(double, double, double);
