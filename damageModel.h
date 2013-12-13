//damageModel.h
//Andrew Stershic
//DCML, Duke University
//(c) 2013

#include <stdio.h>
#include <assert.h>
#include <string>

#ifndef DAMAGEMODEL_H
#define DAMAGEMODEL_H

class DamageModel{

public:
DamageModel();

DamageModel(std::string inType, double lc);

~DamageModel();

void assignType(std::string inType);

void assignLC(double lc);

double getLC();

int getType();

private:
int type;
double lc;

}; 

#endif

