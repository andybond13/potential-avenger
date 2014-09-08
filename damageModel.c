//damageModel.c
//Andrew Stershic
//DCML, Duke University
//(c) 2013

#include <damageModel.h>
#include <math.h>
#include <iostream>

using namespace std;

DamageModel::DamageModel(){
    type = -1;
};

DamageModel::DamageModel(std::string inType, double inLC) {
    type = -1;
    assignType(inType);
    assignLC(inLC);
    assert(type != -1);
};

DamageModel::~DamageModel(){};

void DamageModel::assignType(std::string inType) {
    if (inType.compare("Linear") == 0) {
        type = 1;
    } else if (inType.compare("Parabolic") == 0) {
        type = 2;
    } else if (inType.compare("Cubic") == 0) {
        type = 3;
    } else if (inType.compare("Ess") == 0) {
        type = 0;
    } else if (inType.compare("SquareRoot") == 0) {
        type = 4;
    } else {
        assert(1 == 0);
    }
};

void DamageModel::assignLC(double inLC) {
    lc = inLC;
};

double DamageModel::dval(double phi) {
    double x = phi/lc;
    if (x < 0) return 0;
    if (x > 1) return 1;

    switch (type) {
    case 0:
        if (x <= 0.5) {
            return 2*x*x;
        } else {
            return -2*x*x + 4*x - 1;
        }
    case 1:
        return x;
    case 2:
        return 2*x-x*x;
    case 3:
        return 3*x - 3*x*x + x*x*x;
    case 4:
        return sqrt(x);
    }
    assert(1 == 0);
    return -1;
};

double DamageModel::dp(double phi) {
    //this function is the derivative of d with respect to phi
   
    double x = phi/lc;
 
    if (x < 0) return 0;
    if (x > 1) return 0;

    switch (type) {
    case 0:
        if (x <= 0.5) {
            return 4*x/lc;
        } else {
            return (-4 * x + 4)/lc;
        }
    case 1:
        return 1/lc;
    case 2:
        return 2/lc * (1-x);
    case 3:
        return 1/lc * (3-6*x+3*x*x);
    case 4:
        return 1/lc * 2*x;
    }
    assert(1 == 0);
    return -1;
};

double DamageModel::dpp(double phi) {
    //second derivative of d with respect to phi
    double x = phi/lc;

    if (x < 0) return 0;
    if (x > 1) return 0;

    switch (type) {
    case 0:
        if (x <= 0.5) {
            return 4/(lc*lc);
        } else {
            return -4/(lc*lc);
        }
    case 1:
        return 0;
    case 2:
        return -2/(lc*lc);
    case 3:
        return (1/lc)*(1/lc) * 6 * (-1+x);
    case 4:
        return (1/lc)*(1/lc)*2;
    }
    assert(1 == 0);
    return -1;
};

double DamageModel::phi(double d) {

   if (d < 0) assert(1 == 0);
   if (d == 0) return 0;

   switch (type) {
   case 0:
       if (d <= 0.5) {
           return lc * sqrt(d/2.0);
       } else {
           return lc * (1.0 - sqrt( (1.0 - d)/2.0 ));
       }
   case 1:
       return lc * d;
   case 2:
       return lc * (1.0 - sqrt(1.0 - d));
   case 3:
       return lc * (1.0 - pow(1.0 - d, 1./3.));
   case 4:
       return lc * d * d;
   }
   assert(1 == 0);
   return -1;

};

double DamageModel::getLC() {
    return lc;
};

int DamageModel::getType() {
    return type;
};

