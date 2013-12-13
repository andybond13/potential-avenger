//damageModel.c
//Andrew Stershic
//DCML, Duke University
//(c) 2013

#include <damageModel.h>

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
    if (inType.compare("Linear")) {
        type = 1;
    } else if (inType.compare("Parabolic") == 0) {
        type = 2;
    } else if (inType.compare("Cubic") == 0) {
        type = 3;
    } else if (inType.compare("Ess") == 0) {
        type = 0;
    } else {
        assert(1 == 0);
    }
};

void DamageModel::assignLC(double inLC) {
    lc = inLC;
};

double DamageModel::getLC() {
    return lc;
};

int DamageModel::getType() {
    return type;
};

