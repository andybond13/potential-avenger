//potential-avenger.h
//Andrew Stershic
//DCML, Duke University
//(c) 2013

#include <stdio.h>
#include <vector>
#include <objects.h>
#include <damageModel.h>

#ifndef POTENTIAL_AVENGER_H
#define POTENTIAL_AVENGER_H

class PotentialAvenger{

public:
//dyndimresttls.m
void run();

PotentialAvenger(double& in0, double& in1, double& in2, unsigned& in3, double& in4, unsigned& in5, unsigned& in6);
~PotentialAvenger();


private:
double strain_rate;
double ts_refine;
double end_t;
unsigned Nelt;
double lc;
unsigned intOrder;
unsigned printVTK;

void printRunInfo();
void printToFile(unsigned t, std::vector<double>& x, std::vector<double>& u, std::vector<double>& v, std::vector<double>& phi, std::vector<double>& d, std::vector<double>& s, std::vector<double>& e, double& Yc);

//nucleate.m
void nucleate(double t, const std::vector<double>& x, std::vector<double*> phi, const std::vector<double>& xnuc, const std::vector<double>& phinuc, std::vector<Segment>& newSegment);

//findFragments.m
void findFragments(DamageModel& dm, std::vector<Segment>& newSegment, const std::vector<double>& phi, unsigned& nfrags, std::vector<Fragment>& fragmentList);

//checkFailureCriteria.m
void checkFailureCriteria(double t, const std::vector<double>& x, std::vector<double*> phi, std::vector<double>& criterion, std::string elemOrNodal, const std::vector<double>& qty, bool absOrAsIs, bool phiPos, double failvalue, std::vector<Segment>& newSegment);

//analyzeDamage.m
void analyzeDamage(const std::vector<double>& x, std::vector<double*> phi, const double h, std::vector<Segment>& newSegment);

};



#endif

