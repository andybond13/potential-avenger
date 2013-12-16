//potential-avenger.h
//Andrew Stershic
//DCML, Duke University
//(c) 2013

#include <stdio.h>
#include <vector>

#ifndef POTENTIAL_AVENGER_H
#define POTENTIAL_AVENGER_H

class PotentialAvenger{

private:
//nucleate.m
void nucleate(double t, std::vector<double>& x, std::vector<double>& phi, std::vector<double>& xnuc, std::vector<double>& phinuc);

//findFragments.m
void findFragments(std::vector<double>& x, std::vector<double>& phi, std::vector<double>& d, unsigned& nfrags, std::vector<std::pair<unsigned,unsigned> >& fragmentList);

//checkFailureCriteria.m
void checkFailureCriteria(double t, std::vector<double>& x, std::vector<double>& phi, std::vector<double>& criterion, std::string elemOrNodal, std::vector<double>& qty, bool absOrAsIs, bool phiPos, double failvalue);

};


#endif

