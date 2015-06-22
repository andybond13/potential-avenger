//objects.h
//Andrew Stershic
//DCML, Duke University
//(c) 2013,2014

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>

#ifndef OBJECTS_H
#define OBJECTS_H

class Segment{

    public:
    Segment();
    Segment(double xpeak, double phipeak, int slope);
    ~Segment();

    unsigned size();
    int begin();
    int begin() const;
    unsigned second();
    unsigned end();
    unsigned penult();
    double length();
    std::vector<int> indices;
    double xpeak;   //listmax
    double xmin;
    double phipeak; //phimax
    double phimin;
    double slope;
    double YbarmYc;
    void setPeak(const std::vector<double>& x, const std::vector<double>& phi);

    private:
};
/*
class Fragment{
    public:
    Fragment();
    Fragment(const Fragment& f);
    ~Fragment();
    
    double length();
    int begin();
    unsigned end();
    void add(Segment* in);
    void clear();
    double localLength;

    private:
    std::vector<Segment*> fragSegs; 
};
*/
class SegmentComparer {
	public:
	bool operator() (Segment* a, Segment* b);
};

#endif

