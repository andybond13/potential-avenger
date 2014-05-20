//objects.h
//Andrew Stershic
//DCML, Duke University
//(c) 2013

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
    unsigned begin();
    unsigned begin() const;
    unsigned second();
    unsigned end();
    unsigned penult();
    double length();
    std::vector<unsigned> indices;
    double xpeak;   //listmax
    double xmin;
    double phipeak; //phimax
    double phimin;
    double slope;
    void setPeak(const std::vector<double>& x, const std::vector<double>& phi);
    bool operator<(const Segment& in) const;

    private:
};

class Fragment{
    public:
    Fragment();
    Fragment(const Fragment& f);
    ~Fragment();
    
    double length();
    unsigned begin();
    unsigned end();
    void add(Segment* in);
    void clear();
    double localLength;

    private:
    std::vector<Segment*> fragSegs; 
};

Fragment::Fragment(){
    localLength = 0.0;
};

Fragment::Fragment(const Fragment& f){
    //copy constructor
    fragSegs.clear();
    localLength = f.localLength;
    for (unsigned i = 0; i < f.fragSegs.size(); ++i) fragSegs.push_back(f.fragSegs[i]);
};

Fragment::~Fragment(){};

double Fragment::length() {
    double len = 0;


    for (unsigned i = 0; i < fragSegs.size(); ++i) len += fragSegs[i]->length(); 
    return len + localLength;
}

unsigned Fragment::begin() {
    std::vector<unsigned> begins;
    for (unsigned i = 0; i < fragSegs.size(); ++i) begins.push_back(fragSegs[i]->begin());
    std::sort(begins.begin(),begins.end());
    return begins[0];
}

unsigned Fragment::end() {
    std::vector<unsigned> ends;
    for (unsigned i = 0; i < fragSegs.size(); ++i) ends.push_back(fragSegs[i]->end());
    std::sort(ends.begin(),ends.end());
    return ends[ends.size()-1];
}

void Fragment::add(Segment* in) {
    fragSegs.push_back(in);
    return;
}

void Fragment::clear() {
    fragSegs.clear();
    return;
}

Segment::Segment(double in1, double in2, int in3) {
    xpeak = in1;
    phipeak = in2;
    slope = in3;
    indices.clear();
};

Segment::Segment(){};
Segment::~Segment(){};

double Segment::length() {
    if (size() == 0) return 0;
    unsigned a = begin();
    unsigned b = end();
    double diff = static_cast<double>(b-a+1);
    //if (a == 0) diff -= 0.5;
    return diff;
}

bool Segment::operator<( const Segment& in ) const {
    //sort first by xpeak then by slope 
    //unsigned a = begin();
    //unsigned b = in.begin();
    if (xpeak != in.xpeak) return (xpeak < in.xpeak);
    return (slope < in.slope);
    	//return (a < b); 
}

unsigned Segment::begin() {
    std::sort(indices.begin(),indices.end());
    return indices[0];
}

unsigned Segment::begin() const{
    assert(indices.size() > 0);
    unsigned min = indices[0];
    for (unsigned i = 0; i < indices.size(); ++i) {
        if (indices[i] < min) min = indices[i];
    }
    return min;
}

unsigned Segment::second() {
    std::sort(indices.begin(),indices.end());
    return indices[1];
}

unsigned Segment::end() {
    std::sort(indices.begin(),indices.end());
    return indices[indices.size()-1];
}

unsigned Segment::penult() {
    std::sort(indices.begin(),indices.end());
    return indices[indices.size()-2];
}

unsigned Segment::size() {return indices.size();}

void Segment::setPeak(const std::vector<double>& x, const std::vector<double>& phi) {
    double phimax = -2;
    double xmax = 0;
    for (unsigned i = 0; i < indices.size(); ++i) {
        unsigned index = indices[i];
        if (phi[index] > phimax) {
            phimax = phi[index];
            xmax = x[index];
        }
    }
    assert(phimax > -2);
    this->xpeak = xmax;
    this->phipeak = phimax;

    double phimin = phimax;
    double xmin = xmax;
    for (unsigned i = 0; i < indices.size(); ++i) {
        unsigned index = indices[i];
        if (phi[index] < phimin) {
            phimin = phi[index];
            xmin = x[index];
        }
    }
    assert(phimin <= phimax);
    this->xmin = xmin;
    this->phimin = phimin;


    return;
}

#endif

