//objects.c
//Andrew Stershic
//DCML, Duke University
//(c) 2013,2014

#include <objects.h>

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

int Fragment::begin() {
    std::vector<unsigned> begins;
    for (unsigned i = 0; i < fragSegs.size(); ++i) begins.push_back(fragSegs[i]->begin());
    std::sort(begins.begin(),begins.end());
    if (begins.size() == 0) return -1;
    return begins.front();
}

unsigned Fragment::end() {
    std::vector<unsigned> ends;
    for (unsigned i = 0; i < fragSegs.size(); ++i) ends.push_back(fragSegs[i]->end());
    std::sort(ends.begin(),ends.end());
    if (ends.size() == 0) return 0;
    else return ends.back();
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
    if (a == 0) diff -= 0.5;
    return diff;
}

int Segment::begin() {
    std::sort(indices.begin(),indices.end());
    if (indices.size() == 0) return -1;
    return indices[0];
}

int Segment::begin() const{
    assert(indices.size() > 0);
    int min = indices[0];
    for (unsigned i = 0; i < indices.size(); ++i) {
        if (indices[i] < min) min = indices[i];
    }
    return min;
}

unsigned Segment::second() {
    std::sort(indices.begin(),indices.end());
    assert(indices.size() >= 2);
    return indices[1];
}

unsigned Segment::end() {
    std::sort(indices.begin(),indices.end());
    return indices.back();
}

unsigned Segment::penult() {
    std::sort(indices.begin(),indices.end());
    assert(indices.size() >= 2);
    return indices[indices.size()-2];
}

unsigned Segment::size() {
	if (indices.size() == 1) {
		if (indices[0] == -1) return 0;
    }
	return indices.size();
}

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

bool compareSegments(Segment* a, Segment* b) {
   //sort first by xpeak then by slope 
   //unsigned a = begin();
   //unsigned b = in.begin();
   if (a->xpeak != b->xpeak) return (a->xpeak < b->xpeak);
   return (a->slope < b->slope);
       //return (a < b); 
}

