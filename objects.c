//objects.c
//Andrew Stershic
//DCML, Duke University
//(c) 2013,2014

#include <objects.h>

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
	if (!is_sorted(indices.begin(),indices.end()))    std::sort(indices.begin(),indices.end());
    if (indices.size() == 0) return -1;
    return indices[0];
}

int Segment::begin() const{
	unsigned indSize = indices.size();
    assert(indSize > 0);
    int min = indices[0];
    for (unsigned i = 1; i < indSize; ++i) {
        if (indices[i] < min) min = indices[i];
    }
    return min;
}

unsigned Segment::second() {
	if (!is_sorted(indices.begin(),indices.end()))     std::sort(indices.begin(),indices.end());
    assert(indices.size() >= 2);
    return indices[1];
}

unsigned Segment::end() {
	if (!is_sorted(indices.begin(),indices.end()))    std::sort(indices.begin(),indices.end());
    return indices.back();
}

unsigned Segment::penult() {
	if (!is_sorted(indices.begin(),indices.end()))     std::sort(indices.begin(),indices.end());
    assert(indices.size() >= 2);
    return indices[indices.size()-2];
}


void Segment::beginEnd(int& min, int& max) {
	if (!is_sorted(indices.begin(),indices.end()))     std::sort(indices.begin(),indices.end());
	min = indices[0];
	max = indices.back();
	return;
}

void Segment::beginEnd(unsigned& min, unsigned& max) {
	if (!is_sorted(indices.begin(),indices.end()))     std::sort(indices.begin(),indices.end());
	min = indices[0];
	max = indices.back();
	return;
}

unsigned Segment::size() {
	unsigned indSize = indices.size();
	if (indSize == 1) {
		if (indices[0] == -1) return 0;
    }
	return indSize;
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

bool SegmentComparer::operator() (Segment* a, Segment* b) {
	return (a->begin() < b->begin());
}

