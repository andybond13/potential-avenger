//potential-avenger.h
//Andrew Stershic
//DCML, Duke University
//(c) 2013

#include <stdio.h>
#include <vector>

#ifndef POTENTIAL_AVENGER_H
#define POTENTIAL_AVENGER_H

class PotentialAvenger{

public:
//dyndimresttls.m
void run();

PotentialAvenger(double in0, double in1, double in2, unsigned in3, double in4, unsigned in5);
~PotentialAvenger();


private:
double strain_rate;
double ts_refine;
double end_t;
unsigned Nelt;
double lc;
unsigned intOrder;

void printRunInfo();

//nucleate.m
void nucleate(double t, const std::vector<double>& x, std::vector<double*> phi, const std::vector<double>& xnuc, const std::vector<double>& phinuc);

//findFragments.m
void findFragments(const std::vector<double>& x, const std::vector<double>& phi, const std::vector<double>& d, unsigned& nfrags, std::vector<std::pair<unsigned,unsigned> >& fragmentList);

//checkFailureCriteria.m
void checkFailureCriteria(double t, const std::vector<double>& x, std::vector<double*> phi, std::vector<double>& criterion, std::string elemOrNodal, const std::vector<double>& qty, bool absOrAsIs, bool phiPos, double failvalue);

//analyzeDamage.m
void analyzeDamage(const std::vector<double>& x, std::vector<double*> phi, const double h, std::vector<std::vector<unsigned> >& newSegment);

};

class Matrix
{
//http://stackoverflow.com/questions/2076624/c-matrix-class
public:
    Matrix(size_t rows, size_t cols);
    double& operator()(size_t i, size_t j);
    double operator()(size_t i, size_t j) const;
    std::vector<double*> row(size_t i);
    std::vector<double> rowC(size_t i) const;
    std::vector<double*> col(size_t i);
    std::vector<double> colC(size_t i) const;

    unsigned nRows();
    unsigned nCols();    

private:
    size_t mRows;
    size_t mCols;
    std::vector<double> mData;
};

Matrix::Matrix(size_t rows, size_t cols)
: mRows(rows),
  mCols(cols),
  mData(rows * cols)
{
}

unsigned Matrix::nRows() {
    return mRows;
}

unsigned Matrix::nCols() {
    return mCols;
}

double& Matrix::operator()(size_t i, size_t j)
{
    return mData[i * mCols + j];
}

double Matrix::operator()(size_t i, size_t j) const
{
    return mData[i * mCols + j];
}

std::vector<double*> Matrix::row(size_t i)
{
    //get row
    std::vector<double*> data;
    for (unsigned j = 0; j < mCols; ++j) data.push_back(&mData[i*mCols + j]);
    return data;
}

std::vector<double> Matrix::rowC(size_t i) const
{
    //get row
    std::vector<double> data;
    for (unsigned j = 0; j < mCols; ++j) data.push_back(mData[i*mCols + j]);
    return data;
}

std::vector<double*> Matrix::col(size_t j)
{
    //get col
    std::vector<double*> data;
    for (unsigned i = 0; i < mRows; ++i) data.push_back(&mData[i*mCols + j]);
    return data;
}

std::vector<double> Matrix::colC(size_t j) const
{
    //get col
    std::vector<double> data;
    for (unsigned i = 0; j < mRows; ++i) data.push_back(mData[i*mCols + j]);
    return data;
}


#endif

