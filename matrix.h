//matrix.h
//Andrew Stershic
//DCML, Duke University
//(c) 2013

#include <stdio.h>
#include <vector>

#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
//http://stackoverflow.com/questions/2076624/c-matrix-class
public:
    Matrix(size_t rows, size_t cols);
    ~Matrix();
    double& operator()(size_t i, size_t j);
    double operator()(size_t i, size_t j) const;
    std::vector<double*> row(size_t i);
    std::vector<double> rowC(size_t i) const;
    void resizeRow(size_t i, size_t j);

    unsigned nRows();
    unsigned nCols(size_t i);    

private:
    size_t mRows;
    //size_t mCols;
    //std::vector<double> mData;
    std::vector<std::vector<double> > mMatrix;
};

Matrix::Matrix(size_t rows, size_t cols)
: mRows(rows)
{
    mMatrix = std::vector<std::vector<double> >(rows,std::vector<double>(cols,0));
}

Matrix::~Matrix(){};

unsigned Matrix::nRows() {
    return mRows;
}

unsigned Matrix::nCols(size_t i) {
    return mMatrix[i].size();
}

double& Matrix::operator()(size_t i, size_t j)
{
    return mMatrix[i][j];//mData[i * mCols + j];
}

double Matrix::operator()(size_t i, size_t j) const
{
    return mMatrix[i][j];
}

void Matrix::resizeRow(size_t i, size_t j)
{
    mMatrix[i].resize(j);
    return;
}

std::vector<double*> Matrix::row(size_t i)
{
    //get row
    std::vector<double*> data;
    for (unsigned j = 0; j < nCols(i); ++j) data.push_back(&mMatrix[i][j]);
    return data;
}

std::vector<double> Matrix::rowC(size_t i) const
{
    //get row
    //std::vector<double> data;
    //for (unsigned j = 0; j < mCols; ++j) data.push_back(mData[i*mCols + j]);
    return mMatrix[i];
}


#endif

