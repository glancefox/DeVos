
#ifndef	MATRIX_H
#define MATRIX_H



#pragma once

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <assert.h>

//#include <afx.h>
using namespace  std;


class CMatrix
{
public:
	CMatrix(void);
	CMatrix(int nSize);
	CMatrix(int nRows, int nCols);
	CMatrix(int nSize, double value[]);
	CMatrix(const CMatrix& other);
	CMatrix(int nRows, int nCols, double value[]);
	~CMatrix(void);

protected:
	int m_NumRows;
	int m_NumColumns;
	double* m_pData;//matrix data buffer pointer


public:

	bool Init(int nRows, int nCols);
	void SetData(double value[]);

	bool MakeUnitMatrix(int nSize);
	bool SetElement(int nRow, int nCol, double value);
	double GetElement(int nRow, int nCol) const;
	int GetNumColumns(void) const;
	int GetNumRows(void) const;
	double* GetData(void) const;


	CMatrix& operator=(const CMatrix& other);
	bool operator==(const CMatrix& other) const;
	bool operator!=(const CMatrix& other) const;
	CMatrix operator+(const CMatrix& other) const;
	CMatrix operator-(const CMatrix& other) const;
	CMatrix operator*(double value) const;
	CMatrix operator*(const CMatrix& other) const;

	bool Complex_Mul(const CMatrix& AR, const CMatrix& AI, const CMatrix& BR, const CMatrix& BI, CMatrix& CR, CMatrix& CI) const;
	CMatrix Transpose(void) const;
	bool Invert_GaussJordan(void);
	bool Invert_GaussJordan(CMatrix& mtxImag);
	bool Invert_Ssgj(void);
	bool Invert_Trench(void);
	double Det_Gauss(void);
	int Rank_Gauss(void);
	bool Decom_LU(CMatrix& matrix_L, CMatrix& matrix_U);
	bool Det_Cholesky(double* Det);
	bool Decom_QR(CMatrix& matrix_Q);
	bool Decom_UV(CMatrix& matrix_U, CMatrix& matrix_V, double eps);
private:
	void ppp(double a[], double e[], double s[], double v[], int m, int n);
	void sss(double fg[], double cs[]);
public:
	bool General_Invert_UV(CMatrix& matrix_AP, CMatrix& matrix_U, CMatrix& matrix_V, double eps);
};



#endif