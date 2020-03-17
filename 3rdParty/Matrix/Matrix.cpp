

#include "Matrix.h"
using namespace std;



CMatrix::~CMatrix(void)
{
	if (m_pData)
	{
		delete[] m_pData;
		m_pData = NULL;
	}
}

/************************************************************************
//  initiate function
//  argument 1: rows
//	argument 2: columns
//**********************************************************************/
bool CMatrix::Init(int nRows, int nCols)
{
	if (m_pData)
	{
		delete[] m_pData;
		m_pData = NULL;
	}

	m_NumRows = nRows;
	m_NumColumns = nCols;

	int nsize = nRows*nCols;
	if (nsize<0)
	{
		return false;
	}

	m_pData = new double[nsize];//open memory space for the matrix

	if (m_pData == NULL)
	{
		return false;
	}
	//if (IsBadReadPtr(m_pData, sizeof(double)*nsize))
	//{
	//	return false;
	//}
	memset(m_pData, 0, nsize*sizeof(double));

	return true;
}


/************************************************************************
// basic construct function
//**********************************************************************/
CMatrix::CMatrix(void)
	: m_NumRows(0)
	, m_NumColumns(0)
	, m_pData(NULL)
{
	m_NumColumns = 1;
	m_NumRows = 1;
	m_pData = NULL;
	bool bsuccess = Init(m_NumRows, m_NumColumns);
	assert(bsuccess);
}


/************************************************************************
//  construct function--by the assigned size
//  argument 1: rows
//	argument 2: columns
//**********************************************************************/
CMatrix::CMatrix(int nRows, int nCols)
{
	m_NumRows = nRows;
	m_NumColumns = nCols;
	m_pData = NULL;
	bool bsuccess = Init(m_NumRows, m_NumColumns);
	assert(bsuccess);
}

/************************************************************************
//  construct function--by the assigned size and array
//  argument 1: rows
//	argument 2: columns
//**********************************************************************/
CMatrix::CMatrix(int nRows, int nCols, double value[])
{
	m_NumRows = nRows;
	m_NumColumns = nCols;
	m_pData = NULL;
	bool bSuccess = Init(m_NumRows, m_NumColumns);
	
	assert(bSuccess);

	SetData(value);
}

/************************************************************************
//  construct function--by the assigned size
//  argument: nSize- the dimension of the square matrix
//  set the items to 0
//**********************************************************************/
CMatrix::CMatrix(int nSize)
{
	m_NumRows = nSize;
	m_NumColumns = nSize;
	m_pData = NULL;
	bool bSuccess = Init(nSize, nSize);
	assert(bSuccess);
}

/************************************************************************
//  construct function--by the assigned size, construct a square matrix
//  argument: nSize- the dimension of the square matrix
//  argument: value[]- set the data source
//**********************************************************************/
CMatrix::CMatrix(int nSize, double value[])
{
	m_NumRows = nSize;
	m_NumColumns = nSize;
	m_pData = NULL;
	bool bSuccess = Init(nSize, nSize);
	assert(bSuccess);

	SetData(value);
}

/************************************************************************
//  construct function--by the assigned size, construct a square matrix
//  argument: nSize- the dimension of the square matrix
//  argument: value[]- set the data source
//**********************************************************************/
CMatrix::CMatrix(const CMatrix& other)
{
	m_NumColumns = other.GetNumColumns();
	m_NumRows = other.GetNumRows();
	m_pData = NULL;
	bool bSuccess = Init(m_NumRows, m_NumColumns);
	assert(bSuccess);

	memcpy(m_pData, other.m_pData, sizeof(double)*m_NumColumns*m_NumRows);
}

/************************************************************************
//  SerData function : set the matrix data
//  argument:value[] - the matrix data source
//  set the items to 0
//**********************************************************************/
void CMatrix::SetData(double value[])
{
	memset(m_pData, 0, m_NumRows*m_NumColumns*sizeof(double));//set memory data space
	memcpy(m_pData, value, m_NumRows*m_NumColumns*sizeof(double));//copy data into the memory buffer
}

/************************************************************************
//  GetNumColumns function : get the matrix column number
//**********************************************************************/
int CMatrix::GetNumColumns(void) const
{
	return m_NumColumns;
}

/************************************************************************
//  GetNumRows function : get the matrix row number
//**********************************************************************/
int CMatrix::GetNumRows(void) const
{
	return m_NumRows;
}

/************************************************************************
//  MakeUnitMatrix function : make a unit matrix ,dimension: nSize
//**********************************************************************/
bool CMatrix::MakeUnitMatrix(int nSize)
{
	if (!Init(nSize, nSize))
		return false;

	for (int i = 0; i<nSize; ++i)
		for (int j = 0; j<nSize; ++j)
			if (i == j)
				SetElement(i, j, 1);

	return true;
}

/************************************************************************
//  SetElement function : set value to the assigned element in matrix
//**********************************************************************/
bool CMatrix::SetElement(int nRow, int nCol, double value)
{
	if (nCol < 0 || nCol >= m_NumColumns || nRow < 0 || nRow >= m_NumRows)
		return false;						// array bounds error
	if (m_pData == NULL)
		return false;							// bad pointer error

	m_pData[nCol + nRow * m_NumColumns] = value;

	return true;
}

/************************************************************************
//  GetElement function : get value of the assigned element in matrix
//**********************************************************************/
double CMatrix::GetElement(int nRow, int nCol) const
{
	assert(nCol >= 0 && nCol < m_NumColumns && nRow >= 0 && nRow < m_NumRows); // array bounds error
	assert(m_pData);							// bad pointer error
	return m_pData[nCol + nRow * m_NumColumns];
}

/************************************************************************
//  GetData function : get matrix data
//**********************************************************************/
double* CMatrix::GetData(void) const
{
	return m_pData;
}

/************************************************************************
//  operator= function :  initiate the operator =
//  assign data to the matrix
//**********************************************************************/
CMatrix& CMatrix::operator=(const CMatrix& other)
{
	if (&other != this)
	{
		bool bSuccess = Init(other.GetNumRows(), other.GetNumColumns());
		assert(bSuccess);

		// copy the pointer
		memcpy(m_pData, other.m_pData, sizeof(double)*m_NumColumns*m_NumRows);
	}

	// finally return a reference to ourselves
	return *this;
}

/************************************************************************
//  operator== function :  initiate the operator ==
//  to see if the matrix is the same size
//**********************************************************************/
bool CMatrix::operator==(const CMatrix& other) const
{
	// check if the columns and rows are the same size
	if (m_NumColumns != other.GetNumColumns() || m_NumRows != other.GetNumRows())
		return false;

	for (int i = 0; i<m_NumRows; ++i)
	{
		for (int j = 0; j<m_NumColumns; ++j)
		{
			if (GetElement(i, j) != other.GetElement(i, j))
				return false;
		}
	}

	return true;
}

/************************************************************************
//  operator== function :  initiate the operator !=
//  to see if the matrix is  not the same size
//**********************************************************************/
bool CMatrix::operator!=(const CMatrix& other) const
{
	return !(*this == other);
}

/************************************************************************
//  operator+ function :  initiate the operator +
//  matrix addition
//**********************************************************************/
CMatrix CMatrix::operator+(const CMatrix& other) const
{
	//check the size of the matrix
	assert(m_NumColumns == other.GetNumColumns() && m_NumRows == other.GetNumRows());

	// construct the result matrix
	CMatrix	result(*this);		// copy construction
	// matrix addition
	for (int i = 0; i < m_NumRows; ++i)
	{
		for (int j = 0; j < m_NumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) + other.GetElement(i, j));
	}

	return result;
}

/************************************************************************
//  operator- function :  initiate the operator -
//  matrix minus
//**********************************************************************/
CMatrix CMatrix::operator-(const CMatrix& other) const
{
	//check the size of the matrix
	assert(m_NumColumns == other.GetNumColumns() && m_NumRows == other.GetNumRows());

	// construct the result matrix
	CMatrix	result(*this);		// copy ourselves
	// carry out the minus operation
	for (int i = 0; i < m_NumRows; ++i)
	{
		for (int j = 0; j < m_NumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) - other.GetElement(i, j));
	}

	return result;
}

/************************************************************************
//  operator* function :  initiate the operator *
//  matrix multiply a real number
//**********************************************************************/
CMatrix CMatrix::operator*(double value) const
{
	// construct result matrix
	CMatrix	result(*this);		// copy ourselves
	// carry out number multiply operation
	for (int i = 0; i < m_NumRows; ++i)
	{
		for (int j = 0; j < m_NumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) * value);
	}

	return result;
}

/************************************************************************
//  operator* function :  initiate the operator *
//  two real number matrix multiply operation
//**********************************************************************/
CMatrix CMatrix::operator*(const CMatrix& other) const
{
	// check if the two matrix is the same size
	assert(m_NumColumns == other.GetNumRows());

	// construct the object we are going to return
	CMatrix	result(m_NumRows, other.GetNumColumns());
	//**********************************************************
	// matrix multiply :
	//
	// [A][B][C]   [G][H]     [A*G + B*I + C*K][A*H + B*J + C*L]
	// [D][E][F] * [I][J] =   [D*G + E*I + F*K][D*H + E*J + F*L]
	//             [K][L]
	//**********************************************************
	double	value;
	for (int i = 0; i < result.GetNumRows(); ++i)
	{
		for (int j = 0; j < other.GetNumColumns(); ++j)
		{
			value = 0.0;
			for (int k = 0; k < m_NumColumns; ++k)
			{
				value += GetElement(i, k) * other.GetElement(k, j);
			}

			result.SetElement(i, j, value);
		}
	}

	return result;
}



//*******************************************************************************************
//   A*B=C
//	arguments:
//	AR(BR)(CR): REAL part matrix of matrix A(B)(C)
//	AI(BI)(CI): IMAGINE part matrix of matrix A(B)(C)
//*******************************************************************************************
bool CMatrix::Complex_Mul(const CMatrix& AR, const CMatrix& AI, const CMatrix& BR, const CMatrix& BI, CMatrix& CR, CMatrix& CI) const
{
	// check if the same size
	if (AR.GetNumColumns() != AI.GetNumColumns() ||
		AR.GetNumRows() != AI.GetNumRows() ||
		BR.GetNumColumns() != BI.GetNumColumns() ||
		BR.GetNumRows() != BI.GetNumRows() ||
		AR.GetNumColumns() != BR.GetNumRows())
		return false;
	//construct the REAL and IMAGINE matrix of the two complex matrix
	CMatrix mtxCR(AR.GetNumRows(), BR.GetNumColumns()), mtxCI(AR.GetNumRows(), BR.GetNumColumns());
	//complex matrix multiply
	for (int i = 0; i<AR.GetNumRows(); ++i)
	{
		for (int j = 0; j<BR.GetNumColumns(); ++j)
		{
			double vr = 0;
			double vi = 0;
			for (int k = 0; k<AR.GetNumColumns(); ++k)
			{
				double p = AR.GetElement(i, k) * BR.GetElement(k, j);
				double q = AI.GetElement(i, k) * BI.GetElement(k, j);
				double s = (AR.GetElement(i, k) + AI.GetElement(i, k)) * (BR.GetElement(k, j) + BI.GetElement(k, j));
				vr += p - q;
				vi += s - p - q;
			}
			mtxCR.SetElement(i, j, vr);
			mtxCI.SetElement(i, j, vi);
		}
	}

	CR = mtxCR;
	CI = mtxCI;

	return true;
}

/************************************************************************
//  Transpose function :  transpose the matrix
//**********************************************************************/
CMatrix CMatrix::Transpose(void) const
{
	//construct result matrix 
	CMatrix	Trans(m_NumColumns, m_NumRows);

	//transpose every elements
	for (int i = 0; i < m_NumRows; ++i)
	{
		for (int j = 0; j < m_NumColumns; ++j)
			Trans.SetElement(j, i, GetElement(i, j));
	}

	return Trans;
}

/************************************************************************
//Invert_GaussJordan function:get the invert matrix using Gaussian-Jordan method
//REAL matrix invert
//**********************************************************************/
bool CMatrix::Invert_GaussJordan(void)
{
	int *pnRow, *pnCol, i, j, k, l, u, v;
	double d = 0, p = 0;

	//assign memory space
	pnRow = new int[m_NumColumns];
	pnCol = new int[m_NumColumns];
	if (pnRow == NULL || pnCol == NULL)
		return false;

	//eliminate elements
	for (k = 0; k <= m_NumColumns - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= m_NumColumns - 1; i++)
		{
			for (j = k; j <= m_NumColumns - 1; j++)
			{
				l = i*m_NumColumns + j; p = fabs(m_pData[l]);
				if (p>d)
				{
					d = p;
					pnRow[k] = i;
					pnCol[k] = j;
				}
			}
		}

		//failure
		if (d == 0.0)
		{
			delete[] pnRow;
			delete[] pnCol;
			return false;
		}

		if (pnRow[k] != k)//remember the exchange order
		{
			for (j = 0; j <= m_NumColumns - 1; j++)
			{
				u = k*m_NumColumns + j;
				v = pnRow[k] * m_NumColumns + j;
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}

		if (pnCol[k] != k)
		{
			for (i = 0; i <= m_NumColumns - 1; i++)
			{
				u = i*m_NumColumns + k;
				v = i*m_NumColumns + pnCol[k];
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}

		l = k*m_NumColumns + k;
		m_pData[l] = 1.0 / m_pData[l];
		for (j = 0; j <= m_NumColumns - 1; j++)
		{
			if (j != k)
			{
				u = k*m_NumColumns + j;
				m_pData[u] = m_pData[u] * m_pData[l];
			}
		}

		for (i = 0; i <= m_NumColumns - 1; i++)
		{
			if (i != k)
			{
				for (j = 0; j <= m_NumColumns - 1; j++)
				{
					if (j != k)
					{
						u = i*m_NumColumns + j;
						m_pData[u] = m_pData[u] - m_pData[i*m_NumColumns + k] * m_pData[k*m_NumColumns + j];
					}
				}
			}
		}

		for (i = 0; i <= m_NumColumns - 1; i++)
		{
			if (i != k)
			{
				u = i*m_NumColumns + k;
				m_pData[u] = -m_pData[u] * m_pData[l];
			}
		}
	}

	//reorder the exchange order
	for (k = m_NumColumns - 1; k >= 0; k--)
	{
		if (pnCol[k] != k)
		{
			for (j = 0; j <= m_NumColumns - 1; j++)
			{
				u = k*m_NumColumns + j;
				v = pnCol[k] * m_NumColumns + j;
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}

		if (pnRow[k] != k)
		{
			for (i = 0; i <= m_NumColumns - 1; i++)
			{
				u = i*m_NumColumns + k;
				v = i*m_NumColumns + pnRow[k];
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}
	}

	//release memory
	delete[] pnRow;
	delete[] pnCol;

	//success
	return true;
}

/************************************************************************
//Invert_GaussJordan function:get the invert matrix using Gaussian-Jordan method
//complex matrix invert operation
//**********************************************************************/
bool CMatrix::Invert_GaussJordan(CMatrix& mtxImag)
{
	int *pnRow, *pnCol, i, j, k, l, u, v, w;
	double p, q, s, t, d, b;

	//open memory space
	pnRow = new int[m_NumColumns];
	pnCol = new int[m_NumColumns];
	if (pnRow == NULL || pnCol == NULL)
		return false;

	//eliminate elements
	for (k = 0; k <= m_NumColumns - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= m_NumColumns - 1; i++)
		{
			for (j = k; j <= m_NumColumns - 1; j++)
			{
				u = i*m_NumColumns + j;
				p = m_pData[u] * m_pData[u] + mtxImag.m_pData[u] * mtxImag.m_pData[u];
				if (p>d)
				{
					d = p;
					pnRow[k] = i;
					pnCol[k] = j;
				}
			}
		}

		//failure
		if (d == 0.0)
		{
			delete[] pnRow;
			delete[] pnCol;
			return(0);
		}

		if (pnRow[k] != k)
		{
			for (j = 0; j <= m_NumColumns - 1; j++)
			{
				u = k*m_NumColumns + j;
				v = pnRow[k] * m_NumColumns + j;
				t = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = t;
				t = mtxImag.m_pData[u];
				mtxImag.m_pData[u] = mtxImag.m_pData[v];
				mtxImag.m_pData[v] = t;
			}
		}

		if (pnCol[k] != k)
		{
			for (i = 0; i <= m_NumColumns - 1; i++)
			{
				u = i*m_NumColumns + k;
				v = i*m_NumColumns + pnCol[k];
				t = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = t;
				t = mtxImag.m_pData[u];
				mtxImag.m_pData[u] = mtxImag.m_pData[v];
				mtxImag.m_pData[v] = t;
			}
		}

		l = k*m_NumColumns + k;
		m_pData[l] = m_pData[l] / d; mtxImag.m_pData[l] = -mtxImag.m_pData[l] / d;
		for (j = 0; j <= m_NumColumns - 1; j++)
		{
			if (j != k)
			{
				u = k*m_NumColumns + j;
				p = m_pData[u] * m_pData[l];
				q = mtxImag.m_pData[u] * mtxImag.m_pData[l];
				s = (m_pData[u] + mtxImag.m_pData[u])*(m_pData[l] + mtxImag.m_pData[l]);
				m_pData[u] = p - q;
				mtxImag.m_pData[u] = s - p - q;
			}
		}

		for (i = 0; i <= m_NumColumns - 1; i++)
		{
			if (i != k)
			{
				v = i*m_NumColumns + k;
				for (j = 0; j <= m_NumColumns - 1; j++)
				{
					if (j != k)
					{
						u = k*m_NumColumns + j;
						w = i*m_NumColumns + j;
						p = m_pData[u] * m_pData[v];
						q = mtxImag.m_pData[u] * mtxImag.m_pData[v];
						s = (m_pData[u] + mtxImag.m_pData[u])*(m_pData[v] + mtxImag.m_pData[v]);
						t = p - q;
						b = s - p - q;
						m_pData[w] = m_pData[w] - t;
						mtxImag.m_pData[w] = mtxImag.m_pData[w] - b;
					}
				}
			}
		}

		for (i = 0; i <= m_NumColumns - 1; i++)
		{
			if (i != k)
			{
				u = i*m_NumColumns + k;
				p = m_pData[u] * m_pData[l];
				q = mtxImag.m_pData[u] * mtxImag.m_pData[l];
				s = (m_pData[u] + mtxImag.m_pData[u])*(m_pData[l] + mtxImag.m_pData[l]);
				m_pData[u] = q - p;
				mtxImag.m_pData[u] = p + q - s;
			}
		}
	}

	//reorder the major elements exchange order
	for (k = m_NumColumns - 1; k >= 0; k--)
	{
		if (pnCol[k] != k)
		{
			for (j = 0; j <= m_NumColumns - 1; j++)
			{
				u = k*m_NumColumns + j;
				v = pnCol[k] * m_NumColumns + j;
				t = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = t;
				t = mtxImag.m_pData[u];
				mtxImag.m_pData[u] = mtxImag.m_pData[v];
				mtxImag.m_pData[v] = t;
			}
		}

		if (pnRow[k] != k)
		{
			for (i = 0; i <= m_NumColumns - 1; i++)
			{
				u = i*m_NumColumns + k;
				v = i*m_NumColumns + pnRow[k];
				t = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = t;
				t = mtxImag.m_pData[u];
				mtxImag.m_pData[u] = mtxImag.m_pData[v];
				mtxImag.m_pData[v] = t;
			}
		}
	}

	//release memory
	delete[] pnRow;
	delete[] pnCol;

	//if success
	return true;
}

/************************************************************************
//Invert_Ssgj function:get the invert (Symmetric positive definite matrix)
//                      using Gaussian-Jordan method
//**********************************************************************/
bool CMatrix::Invert_Ssgj(void)
{

	int i, j, k, m;
	double w, g, *pTmp;

	//temp memory
	pTmp = new double[m_NumColumns];

	//process by rows
	for (k = 0; k <= m_NumColumns - 1; k++)
	{
		w = m_pData[0];
		if (w == 0.0)
		{
			delete[] pTmp;
			return false;
		}

		m = m_NumColumns - k - 1;
		for (i = 1; i <= m_NumColumns - 1; i++)
		{
			g = m_pData[i*m_NumColumns];
			pTmp[i] = g / w;
			if (i <= m)
				pTmp[i] = -pTmp[i];
			for (j = 1; j <= i; j++)
				m_pData[(i - 1)*m_NumColumns + j - 1] = m_pData[i*m_NumColumns + j] + g*pTmp[j];
		}

		m_pData[m_NumColumns*m_NumColumns - 1] = 1.0 / w;
		for (i = 1; i <= m_NumColumns - 1; i++)
			m_pData[(m_NumColumns - 1)*m_NumColumns + i - 1] = pTmp[i];
	}

	//adjust row and column
	for (i = 0; i <= m_NumColumns - 2; i++)
		for (j = i + 1; j <= m_NumColumns - 1; j++)
			m_pData[i*m_NumColumns + j] = m_pData[j*m_NumColumns + i];

	//release temp memory
	delete[] pTmp;

	return true;
}

/************************************************************************
//Invert_Trench function:get the invert Toeplitz matrix
//                      using Trench method
//**********************************************************************/
bool CMatrix::Invert_Trench(void)
{
	int i, j, k;
	double a, s, *t, *tt, *c, *r, *p;

	//upper triangle elements
	t = new double[m_NumColumns];
	//lower triangle elements
	tt = new double[m_NumColumns];

	//assign the L and U triangle elements
	for (i = 0; i<m_NumColumns; ++i)
	{
		t[i] = GetElement(0, i);
		tt[i] = GetElement(i, 0);
	}

	//temp  memory buffer
	c = new double[m_NumColumns];
	r = new double[m_NumColumns];
	p = new double[m_NumColumns];

	//if it is not a Toeplitz matrix, then return
	if (t[0] == 0.0)
	{
		delete[] t;
		delete[] tt;
		delete[] c;
		delete[] r;
		delete[] p;
		return false;
	}

	a = t[0];
	c[0] = tt[1] / t[0];
	r[0] = t[1] / t[0];

	for (k = 0; k <= m_NumColumns - 3; k++)
	{
		s = 0.0;
		for (j = 1; j <= k + 1; j++)
			s = s + c[k + 1 - j] * tt[j];

		s = (s - tt[k + 2]) / a;
		for (i = 0; i <= k; i++)
			p[i] = c[i] + s*r[k - i];

		c[k + 1] = -s;
		s = 0.0;
		for (j = 1; j <= k + 1; j++)
			s = s + r[k + 1 - j] * t[j];

		s = (s - t[k + 2]) / a;
		for (i = 0; i <= k; i++)
		{
			r[i] = r[i] + s*c[k - i];
			c[k - i] = p[k - i];
		}

		r[k + 1] = -s;
		a = 0.0;
		for (j = 1; j <= k + 2; j++)
			a = a + t[j] * c[j - 1];

		a = t[0] - a;

		//process failure
		if (a == 0.0)
		{
			delete[] t;
			delete[] tt;
			delete[] c;
			delete[] r;
			delete[] p;
			return false;
		}
	}

	m_pData[0] = 1.0 / a;
	for (i = 0; i <= m_NumColumns - 2; i++)
	{
		k = i + 1;
		j = (i + 1)*m_NumColumns;
		m_pData[k] = -r[i] / a;
		m_pData[j] = -c[i] / a;
	}

	for (i = 0; i <= m_NumColumns - 2; i++)
	{
		for (j = 0; j <= m_NumColumns - 2; j++)
		{
			k = (i + 1)*m_NumColumns + j + 1;
			m_pData[k] = m_pData[i*m_NumColumns + j] - c[i] * m_pData[j + 1];
			m_pData[k] = m_pData[k] + c[m_NumColumns - j - 2] * m_pData[m_NumColumns - i - 1];
		}
	}

	//release temp memory
	delete[] t;
	delete[] tt;
	delete[] c;
	delete[] r;
	delete[] p;

	return true;
}

/************************************************************************
//Det_Gauss function: get the delta of a matrix using Gaussian method
//**********************************************************************/
double CMatrix::Det_Gauss(void)
{
	int i, j, k, is, js, l, u, v;
	double f, det, q, d;

	//initiate
	f = 1.0;
	det = 1.0;

	//eliminate elements
	for (k = 0; k <= m_NumColumns - 2; k++)
	{
		q = 0.0;
		for (i = k; i <= m_NumColumns - 1; i++)
		{
			for (j = k; j <= m_NumColumns - 1; j++)
			{
				l = i*m_NumColumns + j;
				d = fabs(m_pData[l]);
				if (d>q)
				{
					q = d;
					is = i;
					js = j;
				}
			}
		}

		if (q == 0.0)
		{
			det = 0.0;
			return(det);
		}

		if (is != k)
		{
			f = -f;
			for (j = k; j <= m_NumColumns - 1; j++)
			{
				u = k*m_NumColumns + j;
				v = is*m_NumColumns + j;
				d = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = d;
			}
		}

		if (js != k)
		{
			f = -f;
			for (i = k; i <= m_NumColumns - 1; i++)
			{
				u = i*m_NumColumns + js;
				v = i*m_NumColumns + k;
				d = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = d;
			}
		}

		l = k*m_NumColumns + k;
		det = det*m_pData[l];
		for (i = k + 1; i <= m_NumColumns - 1; i++)
		{
			d = m_pData[i*m_NumColumns + k] / m_pData[l];
			for (j = k + 1; j <= m_NumColumns - 1; j++)
			{
				u = i*m_NumColumns + j;
				m_pData[u] = m_pData[u] - d*m_pData[k*m_NumColumns + j];
			}
		}
	}

	//operate and get result
	det = f*det*m_pData[m_NumColumns*m_NumColumns - 1];

	return(det);
}

/************************************************************************
//Rank_Gauss function: get the rank of a matrix using Gaussian method
//**********************************************************************/
int CMatrix::Rank_Gauss(void)
{
	int i, j, k, nn, is, js, l, ll, u, v;
	double q, d;

	// rank is less equal to matrix row or column
	nn = m_NumRows;
	if (m_NumRows >= m_NumColumns)
		nn = m_NumColumns;

	k = 0;

	//eliminate elements
	for (l = 0; l <= nn - 1; l++)
	{
		q = 0.0;
		for (i = l; i <= m_NumRows - 1; i++)
		{
			for (j = l; j <= m_NumColumns - 1; j++)
			{
				ll = i*m_NumColumns + j;
				d = fabs(m_pData[ll]);
				if (d>q)
				{
					q = d;
					is = i;
					js = j;
				}
			}
		}

		if (q == 0.0)
			return(k);

		k = k + 1;
		if (is != l)
		{
			for (j = l; j <= m_NumColumns - 1; j++)
			{
				u = l*m_NumColumns + j;
				v = is*m_NumColumns + j;
				d = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = d;
			}
		}
		if (js != l)
		{
			for (i = l; i <= m_NumRows - 1; i++)
			{
				u = i*m_NumColumns + js;
				v = i*m_NumColumns + l;
				d = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = d;
			}
		}

		ll = l*m_NumColumns + l;
		for (i = l + 1; i <= m_NumColumns - 1; i++)
		{
			d = m_pData[i*m_NumColumns + l] / m_pData[ll];
			for (j = l + 1; j <= m_NumColumns - 1; j++)
			{
				u = i*m_NumColumns + j;
				m_pData[u] = m_pData[u] - d*m_pData[l*m_NumColumns + j];
			}
		}
	}

	return(k);
}


/************************************************************************
//Decom_LU function: decompose a matrix to LOWER and UPPER matrix
//						---using Gaussian method
//Symmetric positive definite matrix
//the original matrix would be the Q matrix after the process
//**********************************************************************/
bool CMatrix::Decom_LU(CMatrix& matrix_L, CMatrix& matrix_U)
{
	int i, j, k, w, v, ll;

	//initial the result matrix
	if (!matrix_L.Init(m_NumColumns, m_NumColumns) || !matrix_U.Init(m_NumColumns, m_NumColumns))
		return false;

	for (k = 0; k <= m_NumColumns - 2; k++)
	{
		ll = k*m_NumColumns + k;
		if (m_pData[ll] == 0.0)
			return false;

		for (i = k + 1; i <= m_NumColumns - 1; i++)
		{
			w = i*m_NumColumns + k;
			m_pData[w] = m_pData[w] / m_pData[ll];
		}

		for (i = k + 1; i <= m_NumColumns - 1; i++)
		{
			w = i*m_NumColumns + k;
			for (j = k + 1; j <= m_NumColumns - 1; j++)
			{
				v = i*m_NumColumns + j;
				m_pData[v] = m_pData[v] - m_pData[w] * m_pData[k*m_NumColumns + j];
			}
		}
	}

	for (i = 0; i <= m_NumColumns - 1; i++)
	{
		for (j = 0; j<i; j++)
		{
			w = i*m_NumColumns + j;
			matrix_L.m_pData[w] = m_pData[w];
			matrix_U.m_pData[w] = 0.0;
		}

		w = i*m_NumColumns + i;
		matrix_L.m_pData[w] = 1.0;
		matrix_U.m_pData[w] = m_pData[w];

		for (j = i + 1; j <= m_NumColumns - 1; j++)
		{
			w = i*m_NumColumns + j;
			matrix_L.m_pData[w] = 0.0;
			matrix_U.m_pData[w] = m_pData[w];
		}
	}

	return true;
}

/************************************************************************
//Det_Cholesky function: decompose and get Delta of a Symmetric positive
definite matrix using Cholesky way
argument : det -- the Delta of the matrix
//**********************************************************************/
bool CMatrix::Det_Cholesky(double* Det)
{
	int i, j, k, u, l;
	double d;

	//check if meet the process condition
	if (m_pData[0] <= 0.0)
		return false;

	//Cholesky decomposition

	m_pData[0] = sqrt(m_pData[0]);
	d = m_pData[0];

	for (i = 1; i <= m_NumColumns - 1; i++)
	{
		u = i*m_NumColumns;
		m_pData[u] = m_pData[u] / m_pData[0];
	}

	for (j = 1; j <= m_NumColumns - 1; j++)
	{
		l = j*m_NumColumns + j;
		for (k = 0; k <= j - 1; k++)
		{
			u = j*m_NumColumns + k;
			m_pData[l] = m_pData[l] - m_pData[u] * m_pData[u];
		}

		if (m_pData[l] <= 0.0)
			return false;

		m_pData[l] = sqrt(m_pData[l]);
		d = d*m_pData[l];

		for (i = j + 1; i <= m_NumColumns - 1; i++)
		{
			u = i*m_NumColumns + j;
			for (k = 0; k <= j - 1; k++)
				m_pData[u] = m_pData[u] - m_pData[i*m_NumColumns + k] * m_pData[j*m_NumColumns + k];

			m_pData[u] = m_pData[u] / m_pData[l];
		}
	}

	//rank valuation
	*Det = d*d;

	//Lower triangle matrix
	for (i = 0; i <= m_NumColumns - 2; i++)
		for (j = i + 1; j <= m_NumColumns - 1; j++)
			m_pData[i*m_NumColumns + j] = 0.0;

	return true;
}

/************************************************************************
//Decom_QR function: decompose matrix to matrix_Q and  matrix_R,
argument : matrix_Q -- the Q matrix of result
matrix_R stored in the original matrix
//**********************************************************************/
bool CMatrix::Decom_QR(CMatrix& matrix_Q)
{
	int i, j, k, l, nn, p, jj;
	double u, alpha, w, t;

	if (m_NumRows < m_NumColumns)
		return false;

	//initiate matrix Q
	if (!matrix_Q.Init(m_NumRows, m_NumRows))
		return false;
	//set Diagonal elements  to 1
	for (i = 0; i <= m_NumRows - 1; i++)
	{
		for (j = 0; j <= m_NumRows - 1; j++)
		{
			l = i*m_NumRows + j;
			matrix_Q.m_pData[l] = 0.0;
			if (i == j)
				matrix_Q.m_pData[l] = 1.0;
		}
	}

	//start decomposition
	nn = m_NumColumns;
	if (m_NumRows == m_NumColumns)
		nn = m_NumRows - 1;

	for (k = 0; k <= nn - 1; k++)
	{
		u = 0.0;
		l = k*m_NumColumns + k;
		for (i = k; i <= m_NumRows - 1; i++)
		{
			w = fabs(m_pData[i*m_NumColumns + k]);
			if (w>u)
				u = w;
		}

		alpha = 0.0;
		for (i = k; i <= m_NumRows - 1; i++)
		{
			t = m_pData[i*m_NumColumns + k] / u;
			alpha = alpha + t*t;
		}

		if (m_pData[l]>0.0)
			u = -u;

		alpha = u*sqrt(alpha);
		if (alpha == 0.0)
			return false;

		u = sqrt(2.0*alpha*(alpha - m_pData[l]));
		if ((u + 1.0) != 1.0)
		{
			m_pData[l] = (m_pData[l] - alpha) / u;
			for (i = k + 1; i <= m_NumRows - 1; i++)
			{
				p = i*m_NumColumns + k;
				m_pData[p] = m_pData[p] / u;
			}

			for (j = 0; j <= m_NumRows - 1; j++)
			{
				t = 0.0;
				for (jj = k; jj <= m_NumRows - 1; jj++)
					t = t + m_pData[jj*m_NumColumns + k] * matrix_Q.m_pData[jj*m_NumRows + j];

				for (i = k; i <= m_NumRows - 1; i++)
				{
					p = i*m_NumRows + j;
					matrix_Q.m_pData[p] = matrix_Q.m_pData[p] - 2.0*t*m_pData[i*m_NumColumns + k];
				}
			}

			for (j = k + 1; j <= m_NumColumns - 1; j++)
			{
				t = 0.0;

				for (jj = k; jj <= m_NumRows - 1; jj++)
					t = t + m_pData[jj*m_NumColumns + k] * m_pData[jj*m_NumColumns + j];

				for (i = k; i <= m_NumRows - 1; i++)
				{
					p = i*m_NumColumns + j;
					m_pData[p] = m_pData[p] - 2.0*t*m_pData[i*m_NumColumns + k];
				}
			}

			m_pData[l] = alpha;
			for (i = k + 1; i <= m_NumRows - 1; i++)
				m_pData[i*m_NumColumns + k] = 0.0;
		}
	}
	//adjust elements
	for (i = 0; i <= m_NumRows - 2; i++)
	{
		for (j = i + 1; j <= m_NumRows - 1; j++)
		{
			p = i*m_NumRows + j;
			l = j*m_NumRows + i;
			t = matrix_Q.m_pData[p];
			matrix_Q.m_pData[p] = matrix_Q.m_pData[l];
			matrix_Q.m_pData[l] = t;
		}
	}

	return true;
}

/************************************************************************
//Decom_UV function:  matrix UV decomposition
arguments:
matrix_U: stores the U matrix
matrix_V: stores the V matrix
eps: the close error of the process
//**********************************************************************/
bool CMatrix::Decom_UV(CMatrix& matrix_U, CMatrix& matrix_V, double eps /*= 0.000001*/)
{
	int i, j, k, l, it, ll, kk, ix, iy, mm, nn, iz, m1, ks;
	double d, dd, t, sm, sm1, em1, sk, ek, b, c, shh, fg[2], cs[2];
	double *s, *e, *w;

	int m = m_NumRows;
	int n = m_NumColumns;

	// initiate matrix_U and matrix_V
	if (!matrix_U.Init(m, m) || !matrix_V.Init(n, n))
		return false;

	//temp buffer
	//int ka = max(m, n) + 1;
	int ka = m > n ? m : n + 1;
	
	s = new double[ka];
	e = new double[ka];
	w = new double[ka];

	// 指定迭代次数为60
	it = 60;
	k = n;

	if (m - 1<n)
		k = m - 1;

	l = m;
	if (n - 2<m)
		l = n - 2;
	if (l<0)
		l = 0;

	//iteration and process 
	ll = k;
	if (l>k)
		ll = l;
	if (ll >= 1)
	{
		for (kk = 1; kk <= ll; kk++)
		{
			if (kk <= k)
			{
				d = 0.0;
				for (i = kk; i <= m; i++)
				{
					ix = (i - 1)*n + kk - 1;
					d = d + m_pData[ix] * m_pData[ix];
				}

				s[kk - 1] = sqrt(d);
				if (s[kk - 1] != 0.0)
				{
					ix = (kk - 1)*n + kk - 1;
					if (m_pData[ix] != 0.0)
					{
						s[kk - 1] = fabs(s[kk - 1]);
						if (m_pData[ix]<0.0)
							s[kk - 1] = -s[kk - 1];
					}

					for (i = kk; i <= m; i++)
					{
						iy = (i - 1)*n + kk - 1;
						m_pData[iy] = m_pData[iy] / s[kk - 1];
					}

					m_pData[ix] = 1.0 + m_pData[ix];
				}

				s[kk - 1] = -s[kk - 1];
			}

			if (n >= kk + 1)
			{
				for (j = kk + 1; j <= n; j++)
				{
					if ((kk <= k) && (s[kk - 1] != 0.0))
					{
						d = 0.0;
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1)*n + kk - 1;
							iy = (i - 1)*n + j - 1;
							d = d + m_pData[ix] * m_pData[iy];
						}

						d = -d / m_pData[(kk - 1)*n + kk - 1];
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1)*n + j - 1;
							iy = (i - 1)*n + kk - 1;
							m_pData[ix] = m_pData[ix] + d*m_pData[iy];
						}
					}

					e[j - 1] = m_pData[(kk - 1)*n + j - 1];
				}
			}

			if (kk <= k)
			{
				for (i = kk; i <= m; i++)
				{
					ix = (i - 1)*m + kk - 1;
					iy = (i - 1)*n + kk - 1;
					matrix_U.m_pData[ix] = m_pData[iy];
				}
			}

			if (kk <= l)
			{
				d = 0.0;
				for (i = kk + 1; i <= n; i++)
					d = d + e[i - 1] * e[i - 1];

				e[kk - 1] = sqrt(d);
				if (e[kk - 1] != 0.0)
				{
					if (e[kk] != 0.0)
					{
						e[kk - 1] = fabs(e[kk - 1]);
						if (e[kk]<0.0)
							e[kk - 1] = -e[kk - 1];
					}

					for (i = kk + 1; i <= n; i++)
						e[i - 1] = e[i - 1] / e[kk - 1];

					e[kk] = 1.0 + e[kk];
				}

				e[kk - 1] = -e[kk - 1];
				if ((kk + 1 <= m) && (e[kk - 1] != 0.0))
				{
					for (i = kk + 1; i <= m; i++)
						w[i - 1] = 0.0;

					for (j = kk + 1; j <= n; j++)
						for (i = kk + 1; i <= m; i++)
							w[i - 1] = w[i - 1] + e[j - 1] * m_pData[(i - 1)*n + j - 1];

					for (j = kk + 1; j <= n; j++)
					{
						for (i = kk + 1; i <= m; i++)
						{
							ix = (i - 1)*n + j - 1;
							m_pData[ix] = m_pData[ix] - w[i - 1] * e[j - 1] / e[kk];
						}
					}
				}

				for (i = kk + 1; i <= n; i++)
					matrix_V.m_pData[(i - 1)*n + kk - 1] = e[i - 1];
			}
		}
	}

	mm = n;
	if (m + 1<n)
		mm = m + 1;
	if (k<n)
		s[k] = m_pData[k*n + k];
	if (m<mm)
		s[mm - 1] = 0.0;
	if (l + 1<mm)
		e[l] = m_pData[l*n + mm - 1];

	e[mm - 1] = 0.0;
	nn = m;
	if (m>n)
		nn = n;
	if (nn >= k + 1)
	{
		for (j = k + 1; j <= nn; j++)
		{
			for (i = 1; i <= m; i++)
				matrix_U.m_pData[(i - 1)*m + j - 1] = 0.0;
			matrix_U.m_pData[(j - 1)*m + j - 1] = 1.0;
		}
	}

	if (k >= 1)
	{
		for (ll = 1; ll <= k; ll++)
		{
			kk = k - ll + 1;
			iz = (kk - 1)*m + kk - 1;
			if (s[kk - 1] != 0.0)
			{
				if (nn >= kk + 1)
				{
					for (j = kk + 1; j <= nn; j++)
					{
						d = 0.0;
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1)*m + kk - 1;
							iy = (i - 1)*m + j - 1;
							d = d + matrix_U.m_pData[ix] * matrix_U.m_pData[iy] / matrix_U.m_pData[iz];
						}

						d = -d;
						for (i = kk; i <= m; i++)
						{
							ix = (i - 1)*m + j - 1;
							iy = (i - 1)*m + kk - 1;
							matrix_U.m_pData[ix] = matrix_U.m_pData[ix] + d*matrix_U.m_pData[iy];
						}
					}
				}

				for (i = kk; i <= m; i++)
				{
					ix = (i - 1)*m + kk - 1;
					matrix_U.m_pData[ix] = -matrix_U.m_pData[ix];
				}

				matrix_U.m_pData[iz] = 1.0 + matrix_U.m_pData[iz];
				if (kk - 1 >= 1)
				{
					for (i = 1; i <= kk - 1; i++)
						matrix_U.m_pData[(i - 1)*m + kk - 1] = 0.0;
				}
			}
			else
			{
				for (i = 1; i <= m; i++)
					matrix_U.m_pData[(i - 1)*m + kk - 1] = 0.0;
				matrix_U.m_pData[(kk - 1)*m + kk - 1] = 1.0;
			}
		}
	}

	for (ll = 1; ll <= n; ll++)
	{
		kk = n - ll + 1;
		iz = kk*n + kk - 1;

		if ((kk <= l) && (e[kk - 1] != 0.0))
		{
			for (j = kk + 1; j <= n; j++)
			{
				d = 0.0;
				for (i = kk + 1; i <= n; i++)
				{
					ix = (i - 1)*n + kk - 1;
					iy = (i - 1)*n + j - 1;
					d = d + matrix_V.m_pData[ix] * matrix_V.m_pData[iy] / matrix_V.m_pData[iz];
				}

				d = -d;
				for (i = kk + 1; i <= n; i++)
				{
					ix = (i - 1)*n + j - 1;
					iy = (i - 1)*n + kk - 1;
					matrix_V.m_pData[ix] = matrix_V.m_pData[ix] + d*matrix_V.m_pData[iy];
				}
			}
		}

		for (i = 1; i <= n; i++)
			matrix_V.m_pData[(i - 1)*n + kk - 1] = 0.0;

		matrix_V.m_pData[iz - n] = 1.0;
	}

	for (i = 1; i <= m; i++)
		for (j = 1; j <= n; j++)
			m_pData[(i - 1)*n + j - 1] = 0.0;

	m1 = mm;
	it = 60;
	while (true)
	{
		if (mm == 0)
		{
			ppp(m_pData, e, s, matrix_V.m_pData, m, n);
			return true;
		}
		if (it == 0)
		{
			ppp(m_pData, e, s, matrix_V.m_pData, m, n);
			return false;
		}

		kk = mm - 1;
		while ((kk != 0) && (fabs(e[kk - 1]) != 0.0))
		{
			d = fabs(s[kk - 1]) + fabs(s[kk]);
			dd = fabs(e[kk - 1]);
			if (dd>eps*d)
				kk = kk - 1;
			else
				e[kk - 1] = 0.0;
		}

		if (kk == mm - 1)
		{
			kk = kk + 1;
			if (s[kk - 1]<0.0)
			{
				s[kk - 1] = -s[kk - 1];
				for (i = 1; i <= n; i++)
				{
					ix = (i - 1)*n + kk - 1;
					matrix_V.m_pData[ix] = -matrix_V.m_pData[ix];
				}
			}

			while ((kk != m1) && (s[kk - 1]<s[kk]))
			{
				d = s[kk - 1];
				s[kk - 1] = s[kk];
				s[kk] = d;
				if (kk<n)
				{
					for (i = 1; i <= n; i++)
					{
						ix = (i - 1)*n + kk - 1;
						iy = (i - 1)*n + kk;
						d = matrix_V.m_pData[ix];
						matrix_V.m_pData[ix] = matrix_V.m_pData[iy];
						matrix_V.m_pData[iy] = d;
					}
				}

				if (kk<m)
				{
					for (i = 1; i <= m; i++)
					{
						ix = (i - 1)*m + kk - 1;
						iy = (i - 1)*m + kk;
						d = matrix_U.m_pData[ix];
						matrix_U.m_pData[ix] = matrix_U.m_pData[iy];
						matrix_U.m_pData[iy] = d;
					}
				}

				kk = kk + 1;
			}

			it = 60;
			mm = mm - 1;
		}
		else
		{
			ks = mm;
			while ((ks>kk) && (fabs(s[ks - 1]) != 0.0))
			{
				d = 0.0;
				if (ks != mm)
					d = d + fabs(e[ks - 1]);
				if (ks != kk + 1)
					d = d + fabs(e[ks - 2]);

				dd = fabs(s[ks - 1]);
				if (dd>eps*d)
					ks = ks - 1;
				else
					s[ks - 1] = 0.0;
			}

			if (ks == kk)
			{
				kk = kk + 1;
				d = fabs(s[mm - 1]);
				t = fabs(s[mm - 2]);
				if (t>d)
					d = t;

				t = fabs(e[mm - 2]);
				if (t>d)
					d = t;

				t = fabs(s[kk - 1]);
				if (t>d)
					d = t;

				t = fabs(e[kk - 1]);
				if (t>d)
					d = t;

				sm = s[mm - 1] / d;
				sm1 = s[mm - 2] / d;
				em1 = e[mm - 2] / d;
				sk = s[kk - 1] / d;
				ek = e[kk - 1] / d;
				b = ((sm1 + sm)*(sm1 - sm) + em1*em1) / 2.0;
				c = sm*em1;
				c = c*c;
				shh = 0.0;

				if ((b != 0.0) || (c != 0.0))
				{
					shh = sqrt(b*b + c);
					if (b<0.0)
						shh = -shh;

					shh = c / (b + shh);
				}

				fg[0] = (sk + sm)*(sk - sm) - shh;
				fg[1] = sk*ek;
				for (i = kk; i <= mm - 1; i++)
				{
					sss(fg, cs);
					if (i != kk)
						e[i - 2] = fg[0];

					fg[0] = cs[0] * s[i - 1] + cs[1] * e[i - 1];
					e[i - 1] = cs[0] * e[i - 1] - cs[1] * s[i - 1];
					fg[1] = cs[1] * s[i];
					s[i] = cs[0] * s[i];

					if ((cs[0] != 1.0) || (cs[1] != 0.0))
					{
						for (j = 1; j <= n; j++)
						{
							ix = (j - 1)*n + i - 1;
							iy = (j - 1)*n + i;
							d = cs[0] * matrix_V.m_pData[ix] + cs[1] * matrix_V.m_pData[iy];
							matrix_V.m_pData[iy] = -cs[1] * matrix_V.m_pData[ix] + cs[0] * matrix_V.m_pData[iy];
							matrix_V.m_pData[ix] = d;
						}
					}

					sss(fg, cs);
					s[i - 1] = fg[0];
					fg[0] = cs[0] * e[i - 1] + cs[1] * s[i];
					s[i] = -cs[1] * e[i - 1] + cs[0] * s[i];
					fg[1] = cs[1] * e[i];
					e[i] = cs[0] * e[i];

					if (i<m)
					{
						if ((cs[0] != 1.0) || (cs[1] != 0.0))
						{
							for (j = 1; j <= m; j++)
							{
								ix = (j - 1)*m + i - 1;
								iy = (j - 1)*m + i;
								d = cs[0] * matrix_U.m_pData[ix] + cs[1] * matrix_U.m_pData[iy];
								matrix_U.m_pData[iy] = -cs[1] * matrix_U.m_pData[ix] + cs[0] * matrix_U.m_pData[iy];
								matrix_U.m_pData[ix] = d;
							}
						}
					}
				}

				e[mm - 2] = fg[0];
				it = it - 1;
			}
			else
			{
				if (ks == mm)
				{
					kk = kk + 1;
					fg[1] = e[mm - 2];
					e[mm - 2] = 0.0;
					for (ll = kk; ll <= mm - 1; ll++)
					{
						i = mm + kk - ll - 1;
						fg[0] = s[i - 1];
						sss(fg, cs);
						s[i - 1] = fg[0];
						if (i != kk)
						{
							fg[1] = -cs[1] * e[i - 2];
							e[i - 2] = cs[0] * e[i - 2];
						}

						if ((cs[0] != 1.0) || (cs[1] != 0.0))
						{
							for (j = 1; j <= n; j++)
							{
								ix = (j - 1)*n + i - 1;
								iy = (j - 1)*n + mm - 1;
								d = cs[0] * matrix_V.m_pData[ix] + cs[1] * matrix_V.m_pData[iy];
								matrix_V.m_pData[iy] = -cs[1] * matrix_V.m_pData[ix] + cs[0] * matrix_V.m_pData[iy];
								matrix_V.m_pData[ix] = d;
							}
						}
					}
				}
				else
				{
					kk = ks + 1;
					fg[1] = e[kk - 2];
					e[kk - 2] = 0.0;
					for (i = kk; i <= mm; i++)
					{
						fg[0] = s[i - 1];
						sss(fg, cs);
						s[i - 1] = fg[0];
						fg[1] = -cs[1] * e[i - 1];
						e[i - 1] = cs[0] * e[i - 1];
						if ((cs[0] != 1.0) || (cs[1] != 0.0))
						{
							for (j = 1; j <= m; j++)
							{
								ix = (j - 1)*m + i - 1;
								iy = (j - 1)*m + kk - 2;
								d = cs[0] * matrix_U.m_pData[ix] + cs[1] * matrix_U.m_pData[iy];
								matrix_U.m_pData[iy] = -cs[1] * matrix_U.m_pData[ix] + cs[0] * matrix_U.m_pData[iy];
								matrix_U.m_pData[ix] = d;
							}
						}
					}
				}
			}
		}
	}

	return true;
}

/************************************************************************
//ppp function: called by the Decom_UV function
//**********************************************************************/
void CMatrix::ppp(double a[], double e[], double s[], double v[], int m, int n)
{
	int i, j, p, q;
	double d;

	if (m >= n)
		i = n;
	else
		i = m;

	for (j = 1; j <= i - 1; j++)
	{
		a[(j - 1)*n + j - 1] = s[j - 1];
		a[(j - 1)*n + j] = e[j - 1];
	}

	a[(i - 1)*n + i - 1] = s[i - 1];
	if (m<n)
		a[(i - 1)*n + i] = e[i - 1];

	for (i = 1; i <= n - 1; i++)
	{
		for (j = i + 1; j <= n; j++)
		{
			p = (i - 1)*n + j - 1;
			q = (j - 1)*n + i - 1;
			d = v[p];
			v[p] = v[q];
			v[q] = d;
		}
	}
}


/************************************************************************
//sss function: called by the Decom_UV function
//**********************************************************************/
void CMatrix::sss(double fg[], double cs[])
{
	double r, d;

	if ((fabs(fg[0]) + fabs(fg[1])) == 0.0)
	{
		cs[0] = 1.0;
		cs[1] = 0.0;
		d = 0.0;
	}
	else
	{
		d = sqrt(fg[0] * fg[0] + fg[1] * fg[1]);
		if (fabs(fg[0])>fabs(fg[1]))
		{
			d = fabs(d);
			if (fg[0]<0.0)
				d = -d;
		}
		if (fabs(fg[1]) >= fabs(fg[0]))
		{
			d = fabs(d);
			if (fg[1]<0.0)
				d = -d;
		}

		cs[0] = fg[0] / d;
		cs[1] = fg[1] / d;
	}

	r = 1.0;
	if (fabs(fg[0])>fabs(fg[1]))
		r = cs[1];
	else if (cs[0] != 0.0)
		r = 1.0 / cs[0];

	fg[0] = d;
	fg[1] = r;
}

/************************************************************************
//General_Invert_UV function:get the General Invert matrix
//                      using UV decomposition method
//**********************************************************************/
bool CMatrix::General_Invert_UV(CMatrix& matrix_AP, CMatrix& matrix_U, CMatrix& matrix_V, double eps)
{
	int i, j, k, l, t, p, q, f;

	// call matrix decompose use UV decomposition method
	if (!Decom_UV(matrix_U, matrix_V, eps))
		return false;

	int m = m_NumRows;
	int n = m_NumColumns;

	//initiate general invert matrix
	if (!matrix_AP.Init(n, m))
		return false;

	//compute general invert matrix

	j = n;
	if (m<n)
		j = m;
	j = j - 1;
	k = 0;
	while ((k <= j) && (m_pData[k*n + k] != 0.0))
		k = k + 1;

	k = k - 1;
	for (i = 0; i <= n - 1; i++)
	{
		for (j = 0; j <= m - 1; j++)
		{
			t = i*m + j;
			matrix_AP.m_pData[t] = 0.0;
			for (l = 0; l <= k; l++)
			{
				f = l*n + i;
				p = j*m + l;
				q = l*n + l;
				matrix_AP.m_pData[t] = matrix_AP.m_pData[t] + matrix_V.m_pData[f] * matrix_U.m_pData[p] / m_pData[q];
			}
		}
	}

	return true;
}
