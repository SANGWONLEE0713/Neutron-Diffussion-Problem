#ifndef MY_MATRIX_H
#define MY_MATRIX_H
class matrix {
public :
	int		row;
	int		col;
	double* 	dat;
	static int  inverse(const matrix &am, matrix &bm);	
public :	
	friend matrix operator+(const matrix &am, const matrix &bm);
	friend matrix operator-(const matrix &am, const matrix &bm);
	friend matrix operator*(const double  k, const matrix &am);
	friend matrix operator*(const matrix &am, const matrix &bm);
	void set_value(int r, int c, double value);
	double get_value(int r, int c);
	void set_size(int r, int c);
	matrix& operator=(const matrix &am);
	matrix& operator=(double* arr);
	void print(void);
	void print2(int r, int c);
	void save2(int r, int c, const char* fname);
	matrix();
	matrix(int r, int c);
	matrix(int r, int c, double* pval);
	~matrix(void);
};
#endif

