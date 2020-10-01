#ifndef CHOLESKY_MATRIX_H
#define CHOLESKY_MATRIX_H



class Matrix {
private:
    int      rows_;
    int      columns_;
    double*  content_;
    double** matrix_;

public:
    Matrix(int rows, int columns);
    Matrix(const Matrix& that);
    ~Matrix();

    void readFromFile(FILE* file);
    void readArray(const double* array);

    void swap(Matrix& that);
    Matrix& operator=(Matrix that);
    Matrix operator*(const Matrix& that);
    Matrix operator-(const Matrix& that);
    Matrix operator+(const Matrix& that);

    Matrix transpose();
    double normaVector();

    Matrix cholesky();
    Matrix firstEqCholesky(const Matrix& vectorB);
    Matrix secondEqCholesky(const Matrix& vectorY);

    void DUMP(const char* name) const;
};


double sumSquares(const double* array, int amount);

double sumProducts(const double* array1, const double* array2, int amount);

double sumFirstEqCholesky(const double* L, double** Y, int amount);

double sumSecondEqCholesky(const double* Lt, double** X, int n, int i);

double sumRowsColumns(double** A, double** B, int columns, int i, int j);



#endif //CHOLESKY_MATRIX_H
