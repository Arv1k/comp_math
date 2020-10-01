#include <algorithm>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "matrix.h"


Matrix::Matrix(int rows, int columns): rows_(rows), columns_(columns) {
    content_ = (double*) malloc(rows_ * columns_ * sizeof(content_));
    for (int i = 0; i < rows_ * columns_; i++) {
        content_[i] = 0;
    }

    matrix_ = (double**) calloc(rows, sizeof(matrix_));
    for (int i = 0; i < rows_; i++) {
        matrix_[i] = &content_[i * columns_];
    }
}

Matrix::Matrix(const Matrix& that) {
    rows_    = that.rows_;
    columns_ = that.columns_;

    content_ = (double*) calloc(rows_ * columns_, sizeof(content_));
    memcpy(content_, that.content_,
    sizeof(*that.content_) * that.rows_ * that.columns_);

    matrix_ = (double**) calloc(rows_, sizeof(matrix_));
    memcpy(matrix_, that.matrix_, sizeof(*that.matrix_) * that.rows_);
}

Matrix::~Matrix() {
    rows_    = -1;
    columns_ = -1;

    free(content_);
    free(matrix_);
}

void Matrix::readFromFile(FILE* file) {
    for (int i = 0; i < rows_ * columns_; i++) {
        fscanf(file, "%lg", &content_[i]);
    }
}

void Matrix::readArray(const double* array) {
    memcpy(content_, array, rows_ * columns_ * sizeof(*array));
}

void Matrix::swap(Matrix& that) {
    std::swap(rows_,    that.rows_);
    std::swap(columns_, that.columns_);
    std::swap(content_, that.content_);
    std::swap(matrix_,  that.matrix_);
}

Matrix& Matrix::operator=(Matrix that) {
    swap(that);

    return *this;
}

Matrix Matrix::operator*(const Matrix& that) {
    if (columns_ != that.rows_) {
        throw std::invalid_argument("Not suitable matrices! (multiplication)\n");
    }

    Matrix Multi(rows_, that.columns_);

    for (int i = 0; i < Multi.rows_; i++) {
        for (int j = 0; j < Multi.columns_; j++) {
            Multi.matrix_[i][j] =
            sumRowsColumns(matrix_, that.matrix_, columns_, i, j);
        }
    }

    return Multi;
}

Matrix Matrix::operator-(const Matrix& that) {
    if (columns_ != that.columns_ || rows_ != that.rows_) {
        throw std::invalid_argument("Not suitable matrices! (subtraction)\n");
    }

    Matrix Sub(rows_, columns_);

    for (int i = 0; i < Sub.rows_; i++) {
        for (int j = 0; j < Sub.columns_; j++) {
            Sub.matrix_[i][j] = matrix_[i][j] - that.matrix_[i][j];
        }
    }

    return Sub;
}

Matrix Matrix::operator+(const Matrix& that) {
    if (columns_ != that.columns_ || rows_ != that.rows_) {
        throw std::invalid_argument("Not suitable matrices! (addition)\n");
    }

    Matrix Sub(rows_, columns_);

    for (int i = 0; i < Sub.rows_; i++) {
        for (int j = 0; j < Sub.columns_; j++) {
            Sub.matrix_[i][j] = matrix_[i][j] + that.matrix_[i][j];
        }
    }

    return Sub;
}

Matrix Matrix::transpose() {
    Matrix Lt(rows_, columns_);

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++) {
            Lt.matrix_[i][j] = matrix_[j][i];
        }
    }

    return Lt;
}

double Matrix::normaVector() {
    if (columns_ != 1) {
        throw std::invalid_argument("Not a vector! (norma)\n");
    }

    double sum = 0;

    for (int i = 0; i < rows_; i++) {
        sum += matrix_[i][0] * matrix_[i][0];
    }

    return sqrt(sum);
}

Matrix Matrix::cholesky() {
    if (rows_ != columns_) {
        throw std::invalid_argument("Not suitable matrix! (cholesky)\n");
    }

    Matrix L(rows_, columns_);

    if (matrix_[0][0] < 0) {
        throw std::invalid_argument("Root of a negative number! (cholesky)\n");
    }
    L.matrix_[0][0] = sqrt(matrix_[0][0]);

    if (L.matrix_[0][0] == 0) {
        throw std::invalid_argument("Division by zero! (cholesky)\n");
    }
    for (int i = 1; i < rows_; i++) {
        L.matrix_[i][0] = matrix_[i][0] / L.matrix_[0][0];
    }

    for (int i = 1; i < rows_; i++) {
        if ((matrix_[i][i] - sumSquares(L.matrix_[i], i)) < 0) {
            throw
            std::invalid_argument("Root of a negative number! (cholesky)\n");
        }
        L.matrix_[i][i] =
        sqrt(matrix_[i][i] - sumSquares(L.matrix_[i], i));

        for (int j = i + 1; j < columns_; j++) {
            L.matrix_[j][i] =
            matrix_[j][i] - sumProducts(L.matrix_[i], L.matrix_[j], i);

            if (L.matrix_[i][i] == 0) {
                throw std::invalid_argument("Division by zero! (cholesky)\n");
            }
            L.matrix_[j][i] /= L.matrix_[i][i];
        }
    }

    return L;
}

Matrix Matrix::firstEqCholesky(const Matrix& vectorB) {
    Matrix vectorY(vectorB.rows_, vectorB.columns_);

    if (matrix_[0][0] == 0) {
        throw std::invalid_argument("Division by zero! (cholesky)\n");
    }
    vectorY.matrix_[0][0] = vectorB.matrix_[0][0] / matrix_[0][0];

    for (int i = 1; i < vectorY.rows_; i++) {
        vectorY.matrix_[i][0] =
        vectorB.matrix_[i][0] -
        sumFirstEqCholesky(matrix_[i], vectorY.matrix_, i);

        if (matrix_[i][i] == 0) {
            throw std::invalid_argument("Division by zero! (cholesky)\n");
        }
        vectorY.matrix_[i][0] /= matrix_[i][i];
    }

    return vectorY;
}

Matrix Matrix::secondEqCholesky(const Matrix& vectorY) {
    Matrix vectorX(vectorY.rows_, vectorY.columns_);

    if (matrix_[rows_ - 1][rows_ - 1] == 0) {
        throw std::invalid_argument("Division by zero! (cholesky)\n");
    }
    vectorX.matrix_[rows_ - 1][0] =
    vectorY.matrix_[rows_ - 1][0] / matrix_[rows_ - 1][rows_ - 1];

    for (int i = rows_ - 2; i >= 0; i--) {
        vectorX.matrix_[i][0] =
        vectorY.matrix_[i][0] -
        sumSecondEqCholesky(matrix_[i], vectorX.matrix_, rows_, i);

        if (matrix_[i][i] == 0) {
            throw std::invalid_argument("Division by zero! (cholesky)\n");
        }
        vectorX.matrix_[i][0] /= matrix_[i][i];
    }

    return vectorX;
}

void Matrix::DUMP(const char* name) const {
    printf("%s\n", name);

    printf("rows: %d\n", rows_);
    printf("columns: %d\n", columns_);

    printf("matrix:\n");
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++) {
            printf("%6.3lf ", matrix_[i][j]);
        }

        printf("\n");
    }

    printf("\n");
}

double sumSquares(const double* array, int amount) {
    double sum = 0;

    for (int i = 0; i < amount; i++) {
        sum += array[i] * array[i];
    }

    return sum;
}

double sumProducts(const double* array1, const double* array2, int amount) {
    double sum = 0;

    for (int i = 0; i < amount; i++) {
        sum += array1[i] * array2[i];
    }

    return sum;
}

double sumFirstEqCholesky(const double* L, double** Y, int amount) {
    double sum = 0;

    for (int i = 0; i < amount; i++) {
        sum += L[i] * Y[i][0];
    }

    return sum;
}

double sumSecondEqCholesky(const double* Lt, double** X, int n, int i) {
    double sum = 0;

    for (; i < n; i++) {
        sum += Lt[i] * X[i][0];
    }

    return sum;
}

double sumRowsColumns(double** A, double** B, int columns, int i, int j) {
    double sum = 0;

    for (int k = 0; k < columns; k++) {
        sum += A[i][k] * B[k][j];
    }

    return sum;
}
