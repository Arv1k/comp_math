#include <cstdio>
#include <cstdlib>

#include "matrix.h"


int main(int argc, char** argv) {
    if (argc < 2) {
        printf("No input file!\n");

        exit(-1);
    }

    FILE* inputFile = fopen(argv[1], "rb");

    int N = 0;
    fscanf(inputFile, "%d", &N);

    Matrix A(N, N);
    A.readFromFile(inputFile);

    fclose(inputFile);

    Matrix L = A.cholesky();
    Matrix Lt = L.transpose();

    L.DUMP("L");
    Lt.DUMP("Lt");

    Matrix vectorNum(N, 1);
    auto arr = (double*) calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        arr[i] = i + 1;
    }
    vectorNum.readArray(arr);
    vectorNum.DUMP("vectorNum");

    Matrix vectorB = A * vectorNum;
    vectorB.DUMP("vectorB");

    Matrix vectorY = L.firstEqCholesky(vectorB);
    vectorY.DUMP("vectorY");

    Matrix vectorX = Lt.secondEqCholesky(vectorY);
    vectorX.DUMP("vectorX");

    printf("norma: %.3lg\n", (vectorNum - vectorX).normaVector());
}
