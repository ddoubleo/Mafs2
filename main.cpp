#include <iostream>
#include <cmath>
#include "MATRIX.H"
#include "MATRIX.CPP"

using namespace std;

void setMatrix(double p, N_TYPE a[ndim][ndim]) {
    float temp[8][8] = {
            {0,  -4,  -4,  7,   2,  3,  8,  7},
            {0,  -15, -1,  5,   -3, 6,  6,  -6},
            {-4, 2,   -16, 7,   0,  8,  -7, 6},
            {0,  8,   -5,  -11, 1,  0,  4,  5},
            {8,  6,   -8,  4,   27, -7, -1, 5},
            {-4, -2,  1,   2,   -8, 10, 7,  0},
            {0,  -1,  5,   2,   -8, 2,  -2, 0},
            {0,  -8,  -7,  3,   -7, -4, -8, 5}
    };
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            a[i][j] = temp[i][j];
        }
    }
    a[0][0] = p - 3;

}

void multiplyMatrix(N_TYPE a1[ndim][ndim], N_TYPE a2[ndim][ndim], N_TYPE result[ndim][ndim]) {
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            result[i][j] = 0;
            for (int k = 0; k < ndim; k++) {
                result[i][j] += a1[i][k] * a2[k][j];
            }
        }
    }
}

void printMatrix(N_TYPE a[ndim][ndim]) {
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            if (a[i][j] >= 0) {
                printf("| %e", a[i][j]);
            } else
                printf("|%e", a[i][j]);
        }
        printf("|\n");
    }
}


void subtractMatrix(N_TYPE a[ndim][ndim]) {
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            if (j == i) a[i][j] = 1 - a[i][j];
        }
    }

}

float getNorm(N_TYPE a[ndim][ndim]) {
    float sum = 0;
    float maxSum = 0;
    for (int i = 0; i < ndim; i++) {
        sum = 0;
        for (int j = 0; j < ndim; j++) {
            sum += abs(a[i][j]);
        }
        if (sum > maxSum) maxSum = sum;
    }
    return maxSum;
}

double solveMatrix(N_TYPE A_lu[ndim][ndim], N_TYPE A_rev[ndim][ndim]) {
    VECTOR(z, ndim);
    N_TYPE cond;
    int ipvt[ndim];
    N_TYPE work[ndim];
    decomp(ndim, A_lu, &cond, ipvt, work);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            z[j] = 0;
        }
        z[i] = 1;
        solve(ndim, A_lu, z, ipvt);
        for (int j = 0; j < ndim; j++) {
            A_rev[i][j] = z[j];
        }
    }
    return cond;
}

int main() {
    MATRIX(A);
    MATRIX(A_lu);
    MATRIX(A_rev);
    MATRIX(result);
    double p[5] = {1.0, 0.1, 0.01, 0.0001, 0.000001};
    for (int i = 0; i < 5; i++) {
        setMatrix(p[i], A);
        setMatrix(p[i], A_lu);

        double cond = solveMatrix(A_lu, A_rev);

        printf("\np=%1.6f\n\n", p[i]);

        multiplyMatrix(A_rev, A, result);
        subtractMatrix(result);

        printf("A Matrix\n");
        printMatrix(A);

        printf("\nInverted Matrix\n");
        printMatrix(A_rev);

        printf("\nMatrix R = E - A * A^-1\n");
        printMatrix(result);

        printf("\ncond = %f\n", cond);
        double norm = getNorm(result);
        printf("norm = %e\n", norm);
    }
}
