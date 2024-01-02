#define _GNU_SOURCE
#include "symnmf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#define EPS 1e-4
#define BETA 0.5
#define MAX_ITER 300

void print_vector(double vec[], int d){
    int i;
    for(i = 0; i < d - 1; i++){
        printf("%.4f,", vec[i]);
    }
    printf("%.4f\n", vec[d-1]);
}

double frobenius_norm_squared(double **A, int n, int k){
    int i;
    int j;
    double ret = 0;
    for(i = 0; i < n; i++){
        for (j = 0; j < k; j++){
            ret += A[i][j] * A[i][j];
        }
    }
    return ret;
}

double **matrix_transpose(double **A, int n, int k){
    int i;
    int j;
    double **A_T = malloc(k * sizeof(double*));
    assert(A_T!=NULL);
    for(i = 0; i < k; i++){
        A_T[i] = malloc(n * sizeof(double));
        assert(A_T[i]!=NULL);
        for(j = 0; j < n; j++){
            A_T[i][j] = A[j][i];
        }
    }
    return A_T;
}

void free_matrix(double **A, int n){
    int i;
    for (i = 0; i<n;i++){
        free(A[i]);
    }
    free(A);
}

double euclidian_norm_squared(double *x1, double *x2, int d){
    int i;
    double ret = 0;
    for(i = 0; i < d; i++){
        ret += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    return ret;
}

double **matrix_multiplication(double **A, double **B, int n, int m, int r){
    int i;
    int j;
    int k;
    double **C = malloc(n * sizeof(double*));
    assert(C!=NULL);
    for(i = 0; i<n; i++){
        C[i] = malloc(r * sizeof(double));
        assert(C[i]!=NULL);
        for(j =0;j<r;j++){
            C[i][j] = 0;
            for(k=0; k<m; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

double **sym(double** X, int n, int d) {
    int i;
    int j;
    double **A = (double**)malloc(n * sizeof(double*));
    assert(A!=NULL);
    
    for(i = 0; i < n; i++){
        A[i] = (double*)malloc(n * sizeof(double));
        assert(A[i]!=NULL);
        for(j = 0; j < n; j++){
            if(i==j){
                A[i][j] = 0;
            }
            else {
                A[i][j] = exp(-1 * euclidian_norm_squared(X[i], X[j], d)/2);
            }
        }
    }
    return A;
}

double **ddg(double** X, int n, int d) {
    int i;
    int j;
    int k;
    double **A;
    double **D = malloc(n * sizeof(double*));
    assert(D!=NULL);
    A = sym(X, n, d);
    for(i = 0; i < n; i++){
        D[i] = malloc(n * sizeof(double));
        assert(D[i]!=NULL);
        for(j = 0; j < n; j++){
            if(i!=j){
                D[i][j]=0;
            }
            else{
                D[i][i] = 0;
                for(k = 0; k < n; k++){
                    D[i][i] += A[i][k];
                }
            }
        }
    }
    free_matrix(A,n);
    return D;
}

double **norm(double** X, int n, int d) {
    int i;
    int j;
    double **A;
    double **D;
    double **DA;
    double **ret;
    double **D_minus_sqrt = malloc(n * sizeof(double*));
    assert(D_minus_sqrt!=NULL);
    D = ddg(X,n,d);
    A = sym(X,n,d);
    for(i = 0; i < n; i++){
        D_minus_sqrt[i] = malloc(n * sizeof(double));
        assert(D_minus_sqrt[i]!=NULL);
        for(j = 0; j < n; j++){
            D_minus_sqrt[i][j] = 0;
            if(i==j){
                D_minus_sqrt[i][j] = 1 / sqrt(D[i][j]);
            }
        }
    }
    DA = matrix_multiplication(D_minus_sqrt, A, n, n, n);
    ret = matrix_multiplication(DA, D_minus_sqrt, n, n, n);
    free_matrix(A,n);
    free_matrix(D,n);
    free_matrix(D_minus_sqrt,n);
    free_matrix(DA,n);
    return ret;
}

double **update_H(double **current_H, double **W, int n, int k){
    int i;
    int j;
    double **WH;
    double **H_T;
    double **H_H_T;
    double **H_H_T_H;
    double **ret;
    WH = matrix_multiplication(W, current_H, n, n, k);
    H_T = matrix_transpose(current_H, n, k);
    H_H_T = matrix_multiplication(current_H, H_T, n, k, n);
    H_H_T_H = matrix_multiplication(H_H_T, current_H, n, n, k);
    ret = malloc(n * sizeof(double*));
    assert(ret!=NULL);
    for(i = 0; i < n; i++){
        ret[i] = malloc(k * sizeof(double));
        assert(ret[i]!=NULL);
        for(j = 0; j < k; j++){
            ret[i][j] = current_H[i][j] * (1 - BETA + (BETA * (WH[i][j]/H_H_T_H[i][j])));
        }
    }
    free_matrix(WH, n);
    free_matrix(H_T, k);
    free_matrix(H_H_T, n);
    free_matrix(H_H_T_H, n);
    return ret;
}

double frob_distance(double **A, double **B, int n, int k){
    int i;
    int j;
    double ret;
    double **distance_matrix = malloc(n * sizeof(double*));
    assert(distance_matrix!=NULL);
    for(i = 0; i<n; i++){
        distance_matrix[i] = malloc(k * sizeof(double));
        assert(distance_matrix[i]!=NULL);
        for(j = 0; j < k; j++){
            distance_matrix[i][j] = A[i][j] - B[i][j];
        }
    }
    ret = frobenius_norm_squared(distance_matrix, n, k);
    free_matrix(distance_matrix, n);
    return ret;
}

double **symnmf(double **current_H, double **W, int n, int k) {
    int i;
    double **next_H;
    for(i = 0; i < MAX_ITER; i++){
        next_H = update_H(current_H, W, n, k);
        if(frob_distance(next_H, current_H, n, k) < EPS)
        {
            break;
        }
        free_matrix(current_H, n);
        current_H = next_H;
    }
    if(i < MAX_ITER)
        free_matrix(current_H, n);
    return next_H;
}

int main(int argc, char *argv[]) {

    int n = 0;
    int d = 1;
    int i;
    int j;
    double **X;
    char ch;
    int error_occured = 0;
    char *goal;
    char *filename;
    FILE *file;
    double **result;

    if (argc < 3) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    goal = argv[1];
    filename = argv[2];

    file = fopen(filename, "r");
    if (file == NULL) {
        perror("An Error Has Occurred\n");
        return 1;
    }
    
    while (!feof(file)) {
        ch = fgetc(file);
        if (ch == '\n') {
            n++;
        } else if (n == 0 && ch == ',') {
            d++;
        }
    }
    fseek(file, 0L, SEEK_SET);
    X = (double **)malloc(n * sizeof(double *));
    for (i = 0; i < n; i++) {
        X[i] = (double *)malloc(d * sizeof(double));
        for (j = 0; j < d; j++) {
            if (fscanf(file, "%lf,", &X[i][j]) != 1) {
                printf("An Error Has Occurred\n");
                fclose(file);
                return 1;
            }
        }
    }
    fclose(file);

    if (strcmp(goal, "sym") == 0) {
        result = sym(X, n, d);
    } else if (strcmp(goal, "ddg") == 0) {
        result = ddg(X, n, d);
    } else if (strcmp(goal, "norm") == 0) {
        result = norm(X, n, d);
    } else {
        printf("An Error Has Occurred!\n");
        exit(1);
    }
    
    for(i = 0; i < n; i++){
        print_vector(result[i], n);
    }

    for (i = 0; i < n; i++) {
        free(X[i]);
    }
    free(X);
    for (i = 0; i < n; i++) {
        free(result[i]);
    }
    free(result);
    
    return error_occured;
}
