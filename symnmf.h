#ifndef SYMNMFMODULE_H
#define SYMNMFMODULE_H

double **sym(double** X, int n, int d);
double **ddg(double** X, int n, int d);
double **norm(double** X, int n, int d);
double **symnmf(double **H, double **W, int n, int k);

#endif
