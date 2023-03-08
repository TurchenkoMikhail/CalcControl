#include <iostream>
#include <stdlib.h>
#include <cmath>

double** CreateUnitMatrix(int N) {
  int i, j;
  double** matrix = new double* [N];

  for (i = 0; i < N; ++i)
    matrix[i] = new double[N];
  
  for (j = 0; j < N; ++j) {
    for (i = 0; i < N; ++i)
      matrix[j][i] = (i == j ? 1.0 : 0.0);
  }
  return matrix;
}

void PrintMatrix(double** matrix, int N) {
  int i, j;
  for (j = 0; j < N; ++j) {
    for (i = 0; i < N; ++i)
      printf("%g ", matrix[j][i]);

    printf("\n");
  }
}

void PrintVector(double* v, int N) {
  for (int j = 0; j < N; ++j) 
    printf("%g ", v[j]);
  printf("\n");
}

void DestroyMatrix(double** matrix, int N) {
  int i;
  for (i = 0; i < N; ++i)
    delete[] matrix[i];
  delete[] matrix;
}

double VectorNorm(double* vector, int N) {
  int i;
  double max = 0.0;
  for (i = 0; i < N; ++i) {
    if (max < fabs(vector[i]))
      max = fabs(vector[i]);
  }
  return max;
}

double* VectorDiff(double* v1, double* v2, int N) {
  int i;
  double* vec = new double[N];
  for (int i = 0; i < N; ++i)
    vec[i] = v1[i] - v2[i];
  return vec;
}

double* MatrixMulVector(double** A, double* x, int N) {
  int j, k;
  double temp;
  double* vec = new double[N];

  for (j = 0; j < N; ++j) {
    temp = 0.0;
    for (k = 0; k < N; ++k)
      temp += (A[j][k] * x[k]);
    vec[j] = temp;
  }
  return vec;
}

double** MatrixSum(double** A, double** B, int N) {
  int i, j;
  double** matrix = new double* [N];
  for (i = 0; i < N; ++i)
    matrix[i] = new double[N];

  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) 
      matrix[i][j] = A[i][j] + B[i][j];
  }

  return matrix;
}

double** MatrixDiff(double** A, double** B, int N) {
  int i, j;
  double** matrix = new double* [N];
  for (i = 0; i < N; ++i)
    matrix[i] = new double[N];

  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      matrix[i][j] = A[i][j] - B[i][j];
    }
  }

  return matrix;
}

double* VectorSum(double* A, double* B, int N) {
  int i;
  double* ans = new double[N];

  for (i = 0; i < N; ++i)
    ans[i] = A[i] + B[i];
  return ans;
}

double MatrixNorm(double** A, int N) {
  int i,j;
  double ans = 0.0;
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      if (ans < fabs(A[i][j]))
        ans = fabs(A[i][j]);
    }
  }
  return ans;
}

double** MatrixMulMatrix(double** A, double** B, int N) {
  int i, j, k;
  double** Ñ = new double* [N];
  for (i = 0; i < N; ++i)
    Ñ[i] = new double[N];

  for (j = 0; j < N; ++j) {
    for (i = 0; i < N; ++i)
      Ñ[j][i] = 0.0;
  }

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < N; k++)
        Ñ[i][j] += A[i][k] * B[k][j];
    }
  }
  return Ñ;
}