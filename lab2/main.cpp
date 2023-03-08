#include <iostream>
#include <cmath>

#define ABS(x) ((x)>0.0?(x):-(x))

double f(double x, double a, double p) {
  return (1.0 - a) * x - pow(x, -p);
}

double df(double x, double a, double p) {
  return (1.0 - a) + p * pow(x, -(p + 1.0));
}

void NewtonMethod(double x0, double eps, double a, double p, double q, double ans) {

  double xkprev = x0, xk = x0, xknext = x0;
  int k = 0;

  const int N = 50;
  double Mk_up[N], Mkl_up[N], ek[N], \
    Mkl_down[N], Mk_down[N], x[N], Mk0_up[N];

  double x1 = x0 - f(x0, a, p) / df(x0, a, p);
  double diff0 = ABS(x1 - x0);

  const int L = 2; //for advanced majorants and minorants
  double xkl, xk_temp;

  do {
    ++k;

    xkprev = xk;
    xk = xknext;

    //advanced
    xk_temp = xk;
    for (int l = 1; l <= L; ++l) {
      xkl = xk_temp - f(xk_temp, a, p) / df(xk_temp, a, p);
      xk_temp = xkl;
    }
    Mkl_down[k - 1] = 1.0 / (1.0 + pow(q, L)) * ABS(xk - xkl);
    Mkl_up[k - 1] = 1.0 / (1.0 - pow(q, L)) * ABS(xk - xkl);

    //newton method
    xknext = xk - f(xk, a, p) / df(xk, a, p);

    x[k - 1] = xk;

    //count majorants and minorants
    Mk0_up[k - 1] = pow(q, k) / (1.0 - q) * diff0;
    Mk_up[k - 1] = q / (1.0 - q) * ABS(xk - xkprev);
    ek[k - 1] = ABS(xk - ans);
    Mk_down[k - 1] = 1.0 / (1.0 + q) * ABS(xknext - xk);

  } while (ABS(xknext - xk) > eps);

  printf("k: ");
  for (int i = 1; i <= k; ++i)
    printf("%i ", i);
  printf("\n");
  printf("xk: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", x[i]);
  printf("\n");
  printf("Mk0_up: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", Mk0_up[i]);
  printf("\n");
  printf("Mk_up: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", Mk_up[i]);
  printf("\n");
  printf("Mkl_up: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", Mkl_up[i]);
  printf("\n");
  printf("ek: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", ek[i]);
  printf("\n");
  printf("Mkl_down: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", Mkl_down[i]);
  printf("\n");
  printf("Mk_down: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", Mk_down[i]);
  printf("\n");

}

void Task1() {

  const int SIZE = 3;
  const double x0[SIZE] = { 1.2, 1.23, 1.26 };
  const double p[SIZE] = { 1.0, 2.0, 4.0 };
  const double ans[SIZE] = { sqrt(10.0 / 3.0), pow(10.0 / 3.0, 1.0 / 3.0), pow(10.0 / 3.0, 1.0 / 5.0) };
  const double q[SIZE] = { 0.55, 0.45, 0.91 };

  //debug
  /*
  const int SIZE = 1;
  const double x0[SIZE] = { 1.2 };
  const double p[SIZE] = { 1.0 };
  const double ans[SIZE] = { sqrt(10.0 / 3.0) };
  const double q[SIZE] = { 0.55 };
  */

  double a = 0.7;

  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {

      printf("x0 = %g, p = %g\n", x0[i], p[j]);
      NewtonMethod(x0[i], 1e-9, a, p[j], q[j], ans[j]);
      printf("\n\n\n");
    }
  }

}

double** CreateUnitMatrix(int N);
void PrintMatrix(double** matrix, int N);
void PrintVector(double* v, int N);
double** MatrixMulMatrix(double** A, double** B, int N);
double** MatrixDiff(double** A, double** B, int N);
void DestroyMatrix(double** matrix, int N);
double* MatrixMulVector(double** A, double* x, int N);
double* VectorSum(double* v1, double* v2, int N);
double VectorNorm(double* vector, int N);
double* VectorDiff(double* v1, double* v2, int N);
double MatrixNorm(double** A, int N);

//residual vector
double* R(double* z, double** L, double* b, int N) {

  //ans = Lz - z + b

  double* temp1 = MatrixMulVector(L, z, N);
  double* temp2 = VectorDiff(temp1, z, N);
  double* ans = VectorSum(temp2, b, N);
  delete[] temp1;
  delete[] temp2;
  return ans;
}

void Task2() {
  double eps = 1e-5;
  int k = 0;
  const int SIZE = 500, N = 10;
  double M_up[SIZE], M_down[SIZE], M0_up[SIZE], ek[SIZE], Mkl_up[SIZE], Mkl_down[SIZE];

  double q;

  //fill matrix and right-part vector
  double** A = CreateUnitMatrix(N);
  double* f;
  double* ans = new double[N];

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {

      if (i != j)
        A[i][j] = 0.93 / (i + j + 2);
      else
        A[i][j] = i + 1.0;
    }
  }

  for (int i = 0; i < N; ++i)
    ans[i] = 1.0 + 1.0 / (i + 1.0);

  f = MatrixMulVector(A, ans, N);


  printf("A:\n");
  PrintMatrix(A, N);
  printf("\n\n");
  printf("f:\n");
  PrintVector(f, N);
  printf("\n\n");

  //create iter matrix
  double* b;
  double** L;
  double** temp, ** unit, ** B_;

  unit = CreateUnitMatrix(N);
  B_ = CreateUnitMatrix(N);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (i == j)
        B_[i][j] = 1.0 / A[i][j];
      else
        B_[i][j] = 0.0;
    }
  }

  /*
  printf("B_:\n");
  PrintMatrix(B_, N);
  printf("\n\n");
  */

  temp = MatrixMulMatrix(B_, A, N);
  L = MatrixDiff(unit, temp, N);
  b = MatrixMulVector(B_, f, N);

  //q = MatrixNorm(L, N);
  q = 0.5;

  DestroyMatrix(temp, N);
  DestroyMatrix(unit, N);
  DestroyMatrix(B_, N);


  printf("L:\n");
  PrintMatrix(L, N);
  printf("\n\n");
  printf("b:\n");
  PrintVector(b, N);
  printf("\n\n");


  double* x0 = new double[N];
  for (int i = 0; i < N; ++i)
    x0[i] = 0.0;

  double* temp1;
  double* xk = new double[N], * xknext = new double[N];
  double* residualprev, * residual, norm_first_vector, * error;

  residual = R(x0, L, b, N);
  norm_first_vector = VectorNorm(residual, N);
  delete[] residual;

  for (int i = 0; i < N; ++i) {
    xk[i] = x0[i];
    xknext[i] = x0[i];
  }

  k = 0;
  do {
    ++k;
    xk = xknext;

    temp1 = MatrixMulVector(L, xk, N);
    xknext = VectorSum(temp1, b, N);
    delete[] temp1; temp1 = nullptr;

    //xk[0] = 2.0; xk[1] = 1.5;
    residualprev = R(xk, L, b, N);
    residual = R(xknext, L, b, N);
    error = VectorDiff(xknext, ans, N);

    //fill majorants and minorants
    M0_up[k - 1] = pow(q, k) / (1.0 - q) * norm_first_vector;
    M_up[k - 1] = q / (1.0 - q) * VectorNorm(residualprev, N);
    ek[k - 1] = VectorNorm(error, N); printf("%g %g\n", xknext[0], xknext[1]);
    M_down[k - 1] = 1.0 / (1.0 + q) * VectorNorm(residual, N);

    delete[] residualprev; residualprev = nullptr;
    delete[] residual; residual = nullptr;
    delete[] error; error = nullptr;
    delete[] xk; xk = nullptr;

  } while (M_up[k - 1] > eps || k == 1);

  printf("\n\nk: ");
  for (int i = 1; i <= k; ++i)
    printf("%i ", i);
  printf("\n");
  printf("M0_up: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", M0_up[i]);
  printf("\n");
  printf("M_up: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", M_up[i]);
  printf("\n");
  printf("ek: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", ek[i]);
  printf("\n");
  printf("M_down: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", M_down[i]);
  printf("\n");

  delete[] xk;
  delete[] xknext;
}

double T(double* ujprev, double a, double b, int n, int i) {
  //return 100.0 * (1.0 - cos(t)); - 75.0 * ujprev / 76.0;

  double h = (b - a) / n;
  double ti = h * i;
  double sum1 = 100.0 * (1.0 - cos(ti));
  double sum2 = 0.0;
  
  for (int k = 0; k < i; ++k)
    sum2 += (ujprev[k] + ujprev[k+1])/2.0;
  
  return sum1 - 75.0 * h * sum2;
}

double ExactSol(double t) {
  return 50.0 / 2813.0 * (exp(-75.0 * t) + 75.0 * sin(t) - cos(t));
}

double MaxNorm(double* u1, double* u2, int N) {
  double max = 0.0;
  for (int i = 0; i <= N; ++i) {
    if (max < ABS(u1[i] - u2[i]))
      max = ABS(u1[i] - u2[i]);
  }
  return max;
}

void Task3() {

  const double t0 = 0.0, t1 = 1.0 / 76.0, L = 75.0;
  const double q = L * (t1 - t0);

  double u0 = 0.0;
  int k = 0;
  const double eps = 1e-6;

  const int SIZE = 100, n = 50;
  double M_up[SIZE + 1], M_down[SIZE + 1], M0_up[SIZE + 1], ek[SIZE + 1], Mkl_up[SIZE+1], Mkl_down[SIZE+1];

  double uk[n + 1], ukprev[n + 1], uknext[n + 1], ans[n + 1];
  double ukl[n + 1], ukltemp[n + 1];

  double h = (t1 - t0) / n;
  double t = t0;

  double first_diff;

  //fill 
  t = t0;
  for (int i = 0; i <= n; ++i, t += h) {
    ans[i] = ExactSol(t);

    ukprev[i] = ExactSol(t0);
    uk[i] = ExactSol(t0);
    uknext[i] = ExactSol(t0);
  }

  k = 0;
  do {
    ++k;

    //uk = u0, ukprev = u0, uknext = u0;
    for (int i = 0; i <= n; ++i) {
      ukprev[i] = uk[i];
      uk[i] = uknext[i];
    }

    for (int i = 0; i <= n; ++i, t += h) 
      uknext[i] = T(uk, t0, t1, n, i);

    if (k == 1)
      first_diff = MaxNorm(uknext, uk, n);

    //advanced
    for (int i = 0; i <= n; ++i) {
      ukl[i] = uk[i];
      ukltemp[i] = uk[i];
    }
    for (int l = 1; l <= 2; ++l) {
      for (int i = 0; i <= n; ++i) 
        ukl[i] = T(ukltemp, t0, t1, n, i);

      for (int i = 0; i <= n; ++i)
        ukltemp[i] = ukl[i];
    }

    Mkl_down[k - 1] = 1.0 / (1.0 + pow(q, L)) * MaxNorm(uk, ukl, n);
    Mkl_up[k - 1] = 1.0 / (1.0 - pow(q, L)) * MaxNorm(uk, ukl, n);

    M0_up[k - 1] = pow(q, k) / (1.0 - q) * first_diff;
    M_up[k - 1] = q / (1.0 - q) * MaxNorm(uk, ukprev, n); //ABS(uk - ukprev);
    ek[k - 1] = MaxNorm(uk, ans, n); //ABS(uk - ExactSol(t));
    M_down[k - 1] = 1.0 / (1.0 + q) * MaxNorm(uknext, uk, n);//ABS(uknext - uk);

  } while (ABS(ek[k - 1]) > eps || k == 1);

  printf("\n\nk: ");
  for (int i = 1; i <= k; ++i)
    printf("%i ", i);
  printf("\n");
  printf("M0_up: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", M0_up[i]);
  printf("\n");
  printf("M_up: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", M_up[i]);
  printf("\n");
  printf("Mkl_up: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", Mkl_up[i]);
  printf("\n");
  printf("ek: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", ek[i]);
  printf("\n");
  printf("Mkl_down: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", Mkl_down[i]);
  printf("\n");
  printf("M_down: ");
  for (int i = 0; i < k; ++i)
    printf("%g ", M_down[i]);
  printf("\n");
}

int main() {
  //Task1();
  Task2();
  //Task3();
  return 0;
}