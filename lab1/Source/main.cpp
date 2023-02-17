#include <iostream>
#pragma warning (disable:4996)
#define PI 3.1415926
#define ABS(x) ((x)>0.0?(x):-(x))
int i, j;

double s;

void Task1() {
  double eps = 1.0;
  while ((eps /= 2.0) + 1.0 > 1.0);
  printf("eps = %g", eps);
}

double mathfunc1(double x) {
  return exp(-x) / x;
}

double Sympson(double a, double b, int n, double(*mathfunc)(double)) {
  //sympson
  double sum, sum1, sum2, sum3, sum4;
  double h, x;

  sum1 = 0.0;
  sum2 = 0.0;
  sum3 = 0.0;
  sum4 = 0.0;

  sum = 0.0;
  //n *= 2;
  h = (b - a) / (2.0 * n);

  // = 0 as sin(x) = 0
  sum1 = (mathfunc)(a);
  sum2 = (mathfunc)(b);
  //sum1 = 0.0;
  //sum2 = 0.0;

  x = a + h;
  for (int k = 1; k <= n; k++) {
    sum3 += (mathfunc)(x);
    x += (2.0 * h);
  }
  sum3 *= 4.0;

  x = a + 2.0 * h;
  for (int k = 1; k < n; k++) {
    sum4 += (mathfunc)(x);
    x += (2.0 * h);
  }
  sum4 *= 2.0;

  sum = sum1 + sum2 + sum3 + sum4;
  sum *= (h / 3.0);
  return sum;
}

void Task2() {

  double eps = 1.0;

  double a = eps, b = 10.0;

  int n = 10;

  //sympson
  double sum;

  printf("n =\t");
  for (n = 2; n < 1e8; n *= 4) {
    printf("%i\t", n);
  }
  printf("\n");

  while (eps > 1e-6) {
    eps *= 0.1;
    printf("eps = %g | ", eps);
    a = eps;
    n = 2;
    while (n < 1e8) {

      sum = Sympson(a, b, n, mathfunc1);
      printf("%.10lf ",sum);
      n *= 4;
    }
    printf("\n");
  }

}


double mathfunc2(double x) {
  return 2.0 / PI * (sin(i * x) * sin(j * x) + ((i + j * j) * cos((i + j * j) * x) + 2.0 * sin((i + j * j) * x)) * exp(2.0 * x));
}

double ak1(int k) {
  return 1.0 / pow((2.0 * k + 1), double(i + 1.0) / i);
}

double ak2(int k) {
  return 1.0 / pow(k, double(i + 1.0) / i);
}

double ak3(int k) {
  return 1.0 / pow(2.0, k);
}

double ak(int k) {
  return 1.0 / pow(k, s);
}

double Sum(double(*ak)(int), int n2) {
  double sum = 0.0;
  int k = 0;
  do {
    ++k;
    sum += ak(k);
  } while (k <= n2);
  return sum;
}

double VectorNorm(double* v1, double* v2, int N) {
  double max = 0.0, temp;
  for (int k = 0; k < N; ++k) {
    temp = ABS(v1[k] - v2[k]);
    if (max < temp)
      max = temp;
  }
  return max;
}

double Adams1(double a, double b, unsigned long long n, double y0, double z0);
double FixedEulerMethod(double a, double b, unsigned long long n, double y0, double z0);

void Task3() {
  //int N = 10, n1[] = { 10000, 100000, 1000000 }, n2[] = { 100, 100000, 1000000 }, n4[] = { 100, 100000, 1000000 };
  //double eps3[] = {1e-1, 1e-2, 1e-3};

  int N = 10, n1[] = {100000 }, n2[] = {100000 }, n4[] = {1000000 };
  double eps3[] = {1e-2};

  const double a = 0.0, b = PI;
  double** B = nullptr;
  double* f = nullptr;

  //FILE* file = fopen("file.txt", "w");
  FILE* file = stdin;

  double ans, * ai = nullptr, * ai2 = nullptr;

  B = new double* [N];
  for (int k = 0; k < N; ++k)
    B[k] = new double[N];

  f = new double[N];

  int SIZE = sizeof(n1)/sizeof(n1[0]);

  for (int i1 = 0; i1 < SIZE; ++i1) {
    fprintf(file, "n1 = %i\n", n1[i1]);

    for (int i2 = 0; i2 < SIZE; ++i2) {
      fprintf(file, "n2 = %i\n", n2[i2]);

      for (int i3 = 0; i3 < SIZE; ++i3) {
        fprintf(file, "eps3 = %g\n", eps3[i3]);

        for (int i4 = 0; i4 < SIZE; ++i4) {
          fprintf(file, "n4 = %i ", n4[i4]);

          //Матрица системы
          for (i = 1; i <= N; ++i) {
            for (j = 1; j <= N; ++j)
              B[i - 1][j - 1] = Sympson(a, b, n1[i1], mathfunc2);
          }

          
          printf("\n\n\n");
          for (i = 0; i < N; ++i) {
            for (j = 0; j < N; ++j) {
              printf("%lf ", B[i][j]);
            }
            printf("\n");
          }
          printf("\n\n\n");
          

          //правая часть
          double s;
          for (i = 1; i <= N; ++i) {
            s = double(i + 1.0) / (i);
            f[i - 1] = sqrt(pow(2.0, s) / (pow(2.0, s) - 1.0) * (Sum(ak1, n2[i2]) + 1.0) / Sum(ak2, n2[i2]));
          }

          
          //Решаем СЛАУ
          printf("\n\n");
          for (i = 0; i < N; ++i)
            printf("%lf ", f[i]);
          printf("\n\n");
          

          double max = 0.0;

          //find ||B|| = q
          double q = 0.0;
          for (i = 0; i < N; ++i) {
            for (j = 0; j < N; ++j) {

              if (i != j) {
                max = ABS(B[i][j] / B[i][i]);
                if (q < max)
                  q = max;
              }

            }
          }

          //Решаем СЛАУ методом Якоби
          ai = new double[N];
          for (i = 0; i < N; ++i)
            ai[i] = 135.0; //any temp value, more than ai2[i]

          ai2 = new double[N];
          for (i = 0; i < N; ++i)
            ai2[i] = 0.0;

          //Метода Якоби
          double temp;
          do {

            //copy
            for (i = 0; i < N; ++i)
              ai2[i] = ai[i];

            for (i = 0; i < N; ++i) {
              temp = 0.0;
              for (j = 0; j < N; ++j) {
                if (i != j)
                  temp += B[i][j] * ai2[j];
              }
              temp -= f[i];
              ai[i] = temp / (-B[i][i]);
            }
            
            printf("\n");
            for (i = 0; i < N; ++i)
              printf("%lf ", ai[i]);
            

          } while (VectorNorm(ai, ai2, N) >= (1.0 - q) / q * eps3[i3]);


          //Задача Коши
          //z'' - 9z' - 10z = 0, 0<=x<=8
          //z(0) = 1, z'(0) = max |ai|
          double z0 = 1.0, z_0 = 0.0;
          for (int k = 0; k < N; ++k) {
            temp = ABS(ai[k]);
            if (max < temp)
              max = temp;
          }
          max = -max;
          //max = -1.0;
          printf("z'(0)=%lf\n", max);
          ans = Adams1(0.0, 8.0, n4[i4], 1.0, max);
          //ans = FixedEulerMethod(0.0, 8.0, n, 1.0, max);
          // 
          printf("\nz(8) = %g\n\n", ans);

        }
      }
    }
  }


  for (int k = 0; k < N; ++k)
    delete[] B[k];
  delete[] B;
  delete[] f;

  delete[] ai;
  delete[] ai2;
}

int main() {
  //Task1();
  //Task2();
  Task3();
  return 0;
}