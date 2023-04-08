#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cmath>

#define ABS(x) ((x)>0?(x):-(x))
#define MAX(x,y) ((x)>(y)?(x):(y)
#define PI 3.1415926

typedef enum norm_t {
  GLOBAL,
  LOCAL
}norm_t;

//кусочно-линейные базисные функции
class Equation {
public:
  double* xk; //сетка
  int n; //размер сетки
  double h;
  double a, b;
  double* ak; //ответ
  //phi(x) = sum from i=0 to n ak*phi_k(x)

  Equation(int size) {
    this->a = 0.0;
    this->b = 1.0;
    n = size;
    xk = new double[n + 1];
    h = (b - a) / n;
    for (int i = 0; i <= n; ++i)
      xk[i] = a + h * i;
    ak = nullptr;
  }

  ~Equation() {
    delete[] xk;
    delete[] ak;
  }

  //коэффициенты СЛАУ метода Галеркина-Ритца

   //a_ki = integral from 0 to 1 (a(x) * phi_k(x)' * phi_i(x)') dx
  double A(int k, int i) {

    if (k == 0 || i == 0 || k == n || i == n)
      return 0.0;

    double a = xk[k - 1], b = xk[k], c = xk[k + 1];

    if (i == k - 1) {
      return (a + b + 2.0) / (2.0 * (a - b));
    }
    else if (i == k) {
      return (b + 1.0) / (c - b) - (b + 1.0) / (a - b);
    }
    else if (i == k + 1) {
      return (b + c + 2.0) / (2.0 * (b - c));
    }
    else
      return 0.0;

  }

  //b_k = integral from 0 to 1 (f * phi_k(x)) dx
  double B(int k) {

    if (k == 0 || k == n)
      return 0.0;

    double a = xk[k - 1], b = xk[k], c = xk[k + 1];
    double a2 = a * a, b2 = b * b, c2 = c * c;
    double a3 = a * a * a, b3 = b * b * b, c3 = c * c * c;

    //f(x) = 27x^2+6x-6
    return 9.0 / 4.0 * (-a3 - a2 * b - a * b2 + b2 * c + b * c2 + c3) - a2 - a * b + 3.0 * a + b * c + c2 - 3.0 * c;
  }

  double f(double x) {
    return 27.0 * x * x + 6.0 * x - 6.0;
  }

  double phi(int i, double x) {
    if (i == 0 || i == n)
      return 0.0;

    double a = xk[i - 1], b = xk[i], c = xk[i + 1];

    if (a <= x && x <= b)
      return (x - a) / (b - a);
    else if (b < x && x <= c)
      return (c - x) / (c - b);
    else
      return 0.0;
  }

  double phidiff(int i, double x) {
    if (i == 0 || i == n)
      return 0.0;

    double a = xk[i - 1], b = xk[i], c = xk[i + 1];

    if (a <= x && x <= b)
      return 1.0 / (b - a);
    else if (b < x && x <= c)
      return -1.0 / (c - b);
    else
      return 0.0;
  }

  //numeric solution
  double uh(double x) {
    double ans = 0.0;
    for (int i = 1; i <= n - 1; ++i)
      ans += (ak[i] * phi(i, x));
    return ans;
  }

  //exact solution
  double u(double x) {
    return 3.0 * x * x * (1.0 - x);
  }

  double udiff(double x) {
    return 3.0 * x * (2.0 - 3.0 * x);
  }

  double uhdiff(double x) {
    double ans = 0.0;
    for (int i = 1; i <= n - 1; ++i)
      ans += (ak[i] * phidiff(i, x));
    return ans;
  }

  double uh1diff(double x) {
    double ans = 0.0, a;
    for (int i = 1; i <= n - 1; ++i)
      ans += ((ak[i] * 1.05) * phidiff(i, x));
    return ans;
  }

  double uh2diff(double x) {
    return 0.25 * PI * cos(PI * x);
  }

  double Integral_diffphi_i_diffphi_k(int i, int k, norm_t type, int j) {
    if (k == 0 || i == 0 || k == n || i == n)
      return 0.0;

    if (type == GLOBAL) {
      double a = xk[k - 1], b = xk[k], c = xk[k + 1];

      if (i == k - 1) {
        return -1.0 / (b - a);
      }
      else if (i == k) {
        return 1.0 / (b - a) + 1.0 / (c - b);
      }
      else if (i == k + 1) {
        return -1.0 / (c - b);
      }
      else
        return 0.0;

    }

    else if (type == LOCAL) {
      double b = xk[j], c = xk[j + 1];

      if (i == k - 1 && j == k - 1)
        return -1.0 / (c - b);

      else if (i == k) {

        if (j == k - 1)
          return 1.0 / (c - b);
        else if (j == k)
          return 1.0 / (c - b);
        else
          return 0.0;
      }

      else if (i == k + 1 && j == k)
        return -1.0 / (c - b);

      else
        return 0.0;
    }

  }

  double Integral_diffu_diffphi_k(int k, norm_t type, int j) {

    //[0,1]
    if (type == GLOBAL) {
      double a = xk[k - 1], b = xk[k], c = xk[k + 1];
      return -3.0 * (a - c) * (a + b + c - 1.0);
    }

    //[xi, xi+1]
    else if (type == LOCAL) {
      double b = xk[j], c = xk[j + 1];

      if (j == k - 1)
        return -3.0 * (b * b + b * (c - 1.0) + (c - 1.0) * c);
      else if (j == k)
        return 3.0 * (b * b + b * (c - 1.0) + (c - 1.0) * c);
      else return 0.0;

    }

  }

  double Integral_diffu2(norm_t type, int j) {
    if (type == GLOBAL) { //[0,1]
      return 6.0 / 5.0;
    }
    else if (type == LOCAL)
    {
      double b = xk[j], c = xk[j + 1];

      double b3 = b * b * b, c3 = c * c * c;
      double b4 = b3 * b, c4 = c3 * c;
      double b5 = b4 * b, c5 = c4 * c;

      return 3.0 / 5.0 * (27.0 * (c5 - b5) + 45.0 * (b4 - c4) + 20.0 * (c3 - b3));
    }

  }


  /////[0,1]

  //f^2 dx on [xi, xi+1]
  double Integral_f2(int i) {
    double b = xk[i], c = xk[i + 1];

    double b2 = b * b, c2 = c * c;
    double b3 = b2 * b, c3 = c2 * c;
    double b4 = b3 * b, c4 = c3 * c;
    double b5 = b4 * b, c5 = c4 * c;

    return 729.0 / 5.0 * (c5 - b5) + 81.0 * (c4 - b4) + 96.0 * (b3 - c3) + 36.0 * (b2 - c2) + 36.0 * (c - b);

  }

  //f*diffphi_k(x) on [xi, xi+1]
  double Integral_f_diffphi_k(int k, int i) {
    double b = xk[i], c = xk[i + 1];

    if (i == k - 1)
      return 3.0 * (3.0 * b * b + 3.0 * b * c + b + 3.0 * c * c + c - 2.0);
    else if (i == k)
      return -3.0 * (3.0 * b * b + 3.0 * b * c + b + 3.0 * c * c + c - 2.0);
    else return 0.0;
  }

  //diffphi_k(x)*diffphi_j(x) on [xi, xi+1]
  double Integral_diffphi_k_diffphi_j(int k, int j, int i) {

    double b = xk[i], c = xk[i + 1];

    if (k == j - 1 && i == j - 1)
      return 1.0 / (b - c);

    else if (k == j) {

      if (i == j - 1)
        return 1.0 / (c - b);
      else if (i == j)
        return 1.0 / (c - b);
      else return 0.0;

    }

    else if (k == j + 1 && i == j)
      return 1.0 / (b - c);

    else return 0.0;
  }

  double TrueErrorGalerkin(norm_t type, int j) {
    double ans = 0.0;
    double s1 = 0.0, s2 = 0.0, s3 = 0.0;

    s1 = Integral_diffu2(type, j);

    for (int k = 1; k <= n - 1; ++k)
      s2 += (ak[k] * Integral_diffu_diffphi_k(k, type, j));
    s2 *= (-2.0);

    for (int i = 1; i <= n - 1; ++i) {
      for (int k = i - 1; k <= i + 1; ++k)
        s3 += (ak[i] * ak[k] * Integral_diffphi_i_diffphi_k(i, k, type, j));
    }

    ans = s1 + s2 + s3;
    return sqrt(ans);
  }

  double EstimateErrorGalerkin(int i) {
    double b = xk[i], c = xk[i + 1];
    double ci = (c - b) / PI;
    double s1 = 0.0, s2 = 0.0, s3 = 0.0;

    s1 = Integral_f2(i);

    for (int k = 1; k <= n - 1; ++k)
      s2 += (ak[k] * Integral_f_diffphi_k(k, i));
    s2 *= (2.0);

    for (int j = 1; j <= n - 1; ++j) {
      for (int k = j - 1; k <= j + 1; ++k)
        s3 += (ak[j] * ak[k] * Integral_diffphi_k_diffphi_j(k, j, i));
    }

    return ci * sqrt(s1 + s2 + s3);
  }

  double EstimateErrorGalerkinGlobal() {
    double ri, ans = 0.0;
    for (int i = 0; i < n; ++i) {
      ri = EstimateErrorGalerkin(i);
      ans += (ri * ri);
    }
    return sqrt(ans);
  }

  //uh = 0.25*sin(PI*x)
  double TrueErrorNonGalerkin(norm_t type, int j) {
    if (type == GLOBAL)
      return sqrt(6.0/5.0 - 3.0/PI + PI*PI / 32.0);

    else if (type == LOCAL) {
      double b = xk[j], c = xk[j + 1];
      double s1 = 0.0, s2 = 0.0, s3 = 0.0;

      s1 = Integral_diffu2(type, j);
      s2 = 3.0 * ((PI * PI * b * (3.0 * b - 2.0) - 6.0) * sin(PI * b) + 2.0 * PI * (3.0 * b - 1.0) * cos(PI * b) + \
        (PI * PI * (2.0 - 3.0 * c) * c + 6.0) * sin(PI * c) + 2.0 * PI * (1.0 - 3.0 * c) * cos(PI * c)) / (4.0 * PI * PI);
      s3 = PI / 64.0 * (2.0 * PI * (c - b) - sin(2.0 * PI * b) + sin(2.0 * PI * c));

      return sqrt(s1 + s2 + s3);
    }

  }

  //uh = 0.25*sin(PI*x)
  double EstimateErrorNonGalerkin(int i) {
    double b = xk[i], c = xk[i + 1];
    double ci = (c - b) / PI;

    double s1 = 0.0, s2 = 0.0, s3 = 0.0;
    s1 = Integral_f2(i);
    s2 = 3.0 * ((PI * PI * (2.0 - b * (9.0 * b + 2.0)) + 18.0) * sin(PI * b) - 2.0 * PI * (9.0 * b + 1.0) * cos(PI * b) + \
      (PI * PI * (c * (9.0 * c + 2.0) - 2.0) - 18.0) * sin(PI * c) + 2.0 * PI * (9.0 * c + 1.0) * cos(PI * c)) / (4.0 * PI * PI);
    s3 = PI / 64.0 * (2.0 * PI * (c - b) - sin(2.0 * PI * b) + sin(2.0 * PI * c));
    return ci * sqrt(s1 + s2 + s3);
  }

  //uh = 0.25*sin(PI*x)
  double EstimateErrorNonGalerkinGlobal() {
    double ri, ans = 0.0;
    for (int i = 0; i < n; ++i) {
      ri = EstimateErrorNonGalerkin(i);
      ans += (ri * ri);
    }
    return sqrt(ans);
  }

};

void GalerkinRitzMethod(Equation* eq) {
  double h = (eq->b - eq->a) / eq->n;
  int n = eq->n;
  eq->ak = new double[n + 1];
  // -u'' = f, u(a) = u(b) = 0
  //u_ = sum y_i*phi_i(x), i=1,n
  //phi_i(x) - базисные функции
  //СЛАУ вида sum a_ki * y_i = b_k
  //a_ki = integral from 0 to 1 |phi_k(x)' * phi_i(x)'|dx
  //b_k = integral from 0 to 1 f*phi_k(x) dx
  // 
  //Метод прогонки
  double* R = new double[n + 1];

  double* B = new double[n + 1];
  double* C = new double[n + 1];
  double* D = new double[n + 1];

  for (int j = 0; j <= n; ++j) {
    D[j] = eq->A(j + 1, j);
    C[j] = eq->A(j, j);
    B[j] = eq->A(j - 1, j);
    R[j] = eq->B(j);
  }
  B[0] = 0.0;
  D[n] = 0.0;

  double* d = new double[n + 1], * l = new double[n + 1];

  //Метод прогонки
  d[1] = -D[1] / C[1];
  l[1] = R[1] / C[1];
  for (int j = 2; j < n; ++j) {
    d[j] = -D[j] / (C[j] + B[j] * d[j - 1]);
    l[j] = (R[j] - B[j] * l[j - 1]) / (C[j] + B[j] * d[j - 1]);
  }
  d[n] = 0.0;

  eq->ak[0] = 0.0;
  eq->ak[n] = 0.0;

  for (int j = n - 1; j > 0; --j) {
    eq->ak[j] = d[j] * eq->ak[j + 1] + l[j];
  }

  delete[] B;
  delete[] C;
  delete[] D;
  delete[] R;
  delete[] d;
  delete[] l;
}



int main() {
  //1)

  /*
  int n[] = { 3,4,5,10,20,50,100,200,500, 1000};
  int SIZE = sizeof(n) / sizeof(n[0]);

  //1)

  
  int NONE = -1;

  printf("||u'-uh'||l2\t|R(uh)| <= \nGalerkin uh\nNonGalerkin 1.05*uh\nnNonGalerkin 0.25sin(Pi*x)\n");

  for (int i = 0; i < SIZE; ++i) {
    printf("\nn = %i\n", n[i]);
    double h = 1.0 / (double)n[i];
    Equation* eq = new Equation(n[i]);
    GalerkinRitzMethod(eq);
    printf("%.6lf\t%.6lf\n", eq->TrueErrorGalerkin(GLOBAL, NONE), eq->EstimateErrorGalerkinGlobal());

    //error
    for (int k = 1; k <= n[i] - 1; ++k)
      eq->ak[k] *= 1.05;
    printf("%.6lf\t%.6lf\n", eq->TrueErrorGalerkin(GLOBAL, NONE), eq->EstimateErrorGalerkinGlobal());

    printf("%.6lf\t%.6lf\n", eq->TrueErrorNonGalerkin(GLOBAL, NONE), eq->EstimateErrorNonGalerkinGlobal());
    delete eq;
  }
  */

  //test galerkin method
  
  int n = 100;
  double h = 1.0 / n;
  Equation* eq = new Equation(100);
  GalerkinRitzMethod(eq);
  printf("numeric\t\t\texact\t\t\terr\n\n");
  double x = 0.0, fx;

  x = -0.001;
  for (int i = 0; i <= n; ++i, x += h)
    printf("%lf ", eq->uh(x));
  printf("\n");
  x = -0.001;
  for (int i = 0; i <= n; ++i, x += h)
    printf("%lf ", eq->u(x));
  printf("\n");

  /*
  for (int i = 0; i <= n; ++i, x += h) {
    fx = eq->u(x);
    printf("f(%lf) = %lf\tu(%lf) = %lf\t%lf\n", x, fx, x, eq->u(x), ABS(fx - eq->u(x)));
  }
  */


  //2)

  /*
  double x, h;
  int n = 1000;
  Equation* eq = new Equation(n);
  GalerkinRitzMethod(eq);

  x = 0.0, h = 1.0 / n;
  for (int i = 0; i <= n; ++i, x += h)
    printf("%lf ", x);
  printf("\n\n");


  //uh
  printf("\n\ngalerkin:\ntrue error ");
  for (int i = 0; i < n; ++i)
    printf("%.6lf ", eq->TrueErrorGalerkin(LOCAL, i));
  printf("\n indicator ");
  for (int i = 0; i < n; ++i)
    printf("%.6lf ", eq->EstimateErrorGalerkin(i));
  printf("\n\n");

  //uh1
  for (int i = 0; i <= n; ++i)
    eq->ak[i] *= 1.05;

  printf("galerkin with 5%% error:\ntrue error ");
  for (int i = 0; i < n; ++i)
    printf("%.6lf ", eq->TrueErrorGalerkin(LOCAL, i));
  printf("\n indicator ");
  for (int i = 0; i < n; ++i)
    printf("%.6lf ", eq->EstimateErrorGalerkin(i));
  printf("\n\n");
  //uh2
  printf("Non galerkin:\ntrue error ");
  for (int i = 0; i < n; ++i)
    printf("%.6lf ", eq->TrueErrorNonGalerkin(LOCAL, i));
  printf("\n indicator ");
  for (int i = 0; i < n; ++i)
    printf("%.6lf ", eq->EstimateErrorNonGalerkin(i));
  printf("\n\n");
  delete eq;
  */

  return 0;
}