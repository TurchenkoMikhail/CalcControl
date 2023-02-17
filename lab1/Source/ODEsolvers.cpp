#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ABS(x) ((x)>0?(x):-(x))
#define MAX(x,y) ((x)>(y)?(x):(y)) 

#define P(x) (-9.0)
#define Q(x) (-10.0)

double Adams1(double a, double b, unsigned long long n, double y0, double z0) {
  double yi, yiNext;
  unsigned long long i;
  double h = (b - a) / n;
  double zi, ziNextRough, ziNext;
  double yiNextRough;
  double xi = a;

  yi = y0;
  zi = z0;

  for (i = 0; i < n; i++, xi += h) {
    yiNextRough = yi + h * zi;
    ziNextRough = zi + h * (-P(xi) * zi - Q(xi) * yi);

    yiNext = yi + h * ziNextRough;
    ziNext = zi + h * (-P(xi + h) * ziNextRough - Q(xi + h) * yiNextRough);

    zi = ziNext;
    yi = yiNext;
  }

  return yiNext;
}

double FixedEulerMethod(double a, double b, unsigned long long n, double y0, double z0) {
  unsigned long long i;
  double x0 = a, xi = x0;
  double ff, gg, h;
  double yi, yinext, zi, zinext;

  h = (b - a) / n;
  xi = a;

  yi = y0; //начальное условие
  zi = z0;

  //заполнение массива значений yk
  for (i = 0; i < n; i++, xi += h) {

    ff = -P(xi) * zi - Q(xi) * yi;  // = func(xi, yi, zi)
    gg = zi;

    zinext = zi + h / 2.0 * (ff + (-P(xi+h) * (zi + h*ff) - Q(xi+h) * (yi + h*ff)));
    yinext = yi + h / 2.0 * (gg + (zi + h*gg));

    yi = yinext;
    zi = zinext;
  }

  return yinext;
}