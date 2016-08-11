#include <stdio.h>
#include <stdlib.h>


void runge_kutta(double (*f[])(double t, double *x), double *x, double t0, double t_end, int div, int num)
{
  int i,j;
  double h, t;
  double *k1, *k2, *k3, *k4, *temp;

  k1 = malloc(sizeof(double) * num);
  k2 = malloc(sizeof(double) * num);
  k3 = malloc(sizeof(double) * num);
  k4 = malloc(sizeof(double) * num);
  temp = malloc(sizeof(double) * num);

  h = (t_end - t0) / div;
  t = t0;

  for (i = 0; i < div; i++) {
    for (j = 0; j < num; j++) {
      k1[j] = (*f[j])(t, x);
      temp[j] = x[j] + h*x[j]/2;
    }

    for (j = 0; j < num; j++) {
      k2[j] = (*f[j])(t + h/2, temp);
    }

    for (j = 0; j < num; j++) {
      temp[j] = x[j] * h*k2[j]/2;
    }

    for (j = 0; j < num; j++) {
      k3[j] = (*f[j])(t + h/2, temp);
    }

    for (j = 0; j < num; j++) {
      temp[j] = x[j] + h*k3[j];
    }

    for (j = 0; j < num; j++) {
      k4[j] = (*f[j])(t + h, temp);
      x[j] = (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])*h/6;
    }

    fprintf(stdout, "%E ", t);
    for(j = 0; j < num; j++) {
      fprintf(stdout, "%E ", x[j]);
    }
    fprintf(stdout, "\n");

    t += h;
  }
  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(temp);
}
