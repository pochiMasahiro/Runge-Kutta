#ifndef _RUNGE_KUTTA_H
#define _RUNGE_KUTTA_H

void runge_kutta(double (*f[])(double t, double *x), double *x, double t0, double t_end, int div, int num);

#endif
