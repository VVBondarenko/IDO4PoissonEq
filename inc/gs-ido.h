#ifndef GSIDO_H
#define GSIDO_H
#include <grid.h>

void GSIDO_InitGrid(Grid* Task,
                    double x0, double x1,
                    double y0, double y1,
                    int n);
void GSIDO_Iteration_w_func(Grid* Task, double omega,
                            double (*func)(double, double));
void GSIDO_IterationSet_w_func(Grid* Task, double omega,
                               double (*func)(double, double),
                               int N);

#endif // GSIDO_H
