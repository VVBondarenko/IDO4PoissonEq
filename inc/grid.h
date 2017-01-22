#ifndef GRID_H
#define GRID_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

typedef struct Grid {
    double **U;
    double **dUdx;
    double **dUdy;
    double **dUdxdy;

    int n; //quantity of nodes per side
    double h, tau;
    double x0, x1;
    double y0, y1;
} Grid;
void Grid_Init          (Grid *Task,
                         double x0, double x1,
                         double y0, double y1,
                         int n,
                         double tau);
void Grid_InitByFunction(Grid *Target, double (*func)(double, double),
                         double x0, double x1,
                         double y0, double y1,
                         int n);
void Grid_InitDirihlet  (Grid *Task, double (*f)(double, double));
void Grid_InitDirihlet_w_d(Grid *Task, double (*f)(double, double));
void Grid_InitDirihlet_w_derivatives(Grid *Task, double (*f)(double, double));
double Grid_Interpolate (Grid *Task, double x, double y); //TBD

void Grid_Plot          (Grid *Task);
void Grid_Plot_error    (Grid *Task, double (*exact)(double,double));
void Grid_print_error   (Grid *Task, double (*exact)(double,double));

#endif // GRID_H
