#ifndef GRID_H
#define GRID_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


typedef struct Grid {
    double **U;
    double **dUdx;
    double **dUdy;

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
void Grid_Plot          (Grid *Task);
double Grid_Interpolate (Grid *Task, double x, double y); //TBD

void Grid_CrossIteration(Grid *Task);

void Grid_IDO_InitDeriv (Grid *Task);
void Grid_IDO_Iteration (Grid *Task);
void Grid_IDO_Iteration_with_f(Grid *Task, Grid *func);
#endif // GRID_H
