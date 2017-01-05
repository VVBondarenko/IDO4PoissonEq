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
void Grid_InitDirihlet  (Grid *Task, double (*f)(double, double));
void Grid_Plot          (Grid *Task);
void Grid_CrossIteration(Grid *Task);
void Grid_IDO_Iteration (Grid *Task);
#endif // GRID_H
