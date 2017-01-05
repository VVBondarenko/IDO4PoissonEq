#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Grid {
    double **U;
    double **dUdx;
    double **dUdy;

    int n; //quantity of nodes per side
    double h, tau;
    double x0, x1;
    double y0, y1;
} Grid;
void Grid_Init          (Grid Task,
                         double x0, double x1,
                         double y0, double y1,
                         int n,
                         double tau);
void Grid_InitDirihlet  (Grid Task, double (*f)(double, double));
void Grid_Plot          (Grid Task);
void Grid_Iteration     (Grid Task);

//void iteration(Grid Task);
int main()
{

    return 0;
}
