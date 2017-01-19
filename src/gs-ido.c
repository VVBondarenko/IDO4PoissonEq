#include <gs-ido.h>

void GSIDO_InitGrid(Grid *Task,
                    double x0, double x1,
                    double y0, double y1,
                    int n)
{
    Task->x0 = x0;
    Task->x1 = x1;
    Task->y0 = y0;
    Task->y1 = y1;
    Task->n = n;

    Task->U      = malloc(n * sizeof(double*));
    Task->dUdx   = malloc((n-1) * sizeof(double*));
    Task->dUdy   = malloc(n * sizeof(double*));
    Task->dUdxdy = malloc((n-1) * sizeof(double*));

    int i;
    for(i=0; i<n; i++)
    {
        Task->U[i]       = malloc(n * sizeof(double));
        Task->dUdx[i]    = malloc(n * sizeof(double));
        Task->dUdy[i]    = malloc((n-1) * sizeof(double));
        Task->dUdxdy[i]  = malloc((n-1) * sizeof(double));
    }
}

void GSIDO_Iteration_w_func(Grid *Task, double omega, double (*func)(double, double))
{
    int i,j, n = Task->n;
    double **nU, **ndUdx, **ndUdy, **ndUdxdy;

    nU      = malloc(n * sizeof(double*));
    ndUdx   = malloc((n-1) * sizeof(double*));
    ndUdy   = malloc(n * sizeof(double*));
    ndUdxdy = malloc((n-1) * sizeof(double*));

    for(i=0; i<n; i++)
    {
        nU[i]       = malloc(n * sizeof(double));
        ndUdy[i]    = malloc((n-1) * sizeof(double));
        if(n<n-1)
        {
            ndUdx[i]    = malloc(n * sizeof(double));
            ndUdxdy[i]  = malloc((n-1) * sizeof(double));
        }
    }




}
