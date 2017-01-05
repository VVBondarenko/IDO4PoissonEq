#include <grid.h>


void Grid_Init(Grid Task,
               double x0, double x1,
               double y0, double y1,
               int n, double tau)
{
    //set sizes and steps
    Task.x0 = x0;
    Task.x1 = x1;
    Task.y0 = y0;
    Task,y1 = y1;
    Task.n = n;
    Task.tau = tau;

    //define arrays
    int i;
    Task.U      = malloc(n * sizeof(double*));
    Task.dUdx   = malloc(n * sizeof(double*));
    Task.dUdy   = malloc(n * sizeof(double*));
    for(i=0; i<n; i++)
    {
        Task.U[i]       = malloc(n * sizeof(double));
        Task.dUdx[i]    = malloc(n * sizeof(double));
        Task.dUdy[i]    = malloc(n * sizeof(double));
    }
}

void Grid_InitDirihlet(Grid Task, double (*f)(double, double))
{
    //for all boundary nodes, set values form f to array
    int i;
    for(i = 0; i<Task.n; i++)
    {
        Task.U[i][0] = (*f)(Task.x0+h*i,Task.y0);
        Task.U[i][n-1]=(*f)(Task.x0+h*i,Task.y1);

        Task.U[0][i] = (*f)(Task.x0,Task.y0+h*i);
        Task.U[n-1][i]=(*f)(Task.x1,Task.y0+h*i);
    }
}

void Grid_Plot(Grid Task)
{
    // for all values in U array, print values of its coords and U
    int i,j;
    FILE *output;

    output = fopen("Plot.dat", "w");
//    for(i=Task.x0; i<=Task.x1; i+=Task.h)
//        for(j=Task.y0; j<=Task.y1; j+=Task.h)
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            fprintf(output,"%15.15f %15.15f %15.15f\n",Task.x0+i*Task.h,
                                                       Task.y0+j*Task.h,
                                                       Task.U[i][j]);
}

void Grid_CrossIteration(Grid Task)
{
    int i,j;
    //define a new array for u on (n+1) iteration
    double **nU;
    nU = malloc(Task.n * sizeof(double*));
    for(i=0; i<n; i++)
        nU[i] = malloc(Task.n * sizeof(double));
    //for all values apply classical schematic
    for(i=1; i<n-1; i++)
    {
        for(j=1; j<n-1; j++)
        {
            nU[i][j] = Task.U[i][j] + Task.tau/Task.h/Task.h*
                    (Task.U[i+1][j]+Task.U[i-1][j]
                    +Task.U[i][j+1]+Task.U[i][j-1]
                     -4*Task.U[i][j]);
        }
    }
    //move data to Task.U
    memcpy(Task.U, nU, sizeof(nU)); //check if it will work
    //clean temporary arrays
    for(i=0; i<n; i++)
        free((void *)nU[i]);
    free((void *)nU);
}

void Grid_IDO_Iteration(Grid Task)
{
    int i,j;
    //define new temporary arrays for U, dUdx, dUdy on (n+1) iteration
    double **nU, **ndUdx, **ndUdy;
    nU    = malloc(Task.n * sizeof(double*));
    ndUdx = malloc(Task.n * sizeof(double*));
    ndUdy = malloc(Task.n * sizeof(double*));
    for(i=0; i<n; i++)
    {
        nU[i]    = malloc(Task.n * sizeof(double));
        ndUdx[i] = malloc(Task.n * sizeof(double));
        ndUdy[i] = malloc(Task.n * sizeof(double));
    }
    //for all values in U apply scheme
    double c, C;
    double d, D;
    for(i=1; i<n-1; i++)
    {
        for(j=1; j<n-1; j++)
        {
            nU[i][j] = Task.U[i][j] + Task.tau/Task.h/Task.h*
                    (Task.U[i+1][j]+Task.U[i-1][j]
                    +Task.U[i][j+1]+Task.U[i][j-1]
                     -4*Task.U[i][j]);
        }
    }
    // --||-- in dUdx

    // --||-- in dUdy

    // update the data in Task

    //clean temporary arrays

}
