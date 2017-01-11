#include <grid.h>


void Grid_Init(Grid *Task,
               double x0, double x1,
               double y0, double y1,
               int n, double tau)
{
    //set sizes and steps
    Task->x0 = x0;
    Task->x1 = x1;
    Task->y0 = y0;
    Task->y1 = y1;
    Task->n = n;
    Task->tau = tau;
    Task->h = (x1-x0)/(n-1);
    //define arrays
    int i,j;
    Task->U      = malloc(n * sizeof(double*));
    Task->dUdx   = malloc(n * sizeof(double*));
    Task->dUdy   = malloc(n * sizeof(double*));
    Task->dUdxdy   = malloc(n * sizeof(double*));
    for(i=0; i<n; i++)
    {
        Task->U[i]       = malloc(n * sizeof(double));
        Task->dUdx[i]    = malloc(n * sizeof(double));
        Task->dUdy[i]    = malloc(n * sizeof(double));
        Task->dUdxdy[i]    = malloc(n * sizeof(double));
    }

    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            Task->U[i][j]=0.;
            Task->dUdx[i][j]=0.;
            Task->dUdy[i][j]=0.;
            Task->dUdxdy[i][j]=0.;
        }
    }
}

void Grid_InitDirihlet(Grid *Task, double (*f)(double, double))
{
    //for all boundary nodes, set values form f to array
    int i;
    for(i = 0; i<Task->n; i++)
    {
        Task->U[i][0] = (*f)(Task->x0+Task->h*i,Task->y0);
        Task->U[i][Task->n-1]=(*f)(Task->x0+Task->h*i,Task->y1);

        Task->U[0][i] = (*f)(Task->x0,Task->y0+Task->h*i);
        Task->U[Task->n-1][i]=(*f)(Task->x1,Task->y0+Task->h*i);
    }
}

void Grid_Plot(Grid *Task)
{
    // for all values in U array, print values of its coords and U
    int i,j;
    FILE *output;

    output = fopen("Plot.dat", "w");
//    for(i=Task->x0; i<=Task->x1; i+=Task->h)
//        for(j=Task->y0; j<=Task->y1; j+=Task->h)
    for(i = 0; i < Task->n; i++)
        for(j = 0; j < Task->n; j++)
            fprintf(output,"%15.15f %15.15f %15.15f\n",Task->x0+i*Task->h,
                                                       Task->y0+j*Task->h,
                                                       Task->U[i][j]);
}

void Grid_CrossIteration(Grid *Task, double omega)
{
    int i,j;
    //define a new array for u on (n+1) iteration
    double **nU;
    nU = malloc(Task->n * sizeof(double*));
    for(i=0; i<Task->n; i++)
        nU[i] = malloc(Task->n * sizeof(double));
    //for all values apply classical schematic
    for(i=1; i<Task->n-1; i++)
    {
        for(j=1; j<Task->n-1; j++)
        {
            nU[i][j] = (1-omega)*Task->U[i][j]
                    +omega*0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                                +Task->U[i][j+1]+Task->U[i][j-1]);
        }
    }
    //move data to Task->U
    for(i=1; i<Task->n-1; i++)
    {
        for(j=1; j<Task->n-1; j++)
        {
            Task->U[i][j] = nU[i][j];
        }
    }

    //clean temporary arrays
    for(i=0; i<Task->n; i++)
        free((void *)nU[i]);
    free((void *)nU);
}

void Grid_CrossIteration_w_force(Grid *Task, double omega, Grid *force)
{
    int i,j;
    //define a new array for u on (n+1) iteration
    double **nU;
    nU = malloc(Task->n * sizeof(double*));
    for(i=0; i<Task->n; i++)
        nU[i] = malloc(Task->n * sizeof(double));
    //for all values apply classical schematic
    #pragma omp parallel for shared(nU,omega,Task,force) private(i,j)
    for(i=1; i<Task->n-1; i++)
    {
        for(j=1; j<Task->n-1; j++)
        {
            nU[i][j] = (1-omega)*Task->U[i][j]
                    +omega*0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                                +Task->U[i][j+1]+Task->U[i][j-1]
                                -force->U[i][j]*Task->h*Task->h);
        }
    }
    //move data to Task->U
    for(i=1; i<Task->n-1; i++)
    {
        for(j=1; j<Task->n-1; j++)
        {
            Task->U[i][j] = nU[i][j];
        }
    }

    //clean temporary arrays
    for(i=0; i<Task->n; i++)
        free((void *)nU[i]);
    free((void *)nU);
}

void Grid_IDO_Iteration(Grid *Task)
{
    int i,j;
    //define new temporary arrays for U, dUdx, dUdy on (n+1) iteration
    double **nU, **ndUdx, **ndUdy;
    nU    = malloc(Task->n * sizeof(double*));
    ndUdx = malloc(Task->n * sizeof(double*));
    ndUdy = malloc(Task->n * sizeof(double*));
    for(i=0; i<Task->n; i++)
    {
        nU[i]    = malloc(Task->n * sizeof(double));
        ndUdx[i] = malloc(Task->n * sizeof(double));
        ndUdy[i] = malloc(Task->n * sizeof(double));
    }
    //for all values in U apply scheme
    double c, C;
    double d, D;


    for(i=1; i<Task->n-1; i++)
    {
        for(j=1; j<Task->n-1; j++)
        {
            c=(1.25/Task->h*(Task->U[i+1][j]-Task->U[i-1][j])
                      -0.25*(Task->dUdx[i+1][j]
                          +8*Task->dUdx[i][j]
                            +Task->dUdx[i-1][j]))/Task->h/Task->h;
            C=(1.25/Task->h*(Task->U[i][j+1]-Task->U[i][j-1])
                      -0.25*(Task->dUdy[i][j+1]
                          +8*Task->dUdy[i][j]
                            +Task->dUdy[i][j-1]))/Task->h/Task->h;

            d=(1./Task->h*(Task->U[i+1][j]
                       -2.*Task->U[i][j]
                          +Task->U[i-1][j])
                    -0.25*(Task->dUdx[i+1][j]-Task->dUdx[i-1][j]))/Task->h;

            D=(1./Task->h*(Task->U[i][j+1]
                       -2.*Task->U[i][j]
                          +Task->U[i][j-1])
                    -0.25*(Task->dUdy[i][j+1]-Task->dUdy[i][j-1]))/Task->h;


            nU[i][j] = Task->U[i][j]
                    + 2.*Task->tau*(d+D);

            ndUdx[i][j] = Task->dUdx[i][j]
                    +Task->tau*(6.*c+(Task->dUdx[i][j+1]
                                   -2*Task->dUdx[i][j]+
                                      Task->dUdx[i][j-1])/Task->h/Task->h);
            ndUdy[i][j] = Task->dUdy[i][j]
                    +Task->tau*(6.*C+(Task->dUdy[i+1][j]
                                   -2*Task->dUdy[i][j]+
                                      Task->dUdy[i-1][j])/Task->h/Task->h);
        }
    }

    // update the data in Task
    for(i=1; i<Task->n-1; i++)
    {
        for(j=1; j<Task->n-1; j++)
        {
            Task->U[i][j] = nU[i][j];
            Task->dUdx[i][j] = ndUdx[i][j];
            Task->dUdy[i][j] = ndUdy[i][j];
        }
    }

    //clean temporary arrays
    for(i=0; i<Task->n; i++)
    {
        free((void *)nU[i]);
        free((void *)ndUdx[i]);
        free((void *)ndUdy[i]);
    }
    free((void *)nU);
    free((void *)ndUdx);
    free((void *)ndUdy);
}

void Grid_IDO_IterationOriginal(Grid *Task, double omega)
{
    int i,j;
    //define new temporary arrays for U, dUdx, dUdy on (n+1) iteration
    double **nU, **ndUdx, **ndUdy, **ndUdxdy;
    nU    = malloc(Task->n * sizeof(double*));
    ndUdx = malloc(Task->n * sizeof(double*));
    ndUdy = malloc(Task->n * sizeof(double*));
    ndUdxdy = malloc(Task->n * sizeof(double*));
    for(i=0; i<Task->n; i++)
    {
        nU[i]    = malloc(Task->n * sizeof(double));
        ndUdx[i] = malloc(Task->n * sizeof(double));
        ndUdy[i] = malloc(Task->n * sizeof(double));
        ndUdxdy[i] = malloc(Task->n * sizeof(double));
    }
    //for all values in U apply scheme
    double h = Task->h;
//    double omega = 1.;
    for(i=1; i<Task->n-1; i++)
    {
        for(j=1; j<Task->n-1; j++)
        {
/*            nU[i][j] = (2/h/h*(Task->U[i+1][j]+Task->U[i-1][j]
                              +Task->U[i][j+1]+Task->U[i][j-1])
                       -0.5/h*(Task->dUdx[i+1][j]-Task->dUdx[i-1][j]
                              +Task->dUdy[i][j+1]-Task->dUdy[i][j-1]))/8.*h*h;

            ndUdx[i][j] = (7.5*(Task->U[i+1][j]-Task->U[i-1][j])/h/h/h
                          -1.5*(Task->dUdx[i+1][j]+Task->dUdx[i-1][j])/h/h
                           +2.*(Task->dUdx[i][j+1]+Task->dUdx[i][j-1])/h/h
                          -0.5*(Task->dUdxdy[i][j+1]-Task->dUdxdy[i][j-1])/h)
                    /(12/h/h+4/h/h);
            ndUdy[i][j] = (7.5*(Task->U[i][j+1]-Task->U[i][j-1])/h/h/h
                          -1.5*(Task->dUdy[i][j+1]+Task->dUdy[i][j-1])/h/h
                           +2.*(Task->dUdy[i+1][j]+Task->dUdy[i-1][j])/h/h
                          -0.5*(Task->dUdxdy[i+1][j]-Task->dUdxdy[i-1][j])/h)
                    /(12/h/h+4/h/h);
            ndUdxdy[i][j] = (7.5*(Task->dUdx[i][j+1]-Task->dUdx[i][j-1]
                                 +Task->dUdy[i+1][j]-Task->dUdy[i-1][j])/h/h/h
                            -1.5*(Task->dUdxdy[i][j+1]+Task->dUdxdy[i][j-1]
                                 +Task->dUdxdy[i+1][j]+Task->dUdxdy[i-1][j])/h/h
                            )/24*h*h;*/
            nU[i][j] = (1-omega)*Task->U[i][j]
                +omega*(2/h/h*(Task->U[i+1][j]+Task->U[i-1][j]
                              +Task->U[i][j+1]+Task->U[i][j-1])
                       -0.5/h*(Task->dUdx[i+1][j]-Task->dUdx[i-1][j]
                              +Task->dUdy[i][j+1]-Task->dUdy[i][j-1]))/8.*h*h;

            ndUdx[i][j] = (1-omega)*Task->dUdx[i][j]
                   +omega*(7.5*(Task->U[i+1][j]-Task->U[i-1][j])/h/h/h
                          -1.5*(Task->dUdx[i+1][j]+Task->dUdx[i-1][j])/h/h
                           +2.*(Task->dUdx[i][j+1]+Task->dUdx[i][j-1])/h/h
                          -0.5*(Task->dUdxdy[i][j+1]-Task->dUdxdy[i][j-1])/h)
                    /(12/h/h+4/h/h);
            ndUdy[i][j] = (1-omega)*Task->dUdy[i][j]
                   +omega*(7.5*(Task->U[i][j+1]-Task->U[i][j-1])/h/h/h
                          -1.5*(Task->dUdy[i][j+1]+Task->dUdy[i][j-1])/h/h
                           +2.*(Task->dUdy[i+1][j]+Task->dUdy[i-1][j])/h/h
                          -0.5*(Task->dUdxdy[i+1][j]-Task->dUdxdy[i-1][j])/h)
                    /(12/h/h+4/h/h);
            ndUdxdy[i][j] = (1-omega)*Task->dUdxdy[i][j]
                     +omega*(7.5*(Task->dUdx[i][j+1]-Task->dUdx[i][j-1]
                                 +Task->dUdy[i+1][j]-Task->dUdy[i-1][j])/h/h/h
                            -1.5*(Task->dUdxdy[i][j+1]+Task->dUdxdy[i][j-1]
                                 +Task->dUdxdy[i+1][j]+Task->dUdxdy[i-1][j])/h/h
                            )/24*h*h;
        }
    }

    // update the data in Task
    for(i=1; i<Task->n-1; i++)
    {
        for(j=1; j<Task->n-1; j++)
        {
            Task->U[i][j] = nU[i][j];
            Task->dUdx[i][j] = ndUdx[i][j];
            Task->dUdy[i][j] = ndUdy[i][j];
            Task->dUdxdy[i][j] = ndUdxdy[i][j];
        }
    }

    //clean temporary arrays
    for(i=0; i<Task->n; i++)
    {
        free((void *)nU[i]);
        free((void *)ndUdx[i]);
        free((void *)ndUdy[i]);
        free((void *)ndUdxdy[i]);
    }
    free((void *)nU);
    free((void *)ndUdx);
    free((void *)ndUdy);
    free((void *)ndUdxdy);
}

void Grid_IDO_IterationModified(Grid *Task, double omega)
{
    int i,j,p;
    //define new temporary arrays for U, dUdx, dUdy on (n+1) iteration
    double **nU, **ndUdx, **ndUdy, **ndUdxdy;
    nU    = malloc(Task->n * sizeof(double*));
    ndUdx = malloc(Task->n * sizeof(double*));
    ndUdy = malloc(Task->n * sizeof(double*));
    ndUdxdy = malloc(Task->n * sizeof(double*));
    for(i=0; i<Task->n; i++)
    {
        nU[i]    = malloc(Task->n * sizeof(double));
        ndUdx[i] = malloc(Task->n * sizeof(double));
        ndUdy[i] = malloc(Task->n * sizeof(double));
        ndUdxdy[i] = malloc(Task->n * sizeof(double));
    }
    //for all values in U apply scheme
    double h = Task->h;


    for(p=1; p<Task->n-1; p++)
    {
        i=p;j=1;
        nU[i][j] = 0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                   +Task->U[i][j+1]+Task->U[i][j-1]);
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j])/h;
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1])/h;
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                             -Task->U[i-1][j+1]
                             -Task->U[i+1][j-1]
                             +Task->U[i-1][j-1])/h/h;
        i=p;j=Task->n-2;
        nU[i][j] = 0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                   +Task->U[i][j+1]+Task->U[i][j-1]);
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j])/h;
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1])/h;
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                             -Task->U[i-1][j+1]
                             -Task->U[i+1][j-1]
                             +Task->U[i-1][j-1])/h/h;
        i=1;j=p;
        nU[i][j] = 0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                   +Task->U[i][j+1]+Task->U[i][j-1]);
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j])/h;
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1])/h;
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                             -Task->U[i-1][j+1]
                             -Task->U[i+1][j-1]
                             +Task->U[i-1][j-1])/h/h;
        i=Task->n-2;j=p;
        nU[i][j] = 0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                   +Task->U[i][j+1]+Task->U[i][j-1]);
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j])/h;
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1])/h;
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                             -Task->U[i-1][j+1]
                             -Task->U[i+1][j-1]
                             +Task->U[i-1][j-1])/h/h;

    }

    for(i=2; i<Task->n-2; i++)
    {
        for(j=2; j<Task->n-2; j++)
        {

                nU[i][j] = (1-omega)*Task->U[i][j]
                        +omega*(2/h/h*(Task->U[i+1][j]+Task->U[i-1][j]
                        +Task->U[i][j+1]+Task->U[i][j-1])
                        -0.5/h*(Task->dUdx[i+1][j]-Task->dUdx[i-1][j]
                        +Task->dUdy[i][j+1]-Task->dUdy[i][j-1]))/8.*h*h;

                ndUdx[i][j] = (1-omega)*Task->dUdx[i][j]
                        +omega*(7.5*(Task->U[i+1][j]-Task->U[i-1][j])/h/h/h
                        -1.5*(Task->dUdx[i+1][j]+Task->dUdx[i-1][j])/h/h
                        +2.*(Task->dUdx[i][j+1]+Task->dUdx[i][j-1])/h/h
                        -0.5*(Task->dUdxdy[i][j+1]-Task->dUdxdy[i][j-1])/h)
                        /(12/h/h+4/h/h);
                ndUdy[i][j] = (1-omega)*Task->dUdy[i][j]
                        +omega*(7.5*(Task->U[i][j+1]-Task->U[i][j-1])/h/h/h
                        -1.5*(Task->dUdy[i][j+1]+Task->dUdy[i][j-1])/h/h
                        +2.*(Task->dUdy[i+1][j]+Task->dUdy[i-1][j])/h/h
                        -0.5*(Task->dUdxdy[i+1][j]-Task->dUdxdy[i-1][j])/h)
                        /(12/h/h+4/h/h);
                ndUdxdy[i][j] = (1-omega)*Task->dUdxdy[i][j]
                        +omega*(7.5*(Task->dUdx[i][j+1]-Task->dUdx[i][j-1]
                        +Task->dUdy[i+1][j]-Task->dUdy[i-1][j])/h/h/h
                        -1.5*(Task->dUdxdy[i][j+1]+Task->dUdxdy[i][j-1]
                        +Task->dUdxdy[i+1][j]+Task->dUdxdy[i-1][j])/h/h
                        )/24*h*h;
        }
    }

    // update the data in Task
    for(i=1; i<Task->n-1; i++)
    {
        for(j=1; j<Task->n-1; j++)
        {
            Task->U[i][j] = nU[i][j];
            Task->dUdx[i][j] = ndUdx[i][j];
            Task->dUdy[i][j] = ndUdy[i][j];
            Task->dUdxdy[i][j] = ndUdxdy[i][j];
        }
    }

    //clean temporary arrays
    for(i=0; i<Task->n; i++)
    {
        free((void *)nU[i]);
        free((void *)ndUdx[i]);
        free((void *)ndUdy[i]);
        free((void *)ndUdxdy[i]);
    }
    free((void *)nU);
    free((void *)ndUdx);
    free((void *)ndUdy);
    free((void *)ndUdxdy);
}

void Grid_IDO_IterationModified_w_f(Grid *Task, double omega, Grid *force)
{
    int i,j,p;
    //define new temporary arrays for U, dUdx, dUdy on (n+1) iteration
    double **nU, **ndUdx, **ndUdy, **ndUdxdy;
    nU    = malloc(Task->n * sizeof(double*));
    ndUdx = malloc(Task->n * sizeof(double*));
    ndUdy = malloc(Task->n * sizeof(double*));
    ndUdxdy = malloc(Task->n * sizeof(double*));
    for(i=0; i<Task->n; i++)
    {
        nU[i]    = malloc(Task->n * sizeof(double));
        ndUdx[i] = malloc(Task->n * sizeof(double));
        ndUdy[i] = malloc(Task->n * sizeof(double));
        ndUdxdy[i] = malloc(Task->n * sizeof(double));
    }
    //for all values in U apply scheme
    double h = Task->h;


    for(p=1; p<Task->n-1; p++)
    {
        i=p;j=1;
        nU[i][j] = 0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                   +Task->U[i][j+1]+Task->U[i][j-1]-force->U[i][j]*Task->h*Task->h);
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j]-force->dUdx[i][j]*Task->h*Task->h)/h;
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1]-force->dUdy[i][j]*Task->h*Task->h)/h;
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                             -Task->U[i-1][j+1]
                             -Task->U[i+1][j-1]
                             +Task->U[i-1][j-1]-force->dUdxdy[i][j]*Task->h*Task->h)/h/h;
        i=p;j=Task->n-2;
        nU[i][j] = 0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                   +Task->U[i][j+1]+Task->U[i][j-1]-force->U[i][j]*Task->h*Task->h);
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j]-force->dUdx[i][j]*Task->h*Task->h)/h;
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1]-force->dUdy[i][j]*Task->h*Task->h)/h;
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                             -Task->U[i-1][j+1]
                             -Task->U[i+1][j-1]
                             +Task->U[i-1][j-1]-force->dUdxdy[i][j]*Task->h*Task->h)/h/h;
        i=1;j=p;
        nU[i][j] = 0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                   +Task->U[i][j+1]+Task->U[i][j-1]-force->U[i][j]*Task->h*Task->h);
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j]-force->dUdx[i][j]*Task->h*Task->h)/h;
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1]-force->dUdy[i][j]*Task->h*Task->h)/h;
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                             -Task->U[i-1][j+1]
                             -Task->U[i+1][j-1]
                             +Task->U[i-1][j-1]-force->dUdxdy[i][j]*Task->h*Task->h)/h/h;
        i=Task->n-2;j=p;
        nU[i][j] = 0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                   +Task->U[i][j+1]+Task->U[i][j-1]-force->U[i][j]*Task->h*Task->h);
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j]-force->dUdx[i][j]*Task->h*Task->h)/h;
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1]-force->dUdy[i][j]*Task->h*Task->h)/h;
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                             -Task->U[i-1][j+1]
                             -Task->U[i+1][j-1]
                             +Task->U[i-1][j-1]-force->dUdxdy[i][j]*Task->h*Task->h)/h/h;

    }

    for(i=2; i<Task->n-2; i++)
    {
        for(j=2; j<Task->n-2; j++)
        {

                nU[i][j] = (1-omega)*Task->U[i][j]
                        +omega*(2/h/h*(Task->U[i+1][j]+Task->U[i-1][j]
                        +Task->U[i][j+1]+Task->U[i][j-1])
                        -0.5/h*(Task->dUdx[i+1][j]-Task->dUdx[i-1][j]
                        +Task->dUdy[i][j+1]-Task->dUdy[i][j-1])-force->U[i][j]*Task->h*Task->h)/8.*h*h;

                ndUdx[i][j] = (1-omega)*Task->dUdx[i][j]
                        +omega*(7.5*(Task->U[i+1][j]-Task->U[i-1][j])/h/h/h
                        -1.5*(Task->dUdx[i+1][j]+Task->dUdx[i-1][j])/h/h
                        +2.*(Task->dUdx[i][j+1]+Task->dUdx[i][j-1])/h/h
                        -0.5*(Task->dUdxdy[i][j+1]-Task->dUdxdy[i][j-1])/h-force->dUdx[i][j]*Task->h*Task->h)
                        /(12/h/h+4/h/h);
                ndUdy[i][j] = (1-omega)*Task->dUdy[i][j]
                        +omega*(7.5*(Task->U[i][j+1]-Task->U[i][j-1])/h/h/h
                        -1.5*(Task->dUdy[i][j+1]+Task->dUdy[i][j-1])/h/h
                        +2.*(Task->dUdy[i+1][j]+Task->dUdy[i-1][j])/h/h
                        -0.5*(Task->dUdxdy[i+1][j]-Task->dUdxdy[i-1][j])/h-force->dUdy[i][j]*Task->h*Task->h)
                        /(12/h/h+4/h/h);
                ndUdxdy[i][j] = (1-omega)*Task->dUdxdy[i][j]
                        +omega*(7.5*(Task->dUdx[i][j+1]-Task->dUdx[i][j-1]
                        +Task->dUdy[i+1][j]-Task->dUdy[i-1][j])/h/h/h
                        -1.5*(Task->dUdxdy[i][j+1]+Task->dUdxdy[i][j-1]
                        +Task->dUdxdy[i+1][j]+Task->dUdxdy[i-1][j])/h/h
                        -force->dUdxdy[i][j]*Task->h*Task->h)/24*h*h;
        }
    }

    // update the data in Task
    for(i=1; i<Task->n-1; i++)
    {
        for(j=1; j<Task->n-1; j++)
        {
            Task->U[i][j] = nU[i][j];
            Task->dUdx[i][j] = ndUdx[i][j];
            Task->dUdy[i][j] = ndUdy[i][j];
            Task->dUdxdy[i][j] = ndUdxdy[i][j];
        }
    }

    //clean temporary arrays
    for(i=0; i<Task->n; i++)
    {
        free((void *)nU[i]);
        free((void *)ndUdx[i]);
        free((void *)ndUdy[i]);
        free((void *)ndUdxdy[i]);
    }
    free((void *)nU);
    free((void *)ndUdx);
    free((void *)ndUdy);
    free((void *)ndUdxdy);
}

void Grid_IDO_InitDeriv(Grid *Task)
{
    int i,j;
    for(i=0;i<Task->n;i++)
    {
        for(j=0;j<Task->n;j++)
        {
            //init dudx and dudy
            if(i!=0 && i!=Task->n-1)
            {
//                Task->dUdx[i][j] =
//                        0.5*(Task->U[i+1][j]-Task->U[i-1][j])/Task->h;
            }
            else
            {
                if(i==0)
                    Task->dUdx[0][j]=(Task->U[1][j]-Task->U[0][j])/Task->h;
                if(i==Task->n-1)
                    Task->dUdx[Task->n-1][j]=
                    (Task->U[Task->n-1][j]-Task->U[Task->n-2][j])/Task->h;
            }

            if(j!=0 && j!=Task->n-1)
            {
//                Task->dUdy[i][j] =
//                        0.5*(Task->U[i][j+1]-Task->U[i][j-1])/Task->h;
                if(i==0 || i==Task->n-1)
                    Task->dUdxdy[i][j]=0.5*(Task->dUdx[i][j+1]
                                           -Task->dUdx[i][j-1])/Task->h;

            }
            else
            {
                if(j==0)
                {
                    Task->dUdy[i][0]=(Task->U[i][1]-Task->U[i][0])/Task->h;
                    Task->dUdxdy[i][0]=(Task->dUdx[i][1]
                                     -Task->dUdx[i][0])/Task->h;
                }
                if(j==Task->n-1)
                {
                    Task->dUdy[i][Task->n-1]=
                            (Task->U[i][Task->n-1]-Task->U[i][Task->n-2])/Task->h;
                    Task->dUdxdy[i][Task->n-1]=(Task->dUdx[i][Task->n-1]
                            -Task->dUdx[i][Task->n-2])/Task->h;
                }
            }

        }
    }
}


void Grid_InitByFunction(Grid *Target, double (*func)(double, double),
                         double x0, double x1,
                         double y0, double y1,
                         int n)
{
    Grid_Init(Target,x0,x1,y0,y1,n,0.);

    int i,j;
    double h = Target->h;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            Target->U[i][j] = (*func)(x0+i*h,y0+j*h);
            if(i!=0 && i!=Target->n-1)
            {
                Target->dUdx[i][j] = ((*func)(x0+i*h+h,y0+j*h)
                                     -(*func)(x0+i*h-h,y0+j*h))/h/h*0.5;
            }
            else
            {
                if(i==0)
                    Target->dUdx[i][j] = ((*func)(x0+h,y0+j*h)
                                         -(*func)(x0,y0+j*h))/h/h;
                if(i==Target->n-1)
                    Target->dUdx[i][j] = ((*func)(x1,y0+j*h)
                                         -(*func)(x1-h,y0+j*h))/h/h;
            }


            if(j!=0 && j!=Target->n-1)
            {
                Target->dUdy[i][j] = ((*func)(x0+i*h,y0+j*h+h)
                                     -(*func)(x0+i*h,y0+j*h-h))/h/h*0.5;
            }
            else
            {
                if(j==0)
                    Target->dUdy[i][j] = ((*func)(x0+i*h,y0+h)
                                         -(*func)(x0+i*h,y0))/h/h;
                if(j==Target->n-1)
                    Target->dUdy[i][j] = ((*func)(x0+i*h,y1)
                                         -(*func)(x0+i*h,y1-h))/h/h;
            }
        }
    }
}

void Grid_Plot_error(Grid *Task, double (*exact)(double, double))
{
    // for all values in U array, print values of its coords and U
    int i,j;
    double maxErr = 0.;
    FILE *output;

    output = fopen("errorPlot.dat", "w");
//    for(i=Task->x0; i<=Task->x1; i+=Task->h)
//        for(j=Task->y0; j<=Task->y1; j+=Task->h)
    for(i = 0; i < Task->n; i++)
    {
        for(j = 0; j < Task->n; j++)
        {
            fprintf(output,"%15.15f %15.15f %15.15f\n",Task->x0+i*Task->h,
                    Task->y0+j*Task->h,
                    fabs(Task->U[i][j]
                         -(*exact)(Task->x0+i*Task->h,
                                   Task->y0+j*Task->h)));
            maxErr = fmax(maxErr,fabs(Task->U[i][j]
                                      -(*exact)(Task->x0+i*Task->h,
                                                Task->y0+j*Task->h)));
        }
    }
    printf("%10.10f\n",maxErr);
    fclose(output);
}

void Grid_print_error(Grid *Task, double (*exact)(double, double))
{
    // for all values in U array, print values of its coords and U
    int i,j;
    double maxErr = 0.;
//    FILE *output;

//    output = fopen("errorPlot.dat", "w");
//    for(i=Task->x0; i<=Task->x1; i+=Task->h)
//        for(j=Task->y0; j<=Task->y1; j+=Task->h)
    for(i = 0; i < Task->n; i++)
    {
        for(j = 0; j < Task->n; j++)
        {
            maxErr = fmax(maxErr,fabs(Task->U[i][j]
                                      -(*exact)(Task->x0+i*Task->h,
                                                Task->y0+j*Task->h)));
        }
    }
    printf("%10.10e\n",maxErr);
}

void Grid_InitDirihlet_w_derivatives(Grid *Task, double (*f)(double, double))
{
    int i;
    double h = Task->h;
    double mh = Task->h/1024.;
    for(i = 0; i<Task->n; i++)
    {
        Task->U[i][0] = (*f)(Task->x0+h*i,Task->y0);
        Task->U[i][Task->n-1]=(*f)(Task->x0+h*i,Task->y1);

        Task->U[0][i] = (*f)(Task->x0,Task->y0+h*i);
        Task->U[Task->n-1][i]=(*f)(Task->x1,Task->y0+h*i);

        Task->dUdx[i][0] = 0.5*((*f)(Task->x0+h*i+mh,Task->y0)-
                                (*f)(Task->x0+h*i-mh,Task->y0))/mh;
        Task->dUdx[i][Task->n-1]=0.5*((*f)(Task->x0+h*i+mh,Task->y1)
                                     -(*f)(Task->x0+h*i-mh,Task->y1))/mh;

        Task->dUdx[0][i] = ((*f)(Task->x0+mh,Task->y0+h*i)
                           -(*f)(Task->x0-mh,Task->y0+h*i))/mh*0.5;
        Task->dUdx[Task->n-1][i] = ((*f)(Task->x1+mh,Task->y0+h*i)
                                   -(*f)(Task->x1-mh,Task->y0+h*i))/mh*0.5;

        Task->dUdy[i][0] = 0.5*((*f)(Task->x0+h*i,Task->y0+mh)-
                                (*f)(Task->x0+h*i,Task->y0-mh))/mh;
        Task->dUdy[i][Task->n-1]=0.5*((*f)(Task->x0+h*i,Task->y1+mh)
                                     -(*f)(Task->x0+h*i,Task->y1-mh))/mh;

        Task->dUdy[0][i] = ((*f)(Task->x0,Task->y0+h*i+mh)
                           -(*f)(Task->x0,Task->y0+h*i-mh))/mh*0.5;
        Task->dUdy[Task->n-1][i] = ((*f)(Task->x1,Task->y0+h*i+mh)
                                   -(*f)(Task->x1,Task->y0+h*i-mh))/mh*0.5;


        Task->dUdxdy[i][0] =0.25*((*f)(Task->x0+h*i+mh,Task->y0+mh)-
                                  (*f)(Task->x0+h*i+mh,Task->y0-mh)
                                 -(*f)(Task->x0+h*i-mh,Task->y0+mh)+
                                  (*f)(Task->x0+h*i-mh,Task->y0-mh)
                                  )/mh/mh;
        Task->dUdxdy[i][Task->n-1]=0.25*((*f)(Task->x0+h*i+mh,Task->y1+mh)
                                       -(*f)(Task->x0+h*i+mh,Task->y1-mh)
                                       -(*f)(Task->x0+h*i-mh,Task->y1+mh)
                                       +(*f)(Task->x0+h*i-mh,Task->y1-mh)
                                        )/mh/mh;

        Task->dUdxdy[0][i] = ( (*f)(Task->x0+mh,Task->y0+h*i+mh)
                              -(*f)(Task->x0+mh,Task->y0+h*i-mh)
                              -(*f)(Task->x0-mh,Task->y0+h*i+mh)
                              +(*f)(Task->x0-mh,Task->y0+h*i-mh)
                              )/mh/mh*0.25;
        Task->dUdxdy[Task->n-1][i] = ((*f)(Task->x1+mh,Task->y0+h*i+mh)
                                     -(*f)(Task->x1+mh,Task->y0+h*i-mh)
                                     -(*f)(Task->x1-mh,Task->y0+h*i+mh)
                                     +(*f)(Task->x1-mh,Task->y0+h*i-mh)
                                      )/mh/mh*0.25;


    }
}

void Grid_IDO_IterationSet      (Grid *Task, double omega, int N)
{
    int i;
    for(i=0;i<N;i++)
    {
        Grid_IDO_IterationModified(Task,omega);
    }
}

void Grid_IDO_IterationSet_w_f  (Grid *Task, double omega, Grid *force, int N)
{
    int i;
    for(i=0;i<N;i++)
    {
        Grid_IDO_IterationModified_w_f(Task,omega,force);
    }
}

void Grid_Cross_IterationSet    (Grid *Task, double omega, int N)
{
    int i;
    for(i=0;i<N;i++)
    {
        Grid_CrossIteration(Task,omega);
    }
}

void Grid_Cross_IterationSet_w_f(Grid *Task, double omega, Grid *force, int N)
{
    int i;
    for(i=0;i<N;i++)
    {
        Grid_CrossIteration_w_force(Task,omega,force);
    }
}

void Grid_Cross_IterationSet_w_f_w_autostop(Grid *Task, double omega, Grid *force, int N,
                                            double prev_diff, double *return_diff)
{
    int i,j,n;
    double diff;    //define a new array for u on (n+1) iteration
    double **nU;
    nU = malloc(Task->n * sizeof(double*));
    for(i=0; i<Task->n; i++)
        nU[i] = malloc(Task->n * sizeof(double));
    //for all values apply classical schematic

    for(n=0;n<N;n++)
    {
#pragma omp parallel for shared(nU,omega,Task,force) private(i,j)
        for(i=1; i<Task->n-1; i++)
        {
            for(j=1; j<Task->n-1; j++)
            {
                nU[i][j] = (1-omega)*Task->U[i][j]
                        +omega*0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                        +Task->U[i][j+1]+Task->U[i][j-1]
                        -force->U[i][j]*Task->h*Task->h);
            }
        }

        for(i=1; i<Task->n-1; i++)
            for(j=1; j<Task->n-1; j++)
                diff = fmax(diff, fabs(nU[i][j] - Task->U[i][j]));

        if(diff > prev_diff)
        {
            *return_diff = prev_diff;
            break;
        }
        else
            *return_diff = diff;

        for(i=1; i<Task->n-1; i++)
            for(j=1; j<Task->n-1; j++)
                Task->U[i][j] = nU[i][j];
    }
    printf("%10.10f\t",diff);


    for(i=0; i<Task->n; i++)
        free((void *)nU[i]);
    free((void *)nU);

}
