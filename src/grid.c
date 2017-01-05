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
    int i;
    Task->U      = malloc(n * sizeof(double*));
    Task->dUdx   = malloc(n * sizeof(double*));
    Task->dUdy   = malloc(n * sizeof(double*));
    for(i=0; i<n; i++)
    {
        Task->U[i]       = malloc(n * sizeof(double));
        Task->dUdx[i]    = malloc(n * sizeof(double));
        Task->dUdy[i]    = malloc(n * sizeof(double));
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

void Grid_CrossIteration(Grid *Task)
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
            nU[i][j] = Task->U[i][j] + Task->tau/Task->h/Task->h*
                    (Task->U[i+1][j]+Task->U[i-1][j]
                    +Task->U[i][j+1]+Task->U[i][j-1]
                     -4*Task->U[i][j]);
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
                Task->dUdx[i][j] =
                        0.5*(Task->U[i+1][j]-Task->U[i-1][j])/Task->h;
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
                Task->dUdx[i][j] =
                        0.5*(Task->U[i][j+1]-Task->U[i][j-1])/Task->h;
            }
            else
            {
                if(j==0)
                    Task->dUdx[i][0]=(Task->U[i][1]-Task->U[i][0])/Task->h;
                if(j==Task->n-1)
                    Task->dUdx[i][Task->n-1]=
                    (Task->U[i][Task->n-1]-Task->U[i][Task->n-2])/Task->h;
            }

        }
    }
}

void Grid_IDO_Iteration_with_f(Grid *Task, Grid *func)
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
                    + 2.*Task->tau*(d+D+func->U[i][j]);

            ndUdx[i][j] = Task->dUdx[i][j]
                    +Task->tau*(func->dUdx[i][j]
                               +6.*c+(Task->dUdx[i][j+1]
                                   -2*Task->dUdx[i][j]+
                                      Task->dUdx[i][j-1])/Task->h/Task->h);
            ndUdy[i][j] = Task->dUdy[i][j]
                    +Task->tau*(func->dUdy[i][j]
                               +6.*C+(Task->dUdy[i+1][j]
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
            if(i!=0 && i!=Task->n-1)
            {
                Target->dUdx[i][j] = ((*func)(x0+i*h+h,y0+j*h)
                                     -(*func)(x0+i*h-h,y0+j*h))/h/h*0.5;
            }
            else
            {
                if(i==0)
                    Target->dUdx[i][j] = ((*func)(x0+h,y0+j*h)
                                         -(*func)(x0,y0+j*h))/h/h;
                if(i==Task->n-1)
                    Target->dUdx[i][j] = ((*func)(x1,y0+j*h)
                                         -(*func)(x1-h,y0+j*h))/h/h;
            }


            if(j!=0 && j!=Task->n-1)
            {
                Target->dUdy[i][j] = ((*func)(x0+i*h,y0+j*h+h)
                                     -(*func)(x0+i*h,y0+j*h-h))/h/h*0.5;
            }
            else
            {
                if(j==0)
                    Target->dUdy[i][j] = ((*func)(x0+i*h,y0+h)
                                         -(*func)(x0+i*h,y0))/h/h;
                if(j==Task->n-1)
                    Target->dUdy[i][j] = ((*func)(x0+i*h,y1)
                                         -(*func)(x0+i*h,y1-h))/h/h;
            }
        }
    }
}
