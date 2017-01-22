#include <gs-ido.h>

void GSIDO_Init(Grid *Task,
                    double x0, double x1,
                    double y0, double y1,
                    int n)
{
    Task->x0 = x0;
    Task->x1 = x1;
    Task->y0 = y0;
    Task->y1 = y1;
    Task->n = n;
    Task->h = (x1-x0)/(n-1);

    Task->U      = malloc(n * sizeof(double*));
    Task->dUdx   = malloc((n-1) * sizeof(double*));
    Task->dUdy   = malloc(n * sizeof(double*));
    Task->dUdxdy = malloc((n-1) * sizeof(double*));

    int i;
    for(i=0; i<n; i++)
    {
        Task->U[i]       = malloc(n * sizeof(double));
        Task->dUdy[i]    = malloc((n-1) * sizeof(double));
        if(i<n-1)
        {
            Task->dUdx[i]    = malloc(n * sizeof(double));
            Task->dUdxdy[i]  = malloc((n-1) * sizeof(double));
        }
    }
}

void GSIDO_InitDirihlet(Grid *Task, double (*f)(double, double))
{
    int i, n = Task->n; double h = Task->h;
    for(i = 0; i<n; i++)
    {

        Task->U[i][0] = (*f)(Task->x0+h*i,Task->y0);
        Task->U[i][n-1]=(*f)(Task->x0+h*i,Task->y1);

        Task->U[0][i] = (*f)(Task->x0,Task->y0+h*i);
        Task->U[n-1][i]=(*f)(Task->x1,Task->y0+h*i);
    }

    for(i = 0; i<n-1; i++)
    {

        Task->dUdx[i][0] = ( (*f)(Task->x0+1.5*h+i*h, Task->y0+h/2)
                            -(*f)(Task->x0-0.5*h+i*h, Task->y0+h/2))
                            /h*0.5;
        Task->dUdx[i][n-1]=( (*f)(Task->x0+1.5*h+i*h, Task->y1-h/2)
                            -(*f)(Task->x0-0.5*h+i*h, Task->y1-h/2))
                            /h*0.5;


        Task->dUdy[0][i] = ( (*f)(Task->x0+h/2, Task->y0+1.5*h+i*h)
                            -(*f)(Task->x0+h/2, Task->y0-0.5*h+i*h))
                            /h*0.5;
        Task->dUdy[n-1][i]=( (*f)(Task->x1-h/2, Task->y0+1.5*h+i*h)
                            -(*f)(Task->x1-h/2, Task->y0-0.5*h+i*h))
                            /h*0.5;

    }

}


void GSIDO_early_Iteration_w_func(Grid *Task, double omega, double (*func)(double, double))
{
    int i,j, n = Task->n;
    double **nU, **ndUdx, **ndUdy, **ndUdxdy, h = Task->h;

    if(1)//allocating memory for temp arrays
    {
        nU      = malloc(n * sizeof(double*));
        ndUdx   = malloc((n-1) * sizeof(double*));
        ndUdy   = malloc(n * sizeof(double*));
        ndUdxdy = malloc((n-1) * sizeof(double*));

        for(i=0; i<n; i++)
        {
            nU[i]       = malloc(n * sizeof(double));
            ndUdy[i]    = malloc((n-1) * sizeof(double));
            if(i<n-1)
            {
                ndUdx[i]    = malloc(n * sizeof(double));
                ndUdxdy[i]  = malloc((n-1) * sizeof(double));
            }
        }
    }

    //iteration in U
    for(i = 1; i < n-1; i++)
    {
        for(j = 1; j< n-1; j++)
        {
//            nU[i][j] =
//            ((Task->U[i-1][j]+Task->U[i+1][j]+2*h*Task->dUdx[i-1][j]-2*h*Task->dUdx[i][j])/h/h
//            +(Task->U[i][j-1]+Task->U[i][j+1]+2*h*Task->dUdy[i][j-1]-2*h*Task->dUdy[i][j])/h/h
//                    + (*func)(Task->x0+i*h,Task->y0+j*h))*0.25*h*h;

            nU[i][j] = (-1/h/h*(Task->U[i+1][j]+Task->U[i-1][j])+2/h*(Task->dUdx[i][j]-Task->dUdx[i-1][j])
                        -1/h/h*(Task->U[i][j+1]+Task->U[i][j-1])+2/h*(Task->dUdy[i][j]-Task->dUdy[i][j-1])
                        -(*func)(Task->x0+i*h,Task->y0+j*h) )*(-0.5)*h*h;

        }
    }

//    printf("starting boundary initialization\n");
    //prepare boundary for dUdx and dUdy

    for(i=0;i<n-1;i++)
    {
        Task->dUdx[0][i] =      (Task->U[1][i]-Task->U[0][i])/h*0.5;
        Task->dUdx[n-2][i] =    (Task->U[n-1][i]-Task->U[n-2][i])/h*0.5;
        Task->dUdy[i][0] =      (Task->U[i][1]-Task->U[i][0])/h*0.5;
        Task->dUdy[i][n-2] =    (Task->U[i][n-1]-Task->U[i][n-2])/h*0.5;
        //and dUdxdy
        Task->dUdxdy[0][i] =    (Task->U[1][i+1]-Task->U[1][i]-Task->U[0][i+1]+Task->U[0][i])/h/h;
        Task->dUdxdy[n-2][i] =  (Task->U[n-1][i+1]-Task->U[n-1][i]-Task->U[n-2][i+1]+Task->U[n-2][i])/h/h;
        Task->dUdxdy[i][0] =    (Task->U[i+1][1]-Task->U[i+1][0]-Task->U[i][1]+Task->U[i][0])/h/h;
        Task->dUdxdy[i][n-2] =  (Task->U[i+1][n-1]-Task->U[i+1][n-2]-Task->U[i][n-1]+Task->U[i][n-2])/h/h;
    }


    //iteration in dUdx
    for(i=1;i<n-2;i++)
        for(j=1;j<n-1;j++)
//            ndUdx[i][j] = (1-omega)*Task->dUdx[i][j]
//                    +omega*0.25*(Task->dUdx[i+1][j]+Task->dUdx[i-1][j]
//                    +Task->dUdx[i][j+1]+Task->dUdx[i][j-1]
//                    -((*func)(Task->x0+h*1.5+i*h,Task->y0+j*h)-(*func)(Task->x0-0.5*h+i*h,Task->y0+j*h))*0.5*Task->h);
            ndUdx[i][j] = (24/h/h/h*(Task->U[i+1][j]-Task->U[i][j])
                    -1./h/h*(Task->dUdx[i][j+1]+Task->dUdx[i][j-1])
                    +2/h*(Task->dUdxdy[i][j]-Task->dUdxdy[i][j-1])
                    -((*func)(Task->x0+h*1.5+i*h,Task->y0+j*h)-(*func)(Task->x0-0.5*h+i*h,Task->y0+j*h))*0.5/Task->h)/22*h*h;

    //iteration in dUdy
    for(i=1;i<n-1;i++)
        for(j=1;j<n-2;j++)
//            ndUdy[i][j] = (1-omega)*Task->dUdy[i][j]
//                    +omega*0.25*(Task->dUdy[i+1][j]+Task->dUdx[i-1][j]
//                    +Task->dUdy[i][j+1]+Task->dUdx[i][j-1]
//                    -((*func)(Task->x0+i*h,Task->y0+1.5*h+j*h)-(*func)(Task->x0+i*h,Task->y0-0.5*h+j*h))*0.5*Task->h);
            ndUdy[i][j] = (-1/h/h*(Task->dUdy[i+1][j]+Task->dUdy[i-1][j])
                    +2/h*(Task->dUdxdy[i][j]-Task->dUdxdy[i-1][j])+
                    24/h/h/h*(Task->U[i][j+1]-Task->U[i][j])
                    -((*func)(Task->x0+i*h,Task->y0+1.5*h+j*h)-(*func)(Task->x0+i*h,Task->y0-0.5*h+j*h))*0.5/Task->h)/22*h*h;

    //iteration in dudxdy
    for(i=1;i<n-2;i++)
        for(j=1;j<n-2;j++)
            ndUdxdy[i][j] = (24/h/h/h*(Task->dUdy[i+1][j]-Task->dUdy[i][j])
                    +24/h/h/h*(Task->dUdx[i][j+1]-Task->dUdx[i][j])
                    -((*func)(Task->x0+i*h+1.5*h,Task->y0+1.5*h+j*h)-(*func)(Task->x0+i*h+1.5*h,Task->y0-0.5*h+j*h)
                     -(*func)(Task->x0+i*h-0.5*h,Task->y0+1.5*h+j*h)+(*func)(Task->x0+i*h-0.5*h,Task->y0-0.5*h+j*h))*0.25/h/h
                    )/48*h*h;

//        dfdxdy = 0.25/h/h*(f(x+h,y+h)-f(x+h,y-h)-f(x-h,y+h)+f(x-h,y-h));
//        ((*func)(Task->x0+i*h+1.5*h,Task->y0+1.5*h+j*h)-(*func)(Task->x0+i*h+1.5*h,Task->y0-0.5*h+j*h)
//        -(*func)(Task->x0+i*h-0.5*h,Task->y0+1.5*h+j*h)+(*func)(Task->x0+i*h-0.5*h,Task->y0-0.5*h+j*h))*0.25/h/h

    //copying
    for(i=1;i<n-1;i++)
        for(j=1;j<n-1;j++)
            Task->U[i][j] = nU[i][j];

    for(i=1;i<n-2;i++)
        for(j=1;j<n-1;j++)
            Task->dUdx[i][j] = ndUdx[i][j];

    for(i=1;i<n-1;i++)
        for(j=1;j<n-2;j++)
            Task->dUdy[i][j] = ndUdy[i][j];

    for(i=1;i<n-2;i++)
        for(j=1;j<n-2;j++)
            Task->dUdxdy[i][j] = ndUdxdy[i][j];

    //clean up memory

}

