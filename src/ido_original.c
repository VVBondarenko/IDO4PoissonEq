//
// Created by vbondarenko on 19.01.17.
//

#include <ido_original.h>

void IDO_InitDeriv(Grid *Task)
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

void IDO_Iteration(Grid *Task)
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

void IDO_IterationOriginal(Grid *Task, double omega)
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

void IDO_IterationOriginal_w_f(Grid *Task, double omega, Grid *force)
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
            nU[i][j] = (1-omega)*Task->U[i][j]
                       +omega*(2/h/h*(Task->U[i+1][j]+Task->U[i-1][j]
                                      +Task->U[i][j+1]+Task->U[i][j-1])
                               -0.5/h*(Task->dUdx[i+1][j]-Task->dUdx[i-1][j]
                                       +Task->dUdy[i][j+1]-Task->dUdy[i][j-1])-force->U[i][j])/8.*h*h;

            ndUdx[i][j] = (1-omega)*Task->dUdx[i][j]
                          +omega*(7.5*(Task->U[i+1][j]-Task->U[i-1][j])/h/h/h
                                  -1.5*(Task->dUdx[i+1][j]+Task->dUdx[i-1][j])/h/h
                                  +2.*(Task->dUdx[i][j+1]+Task->dUdx[i][j-1])/h/h
                                  -0.5*(Task->dUdxdy[i][j+1]-Task->dUdxdy[i][j-1])/h-force->dUdx[i][j])
                           /(12/h/h+4/h/h);
            ndUdy[i][j] = (1-omega)*Task->dUdy[i][j]
                          +omega*(7.5*(Task->U[i][j+1]-Task->U[i][j-1])/h/h/h
                                  -1.5*(Task->dUdy[i][j+1]+Task->dUdy[i][j-1])/h/h
                                  +2.*(Task->dUdy[i+1][j]+Task->dUdy[i-1][j])/h/h
                                  -0.5*(Task->dUdxdy[i+1][j]-Task->dUdxdy[i-1][j])/h-force->dUdy[i][j])
                           /(12/h/h+4/h/h);
            ndUdxdy[i][j] = (1-omega)*Task->dUdxdy[i][j]
                            +omega*(7.5*(Task->dUdx[i][j+1]-Task->dUdx[i][j-1]
                                         +Task->dUdy[i+1][j]-Task->dUdy[i-1][j])/h/h/h
                                    -1.5*(Task->dUdxdy[i][j+1]+Task->dUdxdy[i][j-1]
                                          +Task->dUdxdy[i+1][j]+Task->dUdxdy[i-1][j])/h/h
                                    -force->dUdxdy[i][j]
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

void IDO_IterationModified(Grid *Task, double omega)
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

void IDO_IterationModified_w_f(Grid *Task, double omega, Grid *force)
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
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j])/h-force->dUdx[i][j];
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1])/h-force->dUdy[i][j];
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                              -Task->U[i-1][j+1]
                              -Task->U[i+1][j-1]
                              +Task->U[i-1][j-1])/h/h-force->dUdxdy[i][j];
        i=p;j=Task->n-2;
        nU[i][j] = 0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                         +Task->U[i][j+1]+Task->U[i][j-1]-force->U[i][j]*Task->h*Task->h);
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j])/h-force->dUdx[i][j];
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1])/h-force->dUdy[i][j];
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                              -Task->U[i-1][j+1]
                              -Task->U[i+1][j-1]
                              +Task->U[i-1][j-1])/h/h-force->dUdxdy[i][j];
        i=1;j=p;
        nU[i][j] = 0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                         +Task->U[i][j+1]+Task->U[i][j-1]-force->U[i][j]*Task->h*Task->h);
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j])/h-force->dUdx[i][j];
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1])/h-force->dUdy[i][j];
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                              -Task->U[i-1][j+1]
                              -Task->U[i+1][j-1]
                              +Task->U[i-1][j-1])/h/h-force->dUdxdy[i][j];
        i=Task->n-2;j=p;
        nU[i][j] = 0.25*(Task->U[i+1][j]+Task->U[i-1][j]
                         +Task->U[i][j+1]+Task->U[i][j-1]-force->U[i][j]*Task->h*Task->h);
        ndUdx[i][j] = 0.5*(Task->U[i+1][j]-Task->U[i-1][j])/h-force->dUdx[i][j];
        ndUdy[i][j] = 0.5*(Task->U[i][j+1]-Task->U[i][j-1])/h-force->dUdy[i][j];
        ndUdxdy[i][j] = 0.25*(Task->U[i+1][j+1]
                              -Task->U[i-1][j+1]
                              -Task->U[i+1][j-1]
                              +Task->U[i-1][j-1])/h/h-force->dUdxdy[i][j];

    }

    for(i=2; i<Task->n-2; i++)
    {
        for(j=2; j<Task->n-2; j++)
        {

            nU[i][j] = (1-omega)*Task->U[i][j]
                       +omega*(2/h/h*(Task->U[i+1][j]+Task->U[i-1][j]
                                      +Task->U[i][j+1]+Task->U[i][j-1])
                               -0.5/h*(Task->dUdx[i+1][j]-Task->dUdx[i-1][j]
                                       +Task->dUdy[i][j+1]-Task->dUdy[i][j-1])-force->U[i][j])/8.*h*h;

            ndUdx[i][j] = (1-omega)*Task->dUdx[i][j]
                          +omega*(7.5*(Task->U[i+1][j]-Task->U[i-1][j])/h/h/h
                                  -1.5*(Task->dUdx[i+1][j]+Task->dUdx[i-1][j])/h/h
                                  +2.*(Task->dUdx[i][j+1]+Task->dUdx[i][j-1])/h/h
                                  -0.5*(Task->dUdxdy[i][j+1]-Task->dUdxdy[i][j-1])/h-force->dUdx[i][j])
                           /(12/h/h+4/h/h);
            ndUdy[i][j] = (1-omega)*Task->dUdy[i][j]
                          +omega*(7.5*(Task->U[i][j+1]-Task->U[i][j-1])/h/h/h
                                  -1.5*(Task->dUdy[i][j+1]+Task->dUdy[i][j-1])/h/h
                                  +2.*(Task->dUdy[i+1][j]+Task->dUdy[i-1][j])/h/h
                                  -0.5*(Task->dUdxdy[i+1][j]-Task->dUdxdy[i-1][j])/h-force->dUdy[i][j])
                           /(12/h/h+4/h/h);
            ndUdxdy[i][j] = (1-omega)*Task->dUdxdy[i][j]
                            +omega*(7.5*(Task->dUdx[i][j+1]-Task->dUdx[i][j-1]
                                         +Task->dUdy[i+1][j]-Task->dUdy[i-1][j])/h/h/h
                                    -1.5*(Task->dUdxdy[i][j+1]+Task->dUdxdy[i][j-1]
                                          +Task->dUdxdy[i+1][j]+Task->dUdxdy[i-1][j])/h/h
                                    -force->dUdxdy[i][j])/24*h*h;
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

void IDO_IterationSet      (Grid *Task, double omega, int N)
{
    int i;
    for(i=0;i<N;i++)
    {
        IDO_IterationModified(Task,omega);
    }
}

void IDO_Mod_IterationSet_w_f  (Grid *Task, double omega, Grid *force, int N)
{
    int i;
    for(i=0;i<N;i++)
    {
        IDO_IterationModified_w_f(Task,omega,force);
    }
}

void IDO_Ori_IterationSet_w_f  (Grid *Task, double omega, Grid *force, int N)
{
    int i;
    for(i=0;i<N;i++)
    {
        IDO_IterationOriginal_w_f(Task,omega,force);
    }
}
