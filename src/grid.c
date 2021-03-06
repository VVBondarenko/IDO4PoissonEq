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

void Grid_InitDirihlet_w_d(Grid *Task, double (*f)(double, double))
{
    //for all boundary nodes, set values form f to array
    int i; double h = 1./1024.;
    for(i = 0; i<Task->n; i++)
    {
        Task->U[i][0] = (*f)(Task->x0+Task->h*i,Task->y0);
        Task->U[i][Task->n-1]=(*f)(Task->x0+Task->h*i,Task->y1);

        Task->dUdx[i][0] = ((*f)(Task->x0+Task->h*i+h,Task->y0)
                           -(*f)(Task->x0+Task->h*i-h,Task->y0))*0.5/h;
        Task->dUdx[i][Task->n-1] = ((*f)(Task->x0+Task->h*i+h,Task->y1)
                                   -(*f)(Task->x0+Task->h*i-h,Task->y1))*0.5/h;


        Task->U[0][i] = (*f)(Task->x0,Task->y0+Task->h*i);
        Task->U[Task->n-1][i]=(*f)(Task->x1,Task->y0+Task->h*i);

        Task->dUdy[0][i] = ((*f)(Task->x0,Task->y0+Task->h*i+h)
                           -(*f)(Task->x0,Task->y0+Task->h*i-h))*0.5/h;
        Task->dUdy[Task->n-1][i] = ((*f)(Task->x1,Task->y0+Task->h*i+h)
                                   -(*f)(Task->x1,Task->y0+Task->h*i-h))*0.5/h;
    }
}

void Grid_Plot(Grid *Task)
{
    // for all values in U array, print values of its coords and U
    int i,j;
    FILE *output;

    output = fopen("Plot.dat", "w");
    for(i = 0; i < Task->n; i++)
    {
        for(j = 0; j < Task->n; j++)
            fprintf(output,"%15.15f %15.15f %15.15f\n",Task->x0+i*Task->h,
                                                       Task->y0+j*Task->h,
                                                       Task->U[i][j]);
            fprintf(output,"\n");
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
        fprintf(output,"\n");
    }
    printf("%10.10f\n",maxErr);
    fclose(output);
}

double Grid_print_error(Grid *Task, double (*exact)(double, double))
{
    // for all values in U array, print values of its coords and U
    int i,j;
    double errC = 0., errL1 = 0., errL2 = 0.,
            h = Task->h;
    for(i = 0; i < Task->n; i++)
    {
        for(j = 0; j < Task->n; j++)
        {
            errC = fmax(errC,fabs(Task->U[i][j]-(*exact)(Task->x0+i*h,Task->y0+j*h)));
            errL1 += fabs(Task->U[i][j]-(*exact)(Task->x0+i*h,Task->y0+j*h))*h*h;
            errL2 += pow(fabs(Task->U[i][j]-(*exact)(Task->x0+i*h,Task->y0+j*h)),2)*h*h;
        }
    }
//    printf("%d x %d & %10.10e & %10.10e & %10.10e \\\\ \n",Task->n,Task->n,errC, errL1, sqrt(errL2));
    return errL2;
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

void Grid_Copy(Grid *This, Grid *Source)
{
    int i,j;

    //ToDo: add check for allocated memory
    /*
    for(i=0; i<This->n; i++)
    {
        free((void *)This->U[i]);
        free((void *)This->dUdx[i]);
        free((void *)This->dUdy[i]);
        free((void *)This->dUdxdy[i]);
    }
    free((void *)This->U);
    free((void *)This->dUdx);
    free((void *)This->dUdy);
    free((void *)This->dUdxdy);
*/
    This->U = 0;
    This->dUdx = 0;
    This->dUdy = 0;
    This->dUdxdy=0;

    Grid_Init(This,Source->x0,Source->x1,Source->y0,Source->y1,Source->n,0.);

    for(i=0;i<Source->n;i++)
    {
        for(j=0;j<Source->n;j++)
        {
            This->U[i][j] = Source->U[i][j];
            This->dUdx[i][j] = Source->dUdx[i][j];
            This->dUdy[i][j] = Source->dUdy[i][j];
            This->dUdxdy[i][j]=Source->dUdxdy[i][j];
        }
    }
}

void Grid_Intensify(Grid *This, Grid *Source)
{
    This->U = 0;
    This->dUdx = 0;
    This->dUdy = 0;
    This->dUdxdy=0;

    Grid_Init(This,Source->x0,Source->x1,Source->y0,Source->y1,Source->n*2-1,0.);

    int i, j;
    double a, b, du, dup;

    for(i=0;i<Source->n;i++)
    {
        for(j=0;j<Source->n;j++)
        {
            This->U[2*i][2*j] = Source->U[i][j];
        }
    }

    for(i=0;i<Source->n-1;i++)
    {
        for(j=0;j<Source->n;j++)
        {
            a = (-2./Source->h*(Source->U[i+1][j]-Source->U[i][j])+(Source->dUdx[i+1][j]+Source->dUdx[i][j]))/Source->h/Source->h;
            b = (3./Source->h*(Source->U[i+1][j]-Source->U[i][j])-(Source->dUdx[i+1][j]+2.*Source->dUdx[i][j]))/Source->h;
            This->U[2*i+1][2*j] = Source->U[i][j]+Source->h/2.*(Source->dUdx[i][j]+Source->h/2.*(b+Source->h/2.*a));
        }
    }

    for(i=0;i<Source->n;i++)
    {
        for(j=0;j<Source->n-1;j++)
        {
            a = (-2./Source->h*(Source->U[i][j+1]-Source->U[i][j])+(Source->dUdy[i][j+1]+Source->dUdy[i][j]))/Source->h/Source->h;
            b = (3./Source->h*(Source->U[i][j+1]-Source->U[i][j])-(Source->dUdy[i][j+1]+2.*Source->dUdy[i][j]))/Source->h;
            This->U[2*i][2*j+1] = Source->U[i][j]+Source->h/2.*(Source->dUdy[i][j]+Source->h/2.*(b+Source->h/2.*a));
        }
    }

    for(i=0;i<Source->n-1;i++)
    {
        for(j=0;j<Source->n-1;j++)
        {
            du = (Source->dUdx[i][j]+Source->dUdy[i][j])/sqrt(2);
            dup= (Source->dUdx[i+1][j+1]+Source->dUdy[i+1][j+1])/sqrt(2);
            a = (-2./Source->h*(Source->U[i+1][j+1]-Source->U[i][j])+(dup+du))/Source->h/Source->h;
            b = (3./Source->h*(Source->U[i+1][j+1]-Source->U[i][j])-(dup+2.*du))/Source->h;
            This->U[2*i+1][2*j+1] = Source->U[i][j]+Source->h/2.*(du+Source->h/2.*(b+Source->h/2.*a));
        }
    }


}

void Grid_Intensify_v2(Grid *This, Grid *Source, double (*bf)(double, double))
{
    Grid_Intensify(This, Source);
    Grid_InitDirihlet(This,bf);
}
