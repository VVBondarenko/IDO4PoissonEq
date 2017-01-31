//
// Created by vbondarenko on 19.01.17.
//

#include <cross.h>

void Cross_Iteration(Grid *Task, double omega)
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

void Cross_Iteration_w_force(Grid *Task, double omega, Grid *force)
{
    int i,j;
    //define a new array for u on (n+1) iteration
    double **nU;
    nU = malloc(Task->n * sizeof(double*));
    for(i=0; i<Task->n; i++)
        nU[i] = malloc(Task->n * sizeof(double));
    //for all values apply classical schematic
//#pragma omp parallel for shared(nU,omega,Task,force) private(i,j)
#pragma omp simd
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

void Cross_IterationSet    (Grid *Task, double omega, int N)
{
    int i;
    for(i=0;i<N;i++)
    {
        Cross_Iteration(Task,omega);
    }
}

void Cross_IterationSet_w_f(Grid *Task, double omega, Grid *force, int N)
{
    int i,j, n=Task->n;
    for(i=0;i<N;i++)
    {
        Cross_Iteration_w_force(Task,omega,force);
    }

    double h = Task->h;
    for(i=1;i<n-1;i++)
    {
        for(j=1;j<n-1;j++)
        {
            Task->dUdx[i][j] = (Task->U[i+1][j]-Task->U[i-1][j])*0.5/h;
            Task->dUdy[i][j] = (Task->U[i][j+1]-Task->U[i][j-1])*0.5/h;
            Task->dUdxdy[i][j] = 0.25/h/h*(Task->U[i+1][j+1]-Task->U[i-1][j+1]
                    -Task->U[i+1][j-1]+Task->U[i-1][j-1]);
        }
    }
}

void Cross_IterationSet_w_f_w_autostop(Grid *Task, double omega, Grid *force, int N,
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
//#pragma omp parallel for shared(nU,omega,Task,force) private(i,j)
#pragma omp simd
        for (i = 1; i < Task->n - 1; i++)
        {
            for (j = 1; j < Task->n - 1; j++)
            {
                nU[i][j] = (1 - omega) * Task->U[i][j]
                           + omega * 0.25 * (Task->U[i + 1][j] + Task->U[i - 1][j]
                                             + Task->U[i][j + 1] + Task->U[i][j - 1]
                                             - force->U[i][j] * Task->h * Task->h);
            }
        }

//#pragma omp simd
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
//#pragma omp simd
        for(i=1; i<Task->n-1; i++)
            for(j=1; j<Task->n-1; j++)
                Task->U[i][j] = nU[i][j];
    }
    printf("%10.10f\t",diff);


    for(i=0; i<Task->n; i++)
        free((void *)nU[i]);
    free((void *)nU);

}
