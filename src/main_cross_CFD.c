#include <cross.h>
/*
 * Tasks for the project structure
 *  todo: create UML diagram for this project
 *  todo: optimize classes
 *  todo: move cross and IDO solvers to their own classes
 *
 * Tasks for coding the basis system
 *  todo: implement SG-IDO
 *  todo: implement multigrid cross iterations
 *  todo: implement multigrid SG-IDO
 *
 * Tasks for coding final tests
 *  todo: implement output in format of accuracy tables and graphs
 *  todo: start using unit tests in development process
 */


double f(double t, double p)
{
    if(t<=0.)
        return 0.;
    if(t>0 && t<p)
        return t*t/p/p*(3.-2.*t/p);
    else
        return 1.;
}

double streamBoundary(double x, double y)
{
    return f(x-2./3.*y, 1./3.);
}

void findForce(Grid *force, Grid* stream, Grid* curl, double Reynolds)
{
    int i,j;

    for(i=1;i<force->n-1;i++)
    {
        for(j=1;j<force->n-1;j++)
        {
            force->U[i][j] = Reynolds*0.25*(
       (stream->U[i][j+1]-stream->U[i][j-1])*(curl->U[i+1][j]-curl->U[i-1][j])
      -(stream->U[i+1][j]-stream->U[i-1][j])*(curl->U[i][j+1]-curl->U[i][j-1])
                        )/force->h/force->h;
        }
    }
}

void curl_InitDirichlet(Grid *curl, Grid* stream)
{
    //for all boundary nodes, set values form f to array
    int i,n = curl->n;
    for(i = 1; i<curl->n-1; i++)
    {
        curl->U[i][0] = (stream->U[i][0]-2.*stream->U[i][1]+stream->U[i][2])/stream->h/stream->h;
        curl->U[i][n-1] = (stream->U[i][n-1]-2.*stream->U[i][n-2]+stream->U[i][n-3])/stream->h/stream->h;
        curl->U[0][i] = (stream->U[0][i]-2.*stream->U[1][i]+stream->U[2][i])/stream->h/stream->h;
        curl->U[n-1][i] = (stream->U[n-1][i]-2.*stream->U[n-2][i]+stream->U[n-3][i])/stream->h/stream->h;
    }
}


int main()
{
    omp_set_dynamic(1);
    omp_set_num_threads(8);
    double Re = 5., omega = 0.99;
    int Nodes = 64, IterQ = 5000;
    Grid curl, stream, force;
    Grid_Init(&curl,    0, 1, 0, 1, Nodes, 0.001);
    Grid_Init(&stream,  0, 1, 0, 1, Nodes, 0.001);
    Grid_Init(&force,  0, 1, 0, 1, Nodes, 0.001);
    Grid_InitDirihlet(&stream, &streamBoundary);


    Cross_IterationSet(&stream,1.,IterQ);
    curl_InitDirichlet(&curl,&stream);
    Cross_IterationSet(&curl,1.,IterQ);
    findForce(&force,&stream,&curl,Re);

    double streamDiff = 100., curlDiff = 100.;
    double prevStreamDiff = 1., prevCurlDiff = 100.;
    int i;
    for(i=0;i<30;i++)
    {
        printf("\n%d\t",i);
        Cross_IterationSet_w_f_w_autostop(&stream,omega,&curl,IterQ,1., &streamDiff);
        curl_InitDirichlet(&curl,&stream);
        Cross_IterationSet_w_f_w_autostop(&curl,omega,&force,IterQ, prevCurlDiff, &curlDiff);
        findForce(&force,&stream,&curl,Re);
        if(curlDiff == prevCurlDiff)
            break;
        else
            prevCurlDiff = curlDiff;
    }

    Grid_Plot(&stream);
    return 0;
}