#include <grid.h>


double boundary(double x, double y)
{
//    return cos(x)*cos(y);
    return x*x-y*y;
}

int tester()
{
    Grid test;
    Grid_Init(&test, -M_PI, M_PI, -M_PI, M_PI, 64, 0.001);
//    Grid_InitDirihlet(&test, &boundary); //checked
    Grid_InitDirihlet_w_derivatives(&test, &boundary);
//    Grid_IDO_InitDeriv(&test);

    int i;
    for(i=0;i<5000;i++)
    {
        if(i%5==0)
        {
            printf("%d\t", i);
            Grid_print_error(&test,&boundary);
        }
//        Grid_CrossIteration(&test,1.);
//        Grid_IDO_IterationModified(&test,1.);
        Grid_IDO_IterationOriginal(&test,1.);
    }
    Grid_Plot(&test);
    Grid_Plot_error(&test, &boundary);

    return 0;
}

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
    return f(x-4./5.*y, 1./5.);
}

void findForce(Grid *force, Grid* stream, Grid* curl, double Reynolds)
{
    int i,j;

//    double stream_xx, stream_yy, curl_xx, curl_yy;
    for(i=1;i<force->n-1;i++)
    {
        for(j=1;j<force->n-1;j++)
        {
//            stream_xx = (1./stream->h*(stream->U[i+1][j]
//                                   -2.*stream->U[i][j]
//                                      +stream->U[i-1][j])
//                    -0.25*(stream->dUdx[i+1][j]-stream->dUdx[i-1][j]))/stream->h;
//            stream_yy = (1./stream->h*(stream->U[i][j+1]
//                                   -2.*stream->U[i][j]
//                                      +stream->U[i][j-1])
//                    -0.25*(stream->dUdy[i][j+1]-stream->dUdy[i][j-1]))/stream->h;
//            curl_xx = (1./curl->h*(curl->U[i+1][j]
//                               -2.*curl->U[i][j]
//                                  +curl->U[i-1][j])
//                    -0.25*(curl->dUdx[i+1][j]-curl->dUdx[i-1][j]))/curl->h;
//            curl_yy = (1./curl->h*(curl->U[i][j+1]
//                               -2.*curl->U[i][j]
//                                  +curl->U[i][j-1])
//                    -0.25*(curl->dUdy[i][j+1]-curl->dUdy[i][j-1]))/curl->h;

            force->U[i][j] = Reynolds*0.25*(
       (stream->U[i][j+1]-stream->U[i][j-1])*(curl->U[i+1][j]-curl->U[i-1][j])
      -(stream->U[i+1][j]-stream->U[i-1][j])*(curl->U[i][j+1]-curl->U[i][j-1])
                        )/force->h/force->h;

//            force->dUdx[i][j] = Reynolds*(
//                       stream->dUdxdy[i][j]*curl->dUdx[i][j]+stream->dUdy[i][j]*curl_xx
//                     -(stream_xx*curl->dUdy[i][j]+stream->dUdx[i][j]*curl->dUdxdy[i][j]));

//            force->dUdy[i][j] = Reynolds*(
//                        stream_yy*curl->dUdx[i][j]+stream->dUdy[i][j]*curl->dUdxdy[i][j]
//                     -(stream->dUdxdy[i][j]*curl->dUdy[i][j]+stream->dUdx[i][j]*curl_yy));

//            force->dUdxdy[i][j] = Reynolds*(
//                        0.5/stream->h*(stream->dUdxdy[i][j+1]-stream->dUdxdy[i][j-1])*curl->dUdx[i][j]
//                       +curl_xx*stream_yy+stream->dUdxdy[i][j]*curl->dUdxdy[i][j]+
//                       +stream->dUdy[i][j]*0.5/curl->h*(curl->dUdxdy[i+1][j]-curl->dUdxdy[i-1][j])
//                    -(0.5/stream->h*(stream->dUdxdy[i+1][j]-stream->dUdxdy[i-1][j])*curl->dUdy[i][j]
//                    +curl_yy*stream_xx+stream->dUdxdy[i][j]*curl->dUdxdy[i][j]+
//                    +stream->dUdx[i][j]*0.5/curl->h*(curl->dUdxdy[i][j+1]-curl->dUdxdy[i][j-1])));
        }
    }
}


void findForce_IDO(Grid *force, Grid* stream, Grid* curl, double Reynolds)
{
    int i,j;

    double stream_xx, stream_yy, curl_xx, curl_yy;
    for(i=1;i<force->n-1;i++)
    {
        for(j=1;j<force->n-1;j++)
        {
            stream_xx = (1./stream->h*(stream->U[i+1][j]
                                   -2.*stream->U[i][j]
                                      +stream->U[i-1][j])
                    -0.25*(stream->dUdx[i+1][j]-stream->dUdx[i-1][j]))/stream->h;
            stream_yy = (1./stream->h*(stream->U[i][j+1]
                                   -2.*stream->U[i][j]
                                      +stream->U[i][j-1])
                    -0.25*(stream->dUdy[i][j+1]-stream->dUdy[i][j-1]))/stream->h;
            curl_xx = (1./curl->h*(curl->U[i+1][j]
                               -2.*curl->U[i][j]
                                  +curl->U[i-1][j])
                    -0.25*(curl->dUdx[i+1][j]-curl->dUdx[i-1][j]))/curl->h;
            curl_yy = (1./curl->h*(curl->U[i][j+1]
                               -2.*curl->U[i][j]
                                  +curl->U[i][j-1])
                    -0.25*(curl->dUdy[i][j+1]-curl->dUdy[i][j-1]))/curl->h;

            force->U[i][j] = Reynolds*(
                        stream->dUdy[i][j]*curl->dUdx[i][j]
                       -stream->dUdx[i][j]*curl->dUdy[i][j]);

            force->dUdx[i][j] = Reynolds*(
                       stream->dUdxdy[i][j]*curl->dUdx[i][j]+stream->dUdy[i][j]*curl_xx
                     -(stream_xx*curl->dUdy[i][j]+stream->dUdx[i][j]*curl->dUdxdy[i][j]));

            force->dUdy[i][j] = Reynolds*(
                        stream_yy*curl->dUdx[i][j]+stream->dUdy[i][j]*curl->dUdxdy[i][j]
                     -(stream->dUdxdy[i][j]*curl->dUdy[i][j]+stream->dUdx[i][j]*curl_yy));

            force->dUdxdy[i][j] = Reynolds*(
                        0.5/stream->h*(stream->dUdxdy[i][j+1]-stream->dUdxdy[i][j-1])*curl->dUdx[i][j]
                       +curl_xx*stream_yy+stream->dUdxdy[i][j]*curl->dUdxdy[i][j]+
                       +stream->dUdy[i][j]*0.5/curl->h*(curl->dUdxdy[i+1][j]-curl->dUdxdy[i-1][j])
                    -(0.5/stream->h*(stream->dUdxdy[i+1][j]-stream->dUdxdy[i-1][j])*curl->dUdy[i][j]
                    +curl_yy*stream_xx+stream->dUdxdy[i][j]*curl->dUdxdy[i][j]+
                    +stream->dUdx[i][j]*0.5/curl->h*(curl->dUdxdy[i][j+1]-curl->dUdxdy[i][j-1])));
        }
    }
}


void curl_InitDirihlet(Grid *curl, Grid* stream)
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

int stream_curl_by_cross()
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


    Grid_Cross_IterationSet(&stream,1.,IterQ);
    curl_InitDirihlet(&curl,&stream);
    Grid_Cross_IterationSet(&curl,1.,IterQ);
    findForce(&force,&stream,&curl,Re);

    double streamDiff = 100., curlDiff = 100.;
    double prevStreamDiff = 1., prevCurlDiff = 100.;
    int i;
    for(i=0;i<30;i++)
    {
        printf("\n%d\t",i);
        Grid_Cross_IterationSet_w_f_w_autostop(&stream,omega,&curl,IterQ,1., &streamDiff);
        curl_InitDirihlet(&curl,&stream);
        Grid_Cross_IterationSet_w_f_w_autostop(&curl,omega,&force,IterQ, prevCurlDiff, &curlDiff);
        findForce(&force,&stream,&curl,Re);
        if(curlDiff == prevCurlDiff)
            break;
        else
            prevCurlDiff = curlDiff;
    }

    Grid_Plot(&stream);
    return 0;
}

int main()
{
    omp_set_dynamic(1);
    omp_set_num_threads(8);

    double Re = 1., omega = 0.8;
    int Nodes = 64, IterQ = 7000;
    Grid curl, stream, force;
    Grid_Init(&curl,    0, 1, 0, 1, Nodes, 0.001);
    Grid_Init(&stream,  0, 1, 0, 1, Nodes, 0.001);
    Grid_Init(&force,  0, 1, 0, 1, Nodes, 0.001);
    Grid_InitDirihlet(&stream, &streamBoundary);


    Grid_IDO_IterationSet(&stream,1.,IterQ);
    curl_InitDirihlet(&curl,&stream);
    Grid_IDO_IterationSet(&curl,1.,IterQ);
    findForce_IDO(&force,&stream,&curl,Re);

    double streamDiff = 100., curlDiff = 100.;
    double prevStreamDiff = 1., prevCurlDiff = 100.;
    int i;
    for(i=0;i<30;i++)
    {
        printf("\n%d\t",i);
        Grid_IDO_IterationSet_w_f(&stream,omega,&curl,IterQ);
        curl_InitDirihlet(&curl,&stream);
        Grid_IDO_IterationSet_w_f(&curl,omega,&force,IterQ);
        findForce_IDO(&force,&stream,&curl,Re);

        if(i%5==0)
        {
            Grid_Plot(&stream);
            system("./lines_of_level");
        }
//        if(curlDiff == prevCurlDiff)
//            break;
//        else
//            prevCurlDiff = curlDiff;
    }

    Grid_Plot(&stream);
    return 0;
}
