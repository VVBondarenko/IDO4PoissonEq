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

void curl_InitDirihlet(Grid *curl, Grid* stream)
{
    //for all boundary nodes, set values form f to array
    int i,n = curl->n;
    for(i = 1; i<curl->n-1; i++)
    {
//        curl->U[i][0] = (stream->U[i][2]-2.*stream->U[i][1]-3.*stream->U[i][0]
//                +curl->U[i-1][0]-2.*curl->U[i][0]+curl->U[i+1][0])/curl->h/curl->h;
//        curl->U[i][curl->n-1]=(-stream->U[i][curl->n-3]+2.*stream->U[i][curl->n-2]+3.*stream->U[i][curl->n-1]
//                +curl->U[i-1][curl->n-1]-2.*curl->U[i][curl->n-1]+curl->U[i+1][curl->n-1])/curl->h/curl->h;

//        curl->U[0][i] = (stream->U[2][i]-2.*stream->U[1][i]-3.*stream->U[0][i]
//                +curl->U[0][i-1]-2.*curl->U[0][i]+curl->U[0][i+1])/curl->h/curl->h;
//        curl->U[curl->n-1][i]=(-stream->U[curl->n-3][i]+2.*stream->U[curl->n-2][i]+3.*stream->U[curl->n-1][i]
//                +curl->U[curl->n-1][i-1]-2.*curl->U[curl->n-1][i]+curl->U[curl->n-1][i+1])/curl->h/curl->h;
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

    double Re = 5.;
    int Nodes = 64;
    Grid curl, stream, force;
    Grid_Init(&curl,    0, 1, 0, 1, Nodes, 0.001);
    Grid_Init(&stream,  0, 1, 0, 1, Nodes, 0.001);
    Grid_Init(&force,  0, 1, 0, 1, Nodes, 0.001);
    Grid_InitDirihlet(&stream, &streamBoundary);


    Grid_Cross_IterationSet(&stream,1.,9000);
    curl_InitDirihlet(&curl,&stream);
    Grid_Cross_IterationSet(&curl,1.,9000);
    findForce(&force,&stream,&curl,Re);

    int i;
    for(i=0;i<150;i++)
    {
        printf("%d\n",i);
        Grid_Cross_IterationSet_w_f(&stream,1.,&curl,9000);
        curl_InitDirihlet(&curl,&stream);
        Grid_Cross_IterationSet_w_f(&curl,1.,&force,9000);
        findForce(&force,&stream,&curl,Re);
        if(i%10==0)
        {
            Grid_Plot(&stream);
            system("./lines_of_level &");
        }
    }

    Grid_Plot(&stream);
//    Grid_Plot_error(&stream, &boundary);
    return 0;
}
