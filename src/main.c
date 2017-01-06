#include <grid.h>


double boundary(double x, double y)
{
//    return cos(x)*cos(y);
    return x*x-y*y;
}

int main()
{
    Grid test;
    Grid_Init(&test, -M_PI, M_PI, -M_PI, M_PI, 16, 0.001);
//    Grid_InitDirihlet(&test, &boundary); //checked
    Grid_InitDirihlet_w_derivatives(&test, &boundary);
//    Grid_IDO_InitDeriv(&test);

    int i;
    for(i=0;i<500;i++)
    {
//        Grid_CrossIteration(&test); //checked
//        if(i%100==0)
//            Grid_IDO_InitDeriv(&test);
        if(i%10==0)Grid_print_error(&test,&boundary);
        Grid_IDO_Iteration_2(&test,1.01); //checked
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


//int main()
//{

//    return 0;
//}
