#include <grid.h>


double boundary(double x, double y)
{
    return cos(x)*cos(y);
}

int main()
{
    Grid test;
    Grid_Init(&test, -M_PI, M_PI, -M_PI, M_PI, 64, 0.001);
    Grid_InitDirihlet(&test, &boundary); //checked
    Grid_IDO_InitDeriv(&test);

    int i;
    for(i=0;i<1000;i++)
    {
//        Grid_CrossIteration(&test); //checked
        Grid_IDO_Iteration(&test);
    }
    Grid_Plot(&test);


    return 0;
}
