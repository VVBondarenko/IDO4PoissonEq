#include <grid.h>


double boundary(double x, double y)
{
    return cos(x)*cos(y);
}

int main()
{
    Grid test;
    Grid_Init(&test, -M_PI, M_PI, -M_PI, M_PI, 16, 0.01);
    Grid_InitDirihlet(&test, &boundary); //checked

    int i;
    for(i=0;i<400;i++)
    {
        Grid_CrossIteration(&test); //checked
    }
    Grid_Plot(&test);


    return 0;
}
