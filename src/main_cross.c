//
// Created by vbondarenko on 19.01.17.
//

#include <cross.h>

double boundary(double x, double y)
{
    return x*x-y*y;
}

int main()
{
    Grid test;
    Grid_Init(&test, -M_PI, M_PI, -M_PI, M_PI, 64, 0.001);
    Grid_InitDirihlet_w_derivatives(&test, &boundary);

    Cross_IterationSet(&test,1.,5000);
    Grid_Plot(&test);
    Grid_Plot_error(&test, &boundary);

    return 0;
}
