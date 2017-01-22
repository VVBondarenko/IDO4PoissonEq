#include <cross.h>
#include <ido_original.h>

double boundary(double x, double y)
{
    return x*x-y*y;
}

double zero(double x, double y)
{
    return 0.*x*y;
}

int main(int argc, char **argv)
{

    int iters = atoi(argv[1]);
    Grid test_cross, test_ido, force;

    Grid_InitByFunction(&force,&zero,-M_PI, M_PI, -M_PI, M_PI, 64);
    Grid_Init (&test_cross, -M_PI, M_PI, -M_PI, M_PI, 64, 0.001);
    Grid_Init (&test_ido,   -M_PI, M_PI, -M_PI, M_PI, 64, 0.001);
//    GSIDO_Init(&test_ido,   -M_PI, M_PI, -M_PI, M_PI, 64);

    Grid_InitDirihlet    (&test_cross, &boundary);
    Grid_InitDirihlet_w_d(&test_ido,   &boundary);
//    GSIDO_InitDirihlet(&test_ido, &boundary);

    Cross_IterationSet(&test_cross,1.,iters);
    Grid_print_error(&test_cross, &boundary);

//    IDO_IterationOriginal_w_f(&test_ido,1.,&force);
    IDO_Ori_IterationSet_w_f_2(&test_ido,1.,&force,2*iters);
    Grid_print_error(&test_ido, &boundary);

    Grid_Plot_error(&test_ido,&boundary);
    Grid_Plot(&test_ido);
//    int i;
//    for(i=0;i<iters;i++)
//        GSIDO_early_Iteration_w_func(&test_ido,1.,zero);
//    Grid_print_error(&test_ido, &boundary);

    return 0;
}
