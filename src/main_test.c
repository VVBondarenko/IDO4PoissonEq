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

double f2(double x, double y)
{
    return -2.*sin(x*M_PI)*sin(y*M_PI)*M_PI*M_PI;
}

double e2(double x, double y)
{
    return sin(x*M_PI)*sin(y*M_PI);
}

double f3(double x, double y)
{
    return 2*exp((-1 + x*x)*(-1 + y*y))*
            (-2 + 3*y*y + 2*x*x*x*x*y*y +
              x*x*(3 - 8*y*y + 2*y*y*y*y));
}

double e3(double x, double y)
{
    return exp((x*x-1)*(y*y-1))-1.;
}

int test(int argc, char **argv)
{

    int iters = atoi(argv[1]);
    int sizeN = 4;
    int sizes[] = {8, 16, 32, 64};
    int i;
    printf("test 1: laplace, poly pow-2\n");
    printf("cross\n");
    for(i = 0; i < sizeN; i++)
    {
        Grid test_cross;
        Grid_Init (&test_cross, -M_PI, M_PI, -M_PI, M_PI, sizes[i], 0.001);

        Grid_InitDirihlet    (&test_cross, &boundary);

        Cross_IterationSet(&test_cross,1.,iters);
        Grid_print_error(&test_cross, &boundary);
    }
    printf("IDO:\n");
    for(i = 0; i < sizeN; i++)
    {
        Grid test_ido, zero_force;

        Grid_InitByFunction(&zero_force,&zero,-M_PI, M_PI, -M_PI, M_PI, sizes[i]);
        Grid_Init               (&test_ido,   -M_PI, M_PI, -M_PI, M_PI, sizes[i], 0.001);

        Grid_InitDirihlet_w_d(&test_ido,   &boundary);

        IDO_Ori_IterationSet_w_f_Seidel(&test_ido,1.,&zero_force,iters);
        Grid_print_error(&test_ido, &boundary);
    }

    printf("test 2: jap. sin*sin\n");
    printf("cross:\n");
    for(i = 0; i < sizeN; i++)
    {
        Grid test_cross2, force2;

        Grid_InitByFunction(&force2,&f2,-M_PI, M_PI, -M_PI, M_PI, sizes[i]);
        Grid_Init (&test_cross2,        -M_PI, M_PI, -M_PI, M_PI, sizes[i], 0.001);

        Grid_InitDirihlet    (&test_cross2, &zero);

        Cross_IterationSet_w_f(&test_cross2,1.,&force2,iters);
        Grid_print_error(&test_cross2, &e2);
    }

    printf("IDO:\n");
    for(i = 0; i < sizeN; i++)
    {

        Grid test_ido2, force2;

        Grid_InitByFunction(&force2,&f2,-M_PI, M_PI, -M_PI, M_PI, sizes[i]);
        Grid_Init (&test_ido2,          -M_PI, M_PI, -M_PI, M_PI, sizes[i], 0.001);

        Grid_InitDirihlet_w_d(&test_ido2,   &zero);

        IDO_Ori_IterationSet_w_f_Seidel(&test_ido2,1.,&force2,iters);
        Grid_print_error(&test_ido2, &e2);
    }

    printf("test 3: exp-hat\n");
    printf("cross:\n");
    for(i = 0; i < sizeN; i++)
    {
        Grid test_cross3, force3;

        Grid_InitByFunction(&force3,&f3,-1, 1, -1, 1, sizes[i]);
        Grid_Init (&test_cross3,        -1, 1, -1, 1, sizes[i], 0.001);

        Grid_InitDirihlet    (&test_cross3, &zero);

        Cross_IterationSet_w_f(&test_cross3,1.,&force3,iters);
        Grid_print_error(&test_cross3, &e3);
    }
    printf("IDO:\n");
    for(i = 0; i < sizeN; i++)
    {
        Grid test_ido3, force3;

        Grid_InitByFunction(&force3,&f3,-1, 1, -1, 1, sizes[i]);
        Grid_Init (&test_ido3,          -1, 1, -1, 1, sizes[i], 0.001);

        Grid_InitDirihlet_w_d(&test_ido3,   &zero);

        IDO_Ori_IterationSet_w_f_Seidel(&test_ido3,1.,&force3,iters);
        Grid_print_error(&test_ido3, &e3);
        Grid_Plot(&test_ido3);
    }


//    Grid_Plot_error(&test_ido,&boundary);
//    Grid_Plot(&test_ido2);
//    int i;
//    for(i=0;i<iters;i++)
//        GSIDO_early_Iteration_w_func(&test_ido,1.,zero);
//    Grid_print_error(&test_ido, &boundary);

    return 0;
}


int main()
{
    FILE *err_tab;
    err_tab = fopen("err_of_omega.dat","w");

    int i,N = 129;
    double omega = 0.5;

    Grid initial_cross, force3;

    Grid_InitByFunction(&force3,&f2,-1, 1, -1, 1, N);
    Grid_Init (&initial_cross,      -1, 1, -1, 1, N, 0.001);

    Grid_InitDirihlet_w_d(&initial_cross,   &e2);

    Cross_IterationSet_w_f(&initial_cross,1.,&force3, 10000);
//        Grid_print_error(&test_ido3,&e3);
    IDO_InitNormalD_2(&initial_cross,&force3);


    for(omega; omega<2.;omega+=0.015)
    {
        Grid test_ido3;

        Grid_Copy(&test_ido3,&initial_cross);

        for(i=0;i<20;i++)
        {
            IDO_Ori_IterationSet_w_f_Seidel(&test_ido3,omega,&force3,1);
            fprintf(err_tab,"%d %f %f\n",i,omega,log10(Grid_print_error(&test_ido3,&e2)/4.));
//            Grid_print_error(&test_ido3,&e3);
        }
        fprintf(err_tab,"\n");
    }

    return 0;
}
