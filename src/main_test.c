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

/*========================================
 *
 * Plotting colormap graph for accuracy,
 * depending of relaxation parameter and
 * iteration numbers.
 *
 */
int main_accur_plot()
{
    FILE *err_tab;
    err_tab = fopen("err_of_omega.dat","w");

    int i,N = 129;
    double omega = 1.;

    Grid initial_cross, force3;

    Grid_InitByFunction(&force3,&f2,-1, 1, -1, 1, N);
    Grid_Init (&initial_cross,      -1, 1, -1, 1, N, 0.001);

    Grid_InitDirihlet_w_d(&initial_cross,   &e2);

    Cross_IterationSet_w_f(&initial_cross,1.,&force3, 15000);
//        Grid_print_error(&test_ido3,&e3);
    IDO_InitNormalD_2(&initial_cross,&force3);


    for(omega; omega<2.;omega+=0.015)
    {
        Grid test_ido3;

        Grid_Copy(&test_ido3,&initial_cross);

        for(i=0;i<40;i++)
        {
            IDO_Ori_IterationSet_w_f_Seidel(&test_ido3,omega,&force3,1);
            fprintf(err_tab,"%d %f %f\n",i,omega,log10(Grid_print_error(&test_ido3,&e2)/4.));
//            Grid_print_error(&test_ido3,&e3);
        }
        fprintf(err_tab,"\n");
    }

    return 0;
}

/*
 * =======================================
 *
 * Test for grid intensifier... TBD
 *
 */
int main()
{
    Grid CG1, CG2, CG3, CG4;
    Grid IG1, IG2, IG3, IG4;
    Grid F1, F2, F3, F4;

    Grid_InitByFunction(&F1, &f2, -1,1, -1,1, 9);
    Grid_InitByFunction(&F2, &f2, -1,1, -1,1, 17);
    Grid_InitByFunction(&F3, &f2, -1,1, -1,1, 33);
    Grid_InitByFunction(&F4, &f2, -1,1, -1,1, 65);


    printf("cross iterations section:\n");

    Grid_Init(&CG1, -1,1, -1,1, 9, 0.);
    Cross_IterationSet_w_f(&CG1,1.,&F1,500);
    printf("9x9 on 500 Iters:\t%f\n",log10(Grid_print_error(&CG1,&e2)/4.));

    printf("\n");

    Grid_Intensify(&CG2, &CG1);
    printf("17x17 on interp: \t%f\n",log10(Grid_print_error(&CG2,&e2)/4.));
    Cross_IterationSet_w_f(&CG2,1.,&F2,500);
    printf("17x17 on 500 Iters\t%f\n",log10(Grid_print_error(&CG2,&e2)/4.));

    printf("\n");

    Grid_Intensify(&CG3, &CG2);
    printf("33x33 on interp: \t%f\n",log10(Grid_print_error(&CG3,&e2)/4.));
    Cross_IterationSet_w_f(&CG3,1.,&F3,500);
    printf("33x33 on 500 Iters\t%f\n",log10(Grid_print_error(&CG3,&e2)/4.));

    printf("\n");

    Grid_Intensify(&CG4, &CG3);
    printf("65x65 on interp: \t%f\n",log10(Grid_print_error(&CG4,&e2)/4.));
    Cross_IterationSet_w_f(&CG4,1.,&F4,500);
    printf("65x65 on 500 Iters\t%f\n",log10(Grid_print_error(&CG4,&e2)/4.));

    printf("\n");

    printf("IDO iterations section:\n");

    Grid_Copy(&IG1, &CG1);
    Grid_Copy(&IG2, &CG2);
    Grid_Copy(&IG3, &CG3);
    Grid_Copy(&IG4, &CG4);

    // IDO iterations
    IDO_InitNormalD_2(&IG1,&F1);
    IDO_Ori_IterationSet_w_f_Seidel(&IG1,1.8,&F1,1);
    printf("9x9 on 1   Iters:\t%f\n",log10(Grid_print_error(&IG1,&e2)/4.));

    IDO_InitNormalD_2(&IG2,&F2);
    IDO_Ori_IterationSet_w_f_Seidel(&IG2,1.8,&F2,2);
    printf("17x17 on 1 Iters:\t%f\n",log10(Grid_print_error(&IG2,&e2)/4.));

    IDO_InitNormalD_2(&IG3,&F3);
    IDO_Ori_IterationSet_w_f_Seidel(&IG3,1.8,&F3,5);
    printf("33x33 on 1 Iters:\t%f\n",log10(Grid_print_error(&IG3,&e2)/4.));

    IDO_InitNormalD_2(&IG4,&F4);
    IDO_Ori_IterationSet_w_f_Seidel(&IG4,1.8,&F4,21);
    printf("65x65 on 1 Iters:\t%f\n",log10(Grid_print_error(&IG4,&e2)/4.));

    //add auto iteration identifier... optimize the process.
    //improve accuracy of intepolation in grid intensification
    return 0;
}
