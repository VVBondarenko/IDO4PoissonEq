//
// Created by vbondarenko on 19.01.17.
//

#ifndef IDO4POISSONEQ_CROSS_H
#define IDO4POISSONEQ_CROSS_H

#include <grid.h>

void Cross_Iteration        (Grid *Task, double omega);
void Cross_IterationSet    (Grid *Task, double omega, int N);
void Cross_Iteration_w_force(Grid *Task, double omega, Grid *force);
void Cross_IterationSet_w_f(Grid *Task, double omega, Grid *force, int N);
void Cross_IterationSet_w_f_w_autostop(Grid *Task, double omega, Grid *force, int N,
                                            double prev_diff, double *return_diff);


#endif //IDO4POISSONEQ_CROSS_H
