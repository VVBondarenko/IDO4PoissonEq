//
// Created by vbondarenko on 19.01.17.
//

#ifndef IDO4POISSONEQ_IDO_ORIGINAL_H
#define IDO4POISSONEQ_IDO_ORIGINAL_H

#include <grid.h>

void IDO_InitDeriv  (Grid *Task);
void IDO_InitNormalD(Grid *Task);
void IDO_InitNormalD_2(Grid *Task, Grid *F);

void IDO_Iteration  (Grid *Task);
void IDO_IterationOriginal     (Grid *Task, double omega);
void IDO_IterationOriginal_w_f (Grid *Task, double omega, Grid *force);
void IDO_IterationModified     (Grid *Task, double omega);
void IDO_IterationModified_w_f (Grid *Task, double omega, Grid *force);
void IDO_IterationSet          (Grid *Task, double omega, int N);
void IDO_Mod_IterationSet_w_f  (Grid *Task, double omega, Grid *force, int N);
void IDO_Ori_IterationSet_w_f  (Grid *Task, double omega, Grid *force, int N);
void IDO_Ori_IterationSet_w_f_Jacoby(Grid *Task, double omega, Grid *force, int N);
void IDO_Ori_IterationSet_w_f_Jacoby_adapt(Grid *Task, double omega, Grid *force, int N, double (*f_exact)(double, double));
void IDO_Ori_IterationSet_w_f_Seidel(Grid *Task, double omega, Grid *force, int N);
void IDO_Ori_IterationSet_w_f_Seidel_ILU(Grid *Task, double omega, Grid *force, int N);

#endif //IDO4POISSONEQ_IDO_ORIGINAL_H
