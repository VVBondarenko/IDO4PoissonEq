//
// Created by vbondarenko on 19.01.17.
//

#ifndef IDO4POISSONEQ_IDO_ORIGINAL_H
#define IDO4POISSONEQ_IDO_ORIGINAL_H

#include <grid.h>

void Grid_IDO_InitDeriv (Grid *Task);
void Grid_IDO_Iteration (Grid *Task);
void Grid_IDO_IterationOriginal     (Grid *Task, double omega);
void Grid_IDO_IterationOriginal_w_f (Grid *Task, double omega, Grid *force);
void Grid_IDO_IterationModified     (Grid *Task, double omega);
void Grid_IDO_IterationModified_w_f (Grid *Task, double omega, Grid *force);
void Grid_IDO_IterationSet          (Grid *Task, double omega, int N);
void Grid_IDO_Mod_IterationSet_w_f  (Grid *Task, double omega, Grid *force, int N);
void Grid_IDO_Ori_IterationSet_w_f  (Grid *Task, double omega, Grid *force, int N);

#endif //IDO4POISSONEQ_IDO_ORIGINAL_H
