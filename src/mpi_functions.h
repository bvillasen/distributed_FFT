#ifndef  MPI_ROUTINES_H
#define  MPI_ROUTINES_H

#include <stddef.h>
#include "global.h"


/* MPI reduction wrapper for max(Real)*/
Real ReduceRealMax(Real x);

/* MPI reduction wrapper for min(Real)*/
Real ReduceRealMin(Real x);

/* MPI reduction wrapper for sum(Real)*/
Real ReduceRealSum(Real x);

/* MPI reduction wrapper for avg(Real)*/
Real ReduceRealAvg(Real x, int nproc);








#endif /*MPI_ROUTINES_H*/