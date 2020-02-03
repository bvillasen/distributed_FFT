#include <mpi.h>
#include <math.h>
#include "mpi_functions.h"




/* MPI reduction wrapper for max(Real)*/
Real ReduceRealMax(Real x)
{
  Real in = x;
  Real out;
  Real y;

  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  y = (Real) out;
  return y;
}


/* MPI reduction wrapper for min(Real)*/
Real ReduceRealMin(Real x)
{
  Real in = x;
  Real out;
  Real y;
  
  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  y = (Real) out;
  return y;
}

/* MPI reduction wrapper for sum(Real)*/
Real ReduceRealSum(Real x )
{
  Real in = x;
  Real out;
  Real y;

  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  y = (Real) out ;
  return y;
}

/* MPI reduction wrapper for avg(Real)*/
Real ReduceRealAvg(Real x, int nproc )
{
  Real in = x;
  Real out;
  Real y;

  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  y = (Real) out / nproc;
  return y;
}