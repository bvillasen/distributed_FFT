#include "io.h"
#include <unistd.h>
#include <mpi.h>
#include<hdf5.h>

void print_mpi( int rank, int nprocs, const string &text_out  ){
  
  for (int i=0; i<nprocs; i++ ){
    if ( rank == i ) printf( "%s\n", text_out.c_str() );
    MPI_Barrier(MPI_COMM_WORLD);    
    usleep( 100 );
  }
  
  
  
  
}

void Load_field_from_file( const string &file_name, const string &input_dir, int rank, int nprocs  ){
  
  string out_text = "Loading File: " + input_dir + file_name;
  
  print_mpi( rank, nprocs, out_text );
  
  
  
  
  
  
  
}