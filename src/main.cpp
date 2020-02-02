#include <mpi.h>
#include <stdio.h>
#include <pfft.h>
#include <complex.h>
#include <string>
#include <iostream> 
#include "io.h"
#include <sstream>



using namespace std;

int main(int argc, char** argv) {
  
  
  // Initialize MPI 
  MPI_Init(&argc, &argv);
  // Get the number of processes
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // Get the rank of the process
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  // Print off a hello world message
  // printf("Hello world from processor %s, rank %d out of %d processors\n", processor_name, rank, size);
  
  
  string data_dir, input_dir, output_dir;
  data_dir = "/home/bruno/Desktop/ssd_0/data/";
  input_dir  = data_dir + "cosmo_sims/256_dm_50Mpc/output_files/";
  output_dir = data_dir + "cosmo_sims/256_dm_50Mpc/output_files/data_fft/";
  
  if ( rank == 0 ){
    printf("Distributed FFT \n" );
    printf("InputDir: %s\n",  input_dir.c_str() );
    printf("OutputDir: %s\n", output_dir.c_str() );
  }
  
  
  int n_snapshot = 0;
  ostringstream file_name;
  file_name <<  "_particles.h5." << rank;
  
  
  Load_field_from_file( file_name.str(), input_dir, rank, size   );
  
  
  
  
  // int nprocs[3];
  // ptrdiff_t n_total[3];
  // ptrdiff_t alloc_local;
  // ptrdiff_t local_n_input[3], local_input_start[3];
  // ptrdiff_t local_n_output[3], local_output_start[3];
  // double err, *in;
  // pfft_complex *out;
  // pfft_plan plan_forw=NULL, plan_back=NULL;
  // MPI_Comm comm_cart_3d;
  // 
  // //Set size of processes grid
  // nprocs[0] = 2; nprocs[1] = 2; nprocs[2] = 2;
  // 
  // //Set size of FFT grid 
  // n_total[0] = 32; n_total[1] = 32; n_total[2] = 32;
  // 
  // 
  // 
  // 
  // // Initialize PFFT          
  // pfft_init();
  // 
  // 
  // /* Create three-dimensional process grid of size nprocs[0] x nprocs[1] x nprocs[2], if possible */
  // if( pfft_create_procmesh(3, MPI_COMM_WORLD, nprocs, &comm_cart_3d) ){
  //  pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", nprocs[0]*nprocs[1]*nprocs[2]);
  //  MPI_Finalize();
  //  return 1;
  // }
  // 
  // /* Get parameters of data distribution */
  // alloc_local = pfft_local_size_dft_r2c_3d(n_total, comm_cart_3d, PFFT_TRANSPOSED_NONE,
  //    local_n_input, local_input_start, local_n_output, local_output_start);
  // 
  // /* Allocate memory */
  // in  = pfft_alloc_real(2 * alloc_local);
  // out = pfft_alloc_complex(alloc_local);
  // 
  // /* Plan parallel forward FFT */
  // plan_forw = pfft_plan_dft_r2c_3d(
  //    n_total, in, out, comm_cart_3d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
  // 
  // /* Plan parallel backward FFT */
  // plan_back = pfft_plan_dft_c2r_3d(
  //    n_total, out, in, comm_cart_3d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
  // 
  // /* Initialize input with random numbers */
  // pfft_init_input_real(3, n_total, local_n_input, local_input_start,
  //    in);
  // 
  // /* execute parallel forward FFT */
  // pfft_execute(plan_forw);
  // 
  // /* clear the old input */
  // pfft_clear_input_real(3, n_total, local_n_input, local_input_start,
  //    in);
  // 
  // /* execute parallel backward FFT */
  // pfft_execute(plan_back);
  // 
  // /* Scale data */
  // for(ptrdiff_t l=0; l < local_n_input[0] * local_n_input[1] * local_n_input[2]; l++)
  //  in[l] /= (n_total[0]*n_total[1]*n_total[2]);
  // 
  // /* Print error of back transformed data */
  // MPI_Barrier(MPI_COMM_WORLD);
  // err = pfft_check_output_real(3, n_total, local_n_input, local_input_start, in, comm_cart_3d);
  // pfft_printf(comm_cart_3d, "Error after one forward and backward trafo of size n=(%td, %td, %td):\n", n_total[0], n_total[1], n_total[2]); 
  // pfft_printf(comm_cart_3d, "maxerror = %6.2e;\n", err);
  // 
  // /* free mem and finalize */
  // pfft_destroy_plan(plan_forw);
  // pfft_destroy_plan(plan_back);
  // MPI_Comm_free(&comm_cart_3d);
  // pfft_free(in); pfft_free(out);
  // 
  
  
  MPI_Finalize();
  return 0;



  
}