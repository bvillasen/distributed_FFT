#include <mpi.h>
#include <stdio.h>
#include <pfft.h>
#include <complex.h>
#include <string>
#include <iostream> 
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "io.h"
#include "mpi_functions.h"




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
  // data_dir = "/home/bruno/Desktop/ssd_0/data/";
  // input_dir  = data_dir + "cosmo_sims/256_dm_50Mpc/output_files/";
  // output_dir = data_dir + "cosmo_sims/256_dm_50Mpc/output_files/data_fft/";
  data_dir = "/gpfs/alpine/proj-shared/ast149/";
  input_dir  = data_dir + "cosmo_sims/2048_hydro_50Mpc/output_files_hm12/";
  output_dir = data_dir + "cosmo_sims/2048_hydro_50Mpc/snapshots_hm12/power_spectrum/data_fft/";
  
  if ( rank == 0 ){
    printf("Distributed FFT \n" );
    printf("InputDir: %s\n",  input_dir.c_str() );
    printf("OutputDir: %s\n", output_dir.c_str() );
  }
  
  
  int n_proc_x, n_proc_y, n_proc_z;
  int nx_total, ny_total, nz_total;
  int nx_local, ny_local, nz_local;
  int n_cells_local, n_cells_total;
  Real Lx, Ly, Lz;
  
  Lx = 50000.0;
  Ly = 50000.0;
  Lz = 50000.0;
  
  n_proc_x = 8;
  n_proc_y = 8;
  n_proc_z = 8;
  
  nx_total = 2048;
  ny_total = 2048;
  nz_total = 2048;
  
  nx_local = nx_total / n_proc_x;
  ny_local = ny_total / n_proc_y;
  nz_local = nz_total / n_proc_z;
  
  n_cells_local = nx_local*ny_local*nz_local;
  n_cells_total = nx_total*ny_total*nz_total;
  
  if ( rank == 0 ){
    printf("N Procs: [ %d %d %d]\n", n_proc_x, n_proc_y, n_proc_z );
    printf("Grid Global: [ %d %d %d]\n", nx_total, ny_total, nz_total );
    printf("Grid Global: [ %d %d %d]\n", nx_local, ny_local, nz_local );
  }
  
  //Allocate space for data
  Real *data_field = (Real *) malloc(nx_local*ny_local*nz_local*sizeof(Real)); 
  
  int n_snapshot = 0;
  ostringstream in_file_name;
  in_file_name << n_snapshot << "_particles.h5." << rank;
  string field_name = "density";
  Load_field_from_file( field_name, data_field, nx_local, ny_local, nz_local, in_file_name.str(), input_dir, rank, size   );
  
  Real field_mean_local, field_mean_global;
  field_mean_local = 0;
  for ( int i=0; i<n_cells_local; i++ ) field_mean_local += data_field[i];
  field_mean_local /= n_cells_local;
  field_mean_global = ReduceRealAvg( field_mean_local, size );
  if ( rank == 0 ) printf("Mean %s: %f\n", field_name.c_str(), field_mean_global );
  
  // Compute the overdensity
  for ( int i=0; i<n_cells_local; i++ ) data_field[i] = data_field[i] / field_mean_global;
  
  
  
  int nprocs[3];
  ptrdiff_t n_total[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_n_input[3], local_input_start[3];
  ptrdiff_t local_n_output[3], local_output_start[3];
  double error;
  pfft_complex *in, *out;
  pfft_plan plan_forw=NULL, plan_back=NULL;
  MPI_Comm comm_cart_3d;
  
  //Set size of processes grid
  nprocs[0] = n_proc_z; nprocs[1] = n_proc_y; nprocs[2] = n_proc_x;
  
  //Set size of FFT grid 
  n_total[0] = nz_total; n_total[1] = ny_total; n_total[2] = nx_total;
    
  // Initialize PFFT          
  pfft_init();
  
  
  /* Create three-dimensional process grid of size nprocs[0] x nprocs[1] x nprocs[2], if possible */
  if( pfft_create_procmesh(3, MPI_COMM_WORLD, nprocs, &comm_cart_3d) ){
   pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: This test file only works with %d processes.\n", nprocs[0]*nprocs[1]*nprocs[2]);
   MPI_Finalize();
   return 1;
  }
  
  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_3d(n_total, comm_cart_3d, PFFT_TRANSPOSED_NONE,
     local_n_input, local_input_start, local_n_output, local_output_start);
  
  //Check the dimensions of the input and the transform, make sure that dimensions and offsets  are the same
  bool domain_error = false;
  for (int i=0; i<3; i++ ){
    if ( local_n_input[i] != local_n_output[i]) domain_error = true;
    if ( local_input_start[i] != local_output_start[i]) domain_error = true;
    // if ( rank == 0 ) printf("n_in: %d   n_out: %d   \n", local_n_input[i], local_n_output[i]  );
  }
  //Exit if the dimensions and offsets are not the same
  if ( domain_error ){
   printf("PFFT: FFT doamin error \n" );
   exit(-1);
  }
  
  /* Allocate memory */
  in  = pfft_alloc_complex(alloc_local);
  out = pfft_alloc_complex(alloc_local);
  
  bool compute_backward = false;
  
  /* Plan parallel forward FFT */
  if ( rank == 0 ) printf("Creating FFT Plan Forward \n" );
  plan_forw = pfft_plan_dft_3d(
     n_total, in, out, comm_cart_3d, PFFT_FORWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
  
  if (compute_backward){
    /* Plan parallel backward FFT */
    if ( rank == 0 ) printf("Creating FFT Plan Backward \n" );
    plan_back = pfft_plan_dft_3d(
       n_total, out, in, comm_cart_3d, PFFT_BACKWARD, PFFT_TRANSPOSED_NONE| PFFT_MEASURE| PFFT_DESTROY_INPUT);
  }
  
  // Set the input as the data field
  for ( int i=0; i<n_cells_local; i++ ){
    in[i][0] = data_field[i];
    in[i][1] = 0;
  }
  
  /* execute parallel forward FFT */
  if ( rank == 0 ) printf("Excecuting FFT Forward \n" );
  pfft_execute(plan_forw);
  
  // Compute the amplitude of the FFT
  Real *fft_amp2 = (Real *) malloc(nx_local*ny_local*nz_local*sizeof(Real)); 
  for ( int i=0; i<n_cells_local; i++ ) fft_amp2[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
  
  
  //Compute Wave Numbers kx, ky, kz
  Real *k_mag;
  k_mag = (Real *) malloc(nx_local*ny_local*nz_local*sizeof(Real)); 
  ptrdiff_t k_0, k_1, k_2;
  Real k_x, k_y, k_z, k_magnitude, k_mag_min, k_mag_max;
  int id;
  k_mag_min = 1e100;
  k_mag_max = -1e100;
  for (int k=0; k<nz_local; k++ ){
    for (int j=0; j<ny_local; j++ ){
      for (int i=0; i<nx_local; i++ ){
        id = i + j*nx_local + k*nx_local*ny_local;
        
        k_0 = k + local_input_start[0];
        k_1 = j + local_input_start[1];
        k_2 = i + local_input_start[2];
  
        if ( k_2 >= nx_total/2) k_2 -= nx_total;
        if ( k_0 >= nz_total/2) k_0 -= nz_total;
        if ( k_1 >= ny_total/2) k_1 -= ny_total;
        
        k_x = 2*M_PI*k_2/Lx;
        k_z = 2*M_PI*k_0/Ly;
        k_y = 2*M_PI*k_1/Lz;
        k_magnitude = sqrt( k_x*k_x + k_y*k_y + k_z*k_z);
        k_mag[id] = k_magnitude;
        k_mag_max = fmax( k_mag_max, k_magnitude );
        if (k_magnitude > 0 ) k_mag_min = fmin( k_mag_min, k_magnitude );
      }
    }
  }
  
  if ( compute_backward ){
    /* clear the old input */
    pfft_clear_input_complex_3d(n_total, local_n_input, local_input_start, in);
    
    /* execute parallel backward FFT */
    if ( rank == 0 ) printf("Excecuting FFT Backward \n" );
    pfft_execute(plan_back);
    
    /* Scale data */
    for(ptrdiff_t l=0; l < local_n_input[0] * local_n_input[1] * local_n_input[2]; l++){
      in[l][0] /= (n_total[0]*n_total[1]*n_total[2]);
      in[l][1] /= (n_total[0]*n_total[1]*n_total[2]);
    }
    
    // // /* Print error of back transformed data */
    // Real error_local, error_global;
    // error_local = 0;
    // for ( int i=0; i<n_cells_local; i++ ) error_local += in[i][0] - data_field[i];
    // error_global = ReduceRealSum( error_local );
    // if ( rank == 0 ) printf("Error Global: %f\n", error_global );
  }
  
  //Save to hdf5 file  
  ostringstream out_file_name;
  out_file_name  << n_snapshot << "_data_fft" << ".h5." << rank;
  hid_t   file_id; /* file identifier */
  herr_t  status;

  string out_text = " Writing File: " + output_dir + out_file_name.str();
  print_mpi( out_text, rank, size  );
  // Create a new file using default properties.
  file_id = H5Fcreate( (output_dir + out_file_name.str()).c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
  
  // Write Header  
  hsize_t   attr_dims;
  hid_t     attribute_id, dataspace_id;
  int       int_data[3];
  
  
  // 3D attributes
  attr_dims = 3;
  // Create the data space for the attribute
  dataspace_id = H5Screate_simple(1, &attr_dims, NULL);
  
  int_data[0] = nx_total; int_data[1] = ny_total; int_data[2] = nz_total;
  attribute_id = H5Acreate(file_id, "dims", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Awrite(attribute_id, H5T_NATIVE_INT, int_data);
  status = H5Aclose(attribute_id);
  
  int_data[0] = nx_local; int_data[1] = ny_local; int_data[2] = nz_local;
  attribute_id = H5Acreate(file_id, "dims_local", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Awrite(attribute_id, H5T_NATIVE_INT, int_data);
  status = H5Aclose(attribute_id);
  
  int_data[0] = local_input_start[2]; int_data[1] = local_input_start[1]; int_data[2] = local_input_start[0];
  attribute_id = H5Acreate(file_id, "offset", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Awrite(attribute_id, H5T_NATIVE_INT, int_data);
  status = H5Aclose(attribute_id);
  
  
  field_name = "input";
  Write_field_to_file( field_name, data_field, nx_local, ny_local, nz_local, file_id  );
  
  field_name = "fft_amp2";
  Write_field_to_file( field_name, fft_amp2, nx_local, ny_local, nz_local, file_id  );  
  
  field_name = "k_mag";
  Write_field_to_file( field_name, k_mag, nx_local, ny_local, nz_local, file_id  );
  
  // Single attributes first
  attr_dims = 1;
  dataspace_id = H5Screate_simple(1, &attr_dims, NULL);
  attribute_id = H5Acreate(file_id, "k_mag_min", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &k_mag_min);
  status = H5Aclose(attribute_id);
  
  attribute_id = H5Acreate(file_id, "k_mag_max", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &k_mag_max);
  status = H5Aclose(attribute_id);
  
  // Close the dataspace
  status = H5Sclose(dataspace_id);

  
  // close the file
  status = H5Fclose(file_id);
  if (status < 0) {printf("File write failed.\n"); exit(-1); }

  
  
  
  
  
  
  
  /* free mem and finalize */
  pfft_destroy_plan(plan_forw);
  if ( compute_backward ) pfft_destroy_plan(plan_back);
  MPI_Comm_free(&comm_cart_3d);
  pfft_free(in); pfft_free(out);  
  free( data_field );
  free( k_mag );
  
  MPI_Barrier(MPI_COMM_WORLD);
  if ( rank == 0 ) printf("Finished Successfully\n" );
  
  
  MPI_Finalize();
  return 0;



  
}