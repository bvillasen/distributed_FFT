#include "io.h"
#include<stdlib.h>
#include <unistd.h>
#include <mpi.h>

/* MPI-safe printf routine */
int print_single(const char * __restrict sdata, ...){
  int code = 0;
  va_list ap;
  va_start(ap, sdata);
  code = vfprintf(stdout, sdata, ap);
  va_end(ap);

  return code;
}

void print_mpi( const string &text_out, int rank, int nprocs  ){
  
  for (int i=0; i<nprocs; i++ ){
    if ( rank == i ) printf( "%s\n", text_out.c_str() );
    MPI_Barrier(MPI_COMM_WORLD);    
    usleep( 100 );
  }
  
  
  
  
}

void Load_field_from_file( const string &field_name, Real *data_field, int nx_local, int ny_local, int nz_local,  const string &file_name, const string &input_dir, int rank, int nprocs  ){
  
  
  // if ( rank == 0 ) printf("Loading Field: %s\n", field_name.c_str() );
  // 
  // string out_text = " Loading File: " + input_dir + file_name;
  // print_mpi( out_text, rank, nprocs );
  
  hid_t  file_id;
  herr_t  status;
  // open the file
  file_id = H5Fopen((input_dir + file_name).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    printf("Unable to open input file: %s\n", (input_dir + file_name).c_str() );
    exit(0);
  }
  
  int i, j, k, id, buf_id;
  hid_t     attribute_id, dataset_id; 
  Real      *dataset_buffer;
  
  // need a dataset buffer to remap fastest index
  dataset_buffer = (Real *) malloc(nx_local*ny_local*nz_local*sizeof(Real));
  
  
  string field_to_load = "/" + field_name;
  dataset_id = H5Dopen(file_id, field_to_load.c_str(), H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  status = H5Dclose(dataset_id);
  
  for (k=0; k<nz_local; k++) {
    for (j=0; j<ny_local; j++) {
      for (i=0; i<nx_local; i++) {
        id = i + j*nx_local + k*nx_local*ny_local;
        buf_id = k + j*nz_local + i*nz_local*ny_local;
        data_field[id] = dataset_buffer[buf_id];
      }
    }
  }

  free(dataset_buffer);
  // close the file
  status = H5Fclose(file_id);
    
}

void Write_field_to_file( const string &field_name, Real *data_field, int nx_local, int ny_local, int nz_local, hid_t file_id  ){
  
  // printf(" Writing Field: %s\n", field_name.c_str() );
  int       i, j, k, id, buf_id;
  hid_t     dataset_id, dataspace_id; 
  Real      *dataset_buffer;
  herr_t    status;
  hsize_t   dims[3];
  
  // Create the data space for the datasets
  dims[0] = nx_local;
  dims[1] = ny_local;
  dims[2] = nz_local;
  dataspace_id = H5Screate_simple(3, dims, NULL);
  
  
  dataset_buffer  = (Real *) malloc(nx_local*ny_local*nz_local*sizeof(Real));
  for (k=0; k<nz_local; k++) {
    for (j=0; j<ny_local; j++) {
      for (i=0; i<nx_local; i++) {
        id = i + j*nx_local + k*nx_local*ny_local;
        buf_id = k + j*nz_local + i*nz_local*ny_local;
        dataset_buffer[buf_id] = data_field[id];
      }
    }
  }
  
  // Create a dataset id for density
  dataset_id = H5Dcreate(file_id, field_name.c_str(), H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // Write the density array to file  // NOTE: NEED TO FIX FOR FLOAT REAL!!!
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  // Free the dataset id
  status = H5Dclose(dataset_id);
  
  // Free the dataspace id
  status = H5Sclose(dataspace_id);
  
  
  
  
  
  
  
  
  
  
  free(dataset_buffer);
  
}
