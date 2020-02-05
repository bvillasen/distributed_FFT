#ifndef IO_H
#define IO_H

#include <iostream>
#include <string>
#include <stdio.h>
#include "global.h"
#include<hdf5.h>

using namespace std;

void Load_field_from_file( const string &field_name, Real *data_field, int nx_local, int ny_local, int nz_local, const string &file_name, const string &input_dir, int rank, int nprocs );

void Write_field_to_file( const string &field_name, Real *data_field, int nx_local, int ny_local, int nz_local, hid_t file_id  );


void print_mpi(const string &text_out, int rank, int nprocs  );

int print_single(const char * __restrict sdata, ...)

#endif //IO_H