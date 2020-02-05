EXEC   = distributed_FFT

OPTIMIZE =  -O2  

DIR = ./src
CFILES = $(wildcard $(DIR)/*.c)
CPPFILES = $(wildcard $(DIR)/*.cpp)

OBJS   = $(subst .c,.o,$(CFILES)) $(subst .cpp,.o,$(CPPFILES))  
COBJS   = $(subst .c,.o,$(CFILES)) 
CPPOBJS   = $(subst .cpp,.o,$(CPPFILES)) 


CC	= mpicc
CXX   = mpicxx


#PRECISION = -DPRECISION=1
PRECISION = -DPRECISION=2

# # HDF5 Location in Tornado
# HDF5_INCL = -I/usr/include/hdf5/serial/
# HDF5_LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5

# # HDF5 Location in Summit
# HDF5_INCL = -I$(OLCF_HDF5_ROOT)/include
# HDF5_LIBS = -L$(OLCF_HDF5_ROOT)/lib -lhdf5

# HDF5 Location in Lux
HDF5_INCL = -I/cm/shared/apps/hdf5_18/1.8.20/include
HDF5_LIBS = -L/cm/shared/apps/hdf5_18/1.8.20/lib -lhdf5


INCL   = -I./ $(HDF5_INCL)
LIBS   = -lm $(HDF5_LIBS) 

# # PFFT in Tornado
# FFTW_INCL = -I/home/bruno/apps/fftw-3.3.5/include
# FFTW_LIBS = -L/home/bruno/apps/fftw-3.3.5/lib -lfftw3
# PFFT_INCL = -I/home/bruno/apps/pfft-git/include
# PFFT_LIBS = -L/home/bruno/apps/pfft-git/lib  -lpfft  -lfftw3_mpi -lfftw3

# # PFFT in Summit
# FFTW_INCL = -I/ccs/proj/ast149/code/fftw/include
# FFTW_LIBS = -L/ccs/proj/ast149/code/fftw/lib -lfftw3
# PFFT_INCL = -I/ccs/proj/ast149/code/pfft/include
# PFFT_LIBS = -L/ccs/proj/ast149/code/pfft/lib  -lpfft  -lfftw3_mpi -lfftw3


# # PFFT in Lux
FFTW_INCL = -I/home/brvillas/code/fftw-3.3.8/include
FFTW_LIBS = -L/home/brvillas/code/fftw-3.3.8/lib -lfftw3
PFFT_INCL = -I/home/brvillas/code/pfft/include
PFFT_LIBS = -L/home/brvillas/code/pfft/lib  -lpfft  -lfftw3_mpi -lfftw3


INCL += $(FFTW_INCL) $(PFFT_INCL)
LIBS += $(FFTW_LIBS) $(PFFT_LIBS)






FLAGS = $(PRECISION)
CFLAGS = $(FLAGS)
CXXFLAGS = $(FLAGS)



%.o:	%.c
		$(CC) $(CFLAGS)  $(INCL)  -c $< -o $@ 

%.o:	%.cpp
		$(CXX) $(CXXFLAGS)  $(INCL) -c $< -o $@ 

$(EXEC): $(OBJS) 
	 	 $(CXX) $(OBJS)  $(LIBS) -o $(EXEC)



.PHONY : clean

clean:
	 rm -f $(OBJS)  $(EXEC)