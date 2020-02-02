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

# HDF5 Location
HDF5_INCL = -I/usr/include/hdf5/serial/
HDF5_LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5

INCL   = -I./ $(HDF5_INCL)
LIBS   = -lm $(HDF5_LIBS) 

# PFFT in Tornado
FFTW_INCL = -I/home/bruno/apps/fftw-3.3.5/include
FFTW_LIBS = -L/home/bruno/apps/fftw-3.3.5/lib -lfftw3
PFFT_INCL = -I/home/bruno/apps/pfft-git/include
PFFT_LIBS = -L/home/bruno/apps/pfft-git/lib  -lpfft  -lfftw3_mpi -lfftw3
INCL += $(FFTW_INCL) $(PFFT_INCL)
LIBS += $(FFTW_LIBS) $(PFFT_LIBS)











%.o:	%.c
		$(CC) $(CFLAGS)  $(INCL)  -c $< -o $@ 

%.o:	%.cpp
		$(CXX) $(CXXFLAGS)  $(INCL) -c $< -o $@ 

$(EXEC): $(OBJS) 
	 	 $(CXX) $(OBJS)  $(LIBS) -o $(EXEC)



.PHONY : clean

clean:
	 rm -f $(OBJS)  $(EXEC)