#-----------------------
# SGI-Calhoun
#MYFC = mpif77
#MYCC = cc
#FLAGS = -O3 -xT -ip -no-prec-div -ansi-alias -fno-alias -mcmodel=large
#
#-----------------------
# Ls
#MYFC = mpif90
#MYCC = gcc
#FLAGS = -O3 -openmp -tpp7 -xW -fpe0 -w -ip -Ob2 -pc32 -ftz -fpp1 -vec_report2 -I/usr/local/include
#
#-----------------------
IFORT = ifort
GFORT = gfortran
MYFC = $(IFORT)
MYCC = gcc
#FLAGS = -fast
#FLAGS = -O3 -openmp -tpp7 -xW -fpe0 -w -ip -Ob2 -pc32 -ftz -fpp1 -vec_report2 -I/usr/local/include
#FLAGS = -O1 -openmp

##########################
#######FLAGS = -g -openmp -check all -traceback -I/usr/local/include


# >>>>>>>>>>>>>>>

IFORT_FLAGS = -O3 -g -traceback -openmp -tpp7 -xW -fpe0 -w -ip -Ob2 -pc32 -ftz -fpp1 -I/usr/local/include -fpp
GFORT_FLAGS = -O3 -g -fbacktrace -fopenmp -ffpe-trap=zero,overflow,underflow -w -I/usr/local/include -D ISGFORTRAN -cpp
FLAGS = $(IFORT_FLAGS)
#FLAGS = -g -openmp -traceback -check all

FLAGSavx = -g -g -traceback -openmp -O3 -ip -axCORE_AVX2 -Zp16 -align array32byte -mcmodel=large -unroll-aggressive -fno-alias -prefetch -inline-level=2

# <<<<<<<<<<<<<<<<


DEFS=-DLINUX -O3
#HDF_ROOT=/soft/hdf5/hdf5-1.8.7-intel-serial
HDF_ROOT=/home/tech/porter/src/a3d/hdf5128
ZLIBDIR=/home/tech/porter/src/a3d/hdf5128
HDF5LIBS = $(HDF_ROOT)/lib/libhdf5_fortran.a $(HDF_ROOT)/lib/libhdf5.a $(ZLIBDIR)/lib/libz.a

E3D_BIN_DIR=bin

 # switched_io.o mnt.o BlockTree.o avminmax.o fft.o
OBJ    = e3d.o read_adump_e3d.o switched_io.o mnt.o BlockTree.o avminmax.o sofplot.o vecplot.o fft3d.o principle.o write_tec.o
OBJPDT = e3d.o read_pdt_e3d.o   switched_io.o mnt.o BlockTree.o avminmax.o sofplot.o vecplot.o fft3d.o principle.o write_tec.o
OBJHDF = e3d.o read_hdf5_e3d.o read_hdf5.o switched_io.o mnt.o BlockTree.o avminmax.o sofplot.o vecplot.o fft3d.o principle.o write_tec.o

default: e3d vert2stats

ifort: MYFC = $(IFORT)
ifort: FLAGS = $(IFORT_FLAGS)
ifort: clean default

gfort: MYFC = $(GFORT)
gfort: FLAGS = $(GFORT_FLAGS)
gfort: clean default

gen_plt_file: gen_plt_file.c
	icc -o gen_plt_file gen_plt_file.c

install: e3d vert2stats
	mkdir -p $(E3D_BIN_DIR)
	cp e3d $(E3D_BIN_DIR)
	cp vert2stats $(E3D_BIN_DIR)
	cp ee* $(E3D_BIN_DIR)

#install: e3d e3d_pdt e3d_hdf vert2stats gen_plt_file
#	cp e3d $(E3D_BIN_DIR)
#	cp e3d_pdt $(E3D_BIN_DIR)
#	cp e3d_hdf $(E3D_BIN_DIR)
#	cp vert2stats $(E3D_BIN_DIR)
#	cp gen_plt_file $(E3D_BIN_DIR)

e3d: $(OBJ)
	$(MYFC) $(FLAGS) -o e3d  $(OBJ)

e3d_pdt: $(OBJPDT)
	$(MYFC) $(FLAGS) -o e3d_pdt  $(OBJPDT)

e3d_hdf: $(OBJHDF)
	$(MYFC) $(FLAGS) -o e3d_hdf  $(OBJHDF) $(HDF5LIBS)

# -lsvml -static-libcxa

e3d.o: e3d.f read_adump_e3d.h postproc_info_e3d.h amr_e3d.h
	$(MYFC) $(FLAGS) -c e3d.f

read_adump_e3d.o: read_adump_e3d.f read_adump_e3d.h amr_e3d.h
	$(MYFC) $(FLAGS) -c read_adump_e3d.f

read_pdt_e3d.o: read_pdt_e3d.f read_pdt_e3d.h
	$(MYFC) $(FLAGS) -c read_pdt_e3d.f

read_hdf5_e3d.o: read_hdf5_e3d.f read_hdf5_e3d.h
	$(MYFC) $(FLAGS) -c read_hdf5_e3d.f

read_hdf5.o: read_hdf5.f90
	$(MYFC) $(FLAGS) -c read_hdf5.f90 -I$(HDF_ROOT)/include


sofplot.o: sofplot.f
	 $(MYFC) $(FLAGS) -c sofplot.f

vecplot.o: vecplot.f
	 $(MYFC) $(FLAGS) -c vecplot.f

principle.o: principle.f
	 $(MYFC) $(FLAGS) -c principle.f

fft3d.o: fft3d.f
	 $(MYFC) $(FLAGS) -c fft3d.f

fft.o: fft.f
	$(MYFC) $(FLAGS) -c fft.f

write_tec.o: write_tec.f
	$(MYFC) $(FLAGS) -c write_tec.f

BlockTree.o: BlockTree.c
	$(MYCC) $(DEFS) -c BlockTree.c

avminmax.o: avminmax.f
	$(MYFC) $(FLAGS) -c avminmax.f


switched_io.o: switched_io.c
	$(MYCC) -c switched_io.c

mnt.o: mnt.c
	$(MYCC) -c mnt.c

parse: parse.o
	$(MYFC) $(FLAGS) parse.o -o parse  -lsvml -static-libcxa

parse.o: parse.f
	$(MYFC) $(FLAGS) -c parse.f

vert2stats: vert2stats.f
	ifort -g -o vert2stats vert2stats.f

clean:
	rm -f *.o e3d e3d_pdt
