FCdebug = ifort -g -check all -check bounds -traceback -fpp
IFORT = ifort -O3  -g -traceback -fpp
GFORT = gfortran -O3 -g -fbacktrace -D ISGFORTRAN -cpp
FC = $(IFORT)

all: read_ppm

ifort: FC = $(IFORT)
ifort: clean all

gfort: FC = $(GFORT)
gfort: clean all

read_ppm: read_ppm.f read_ppm.h format_adump_header2.o
	$(FC) -o read_ppm read_ppm.f format_adump_header2.o
	cp read_ppm ~/bin

format_adump_header2.o: format_adump_header2.f
	$(FC) -c format_adump_header2.f

clean:
	rm -f read_ppm *.o
