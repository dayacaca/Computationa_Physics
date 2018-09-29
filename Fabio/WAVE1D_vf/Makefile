#FC = mpif90
FC=gfortran
#FC=/opt/intel/bin/ifort

#FFLAGS = -O3 -fopenmp
FFLAGS = -O3

#LNK = mpif90
LNK=gfortran
#LNK=/opt/intel/bin/ifort

OBJS = evolve.o allocate.o initializing_arrays.o grid.o initial.o rhs.o main.o check_parameters.o save1Ddata_x.o RK2.o RK3.o RK4.o Info_Screen.o BC.o exact.o

MODS = arrays.o global_numbers.o

$(OBJS):	$(MODS)

cafe:  $(OBJS) $(MODS)
		$(LNK) $(FFLAGS) -o xwave $(OBJS) $(MODS) 
	@ mkdir -p xxx
	@ mv xwave xxx

.PHONY:	clean

clean:
	-rm -f *.o *.mod xxx/xwave xxx/*xl

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@
