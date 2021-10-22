comp:=gfortran
pattern:=*.f90
paral:=-fopenmp
opt:=-O2
source:=$(wildcard $(pattern))
   obj:=$(patsubst %.f90, %.o, $(source))
prog: $(obj)
	$(comp) $(opt) $^ -o $@
%.o %.mod : %.f90
	$(comp) -c $(opt) $<
main.o : init_data.mod SIEI.mod coord_grid_nodes.mod
gauss_quad.o : init_data.mod
coord_grid_nodes.o : init_data.mod
SIEI.o : init_data.mod coord_grid_nodes.mod integrator.mod SSE.mod
SSE.o : init_data.mod
clear:
	rm -f *.o *.mod prog
result: prog
	./prog
	rm -f *.o *.mod prog
