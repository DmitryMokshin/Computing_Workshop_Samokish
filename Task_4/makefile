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
main.o : icad.mod finding_eigenvalues.mod fun_for_matrix.mod
sse.o : icad.mod
fun_for_matrix.o : icad.mod
finding_eigenvalues.o : icad.mod sse.mod fun_for_matrix.mod
clear:
	rm -f *.o *.mod prog
results: prog
	./prog
resultall: prog
	./prog all
