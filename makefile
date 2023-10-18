#begin makefile
# use this to debug: FLAGS=  -O0 -fdefault-real-8 -fdefault-double-8 -fcheck=all
OBJS=   charPDE_new.o                                                                
FLAGS=  -O3 -fdefault-real-8 -fdefault-double-8  
FLAGS=  -O0 -fdefault-real-8 -fdefault-double-8 -fcheck=all     
##FLAGS=  -O2 -freal-8-real-16  # this is for quadrpul execute

# The compiler                        
F90=    gfortran

charPDE_new:        charPDE_new.o
	$(F90) -o charPDE_new $(FLAGS) $(OBJS) need.a

charPDE_new.o:       charPDE_new.f90
	$(F90) -c $(FLAGS) charPDE_new.f90

clean:
	rm -f charPDE_new $(OBJS) *.mod *__*




