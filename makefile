# This is an commentary line in a makefile
# Start of the makefile

ifortErr = $(shell which ifort >/dev/null; echo $$?)

ifeq "$(ifortErr)" "0"
FC=ifort
FLAGS= -free -qopenmp -no-wrap-margin
else
FC=gfortran
FLAGS= -free -fopenmp
endif

primospectra: mod1.o main.o mod2.o mod3.o mod3_5.o mod4.o mod5.o spectral.o solve.o
	 $(FC) $(FLAGS)   *.o -o primospectra  

mod1.mod: mod1.o mod1.f90
	 $(FC) $(FLAGS)   -c mod1.f90  
mod1.o: mod1.f90
	 $(FC) $(FLAGS)   -c mod1.f90  

mod2.mod: mod2.o mod2.f90
	 $(FC) $(FLAGS)   -c mod2.f90  
mod2.o: mod2.f90
	 $(FC) $(FLAGS)   -c mod2.f90  

mod3.mod: mod3.o mod3.f90
	 $(FC) $(FLAGS)   -c mod3.f90  
mod3.o: mod3.f90
	 $(FC) $(FLAGS)   -c mod3.f90  

mod3_5.mod: mod3_5.o mod3_5.f90
	 $(FC) $(FLAGS)   -c mod3_5.f90  
mod3_5.o: mod3_5.f90
	 $(FC) $(FLAGS)   -c mod3_5.f90  
	
mod4.mod: mod4.o mod4.f90
	 $(FC) $(FLAGS)   -c mod4.f90  
mod4.o: mod4.f90
	 $(FC) $(FLAGS)   -c mod4.f90  


mod5.mod: mod5.o mod5.f90
	 $(FC) $(FLAGS)   -c mod5.f90  
mod5.o: mod5.f90
	 $(FC) $(FLAGS)   -c mod5.f90  


spectral.mod: spectral.o spectral.f90
	 $(FC) $(FLAGS)   -c spectral.f90  
spectral.o: spectral.f90
	 $(FC) $(FLAGS)   -c spectral.f90  


solve.mod: solve.o solve.f90
	 $(FC) $(FLAGS)   -c solve.f90  
solve.o: mod5.f90
	 $(FC) $(FLAGS)   -c solve.f90  


main.o : mod1.mod mod2.mod mod3.mod mod3_5.mod mod4.mod mod5.mod spectral.mod solve.mod main.f90
	 $(FC) $(FLAGS)   -c main.f90  

cleanall:
	rm -r core *.o primospectra *.mod *.pdf *.png *.dat *.log *.txt *.TXT .fuse* ../data ../plot

clean:
	rm -r core *.o primospectra *.mod *.pdf *.png *.dat *.log *.txt *.TXT .fuse* 


cleanfuse:
	rm -rf ../data/*.fuse*

cleandata:
	 rm -r ../data/ #~/.local/share/Trash
	#sudo ra -fr ~/.local/share/Trash/*

cleanplot:
	rmdir plot

