# This is an commentary line in a makefile
# Start of the makefile
primospectra: mod1.o main.o mod2.o mod3.o mod3_5.o mod4.o mod5.o spectral.o solve.o
	 gfortran -free   *.o -o primospectra  -fopenmp

mod1.mod: mod1.o mod1.f90
	 gfortran -free   -c mod1.f90  -fopenmp
mod1.o: mod1.f90
	 gfortran -free   -c mod1.f90  -fopenmp

mod2.mod: mod2.o mod2.f90
	 gfortran -free   -c mod2.f90  -fopenmp
mod2.o: mod2.f90
	 gfortran -free   -c mod2.f90  -fopenmp

mod3.mod: mod3.o mod3.f90
	 gfortran -free   -c mod3.f90  -fopenmp
mod3.o: mod3.f90
	 gfortran -free   -c mod3.f90  -fopenmp

mod3_5.mod: mod3_5.o mod3_5.f90
	 gfortran -free   -c mod3_5.f90  -fopenmp
mod3_5.o: mod3_5.f90
	 gfortran -free   -c mod3_5.f90  -fopenmp
	
mod4.mod: mod4.o mod4.f90
	 gfortran -free   -c mod4.f90  -fopenmp
mod4.o: mod4.f90
	 gfortran -free   -c mod4.f90  -fopenmp


mod5.mod: mod5.o mod5.f90
	 gfortran -free   -c mod5.f90  -fopenmp
mod5.o: mod5.f90
	 gfortran -free   -c mod5.f90  -fopenmp


spectral.mod: spectral.o spectral.f90
	 gfortran -free   -c spectral.f90  -fopenmp
spectral.o: spectral.f90
	 gfortran -free   -c spectral.f90  -fopenmp


solve.mod: solve.o solve.f90
	 gfortran -free   -c solve.f90  -fopenmp
solve.o: mod5.f90
	 gfortran -free   -c solve.f90  -fopenmp


main.o : mod1.mod mod2.mod mod3.mod mod3_5.mod mod4.mod mod5.mod spectral.mod solve.mod main.f90
	 gfortran -free   -c main.f90  -fopenmp

cleanall:
	rm -r *.o primospectra *.mod *.pdf *.png *.dat *.log *.txt *.TXT .fuse* ../data ../plot pbss/igfortran* pbss/*.txt

clean:
	rm -r *.o primospectra *.mod *.pdf *.png *.dat *.log *.txt *.TXT .fuse* 


cleanfuse:
	rm -rf ../data/*.fuse*

cleandata:
	 rm -r ../data/ #~/.local/share/Trash
	#sudo ra -fr ~/.local/share/Trash/*

cleanplot:
	rmdir plot

