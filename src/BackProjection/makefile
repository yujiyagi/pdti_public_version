F = ifort
Opt1 = -O3 -ipo -mavx -parallel -par-threshold50 -qopenmp -traceback -g
#Opt1 = -O3 -ipo -xAVX -parallel -par-threshold50 -qopenmp -traceback -g
all: ../../bin/RunBP 

../../bin/RunBP :  RunBP.f90 sub.RunBP.f90
	$(F) -o $@ $^ $(Opt1)

