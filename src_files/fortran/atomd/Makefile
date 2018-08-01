include ./../../../Makefile.in

all : gaunt.so gutils.so

gaunt.so: gaunt.f90
	$(CMP) $(F2PL) -c $? -m gaunt 


gutils.so: gutils.cxx
	python setup.py build_ext --inplace


clean:
	rm -rf build gutils.so gaunt.so

