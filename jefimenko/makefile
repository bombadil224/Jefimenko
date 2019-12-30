# latest: constants.so charge_current_integrals.so  simulation.so
latest: charge_current_integrals.so
	touch -am latest

# constants.so: constants.o
#	gcc -shared -o constants.so constants.o `python-config --libs`

# constants.o: constants.c
#	cc `pkg-config --libs --cflags python3` -c constants.c -o constants.o -fPIC

# constants.c: constants.pyx
#	cython3 constants.pyx -a



charge_current_integrals.so: charge_current_integrals.o
	gcc -shared -o charge_current_integrals.so charge_current_integrals.o `python-config --libs`


charge_current_integrals.o: charge_current_integrals.c
	cc `pkg-config --libs --cflags python3` -c charge_current_integrals.c -o charge_current_integrals.o -fPIC

charge_current_integrals.c: charge_current_integrals.pyx
	cython3 charge_current_integrals.pyx -a



# simulation.so: simulation.o
#	gcc -shared -o simulation.so simulation.o `python-config --libs`

# simulation.o: simulation.c
#	cc `pkg-config --libs --cflags python3` -c simulation.c -o simulation.o -fPIC

# simulation.c: simulation.pyx
#	cython3 simulation.pyx -a


clean :
#	rm constants.c
#	rm constants.o
	rm charge_current_integrals.c
	rm charge_current_integrals.o
#	rm simulation.c
#	rm simulation.o
