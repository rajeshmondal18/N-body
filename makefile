LINKLIB=  -L/usr/local/lib/  -fopenmp -lfftw3f_omp -lfftw3f -lm

#LINKLIB= -lm -lc -L/usr/local/lib/ -lfftw3

INCLUDE=-I/usr/local/include/


CFLAGS=-g
CC=gcc


nbody_comp: nbody_comp.o nbody_funcs.o allotarrays.o funcs.o powerspec.o tf_fit.o 
	$(CC) $(CFLAGS) -o nbody_comp nbody_comp.o nbody_funcs.o funcs.o allotarrays.o powerspec.o tf_fit.o $(LINKLIB)
	rm -rf *.o
	rm -rf *~

nbody_comp.o:	nbody_comp.c
	$(CC) -c $(CFLAGS) $(INCLUDE) nbody_comp.c $(LINKLIB)


ionz_main: ionz_main.o nbody_funcs.o allotarrays.o funcs.o powerspec.o tf_fit.o ionz_funcs.o  
	$(CC) $(CFLAGS) -o ionz_main  ionz_main.o nbody_funcs.o allotarrays.o funcs.o powerspec.o tf_fit.o ionz_funcs.o  $(LINKLIB)
	rm -rf *.o
	rm -rf *~

ionz_main.o:	ionz_main.c
	$(CC) -c $(CFLAGS) $(INCLUDE) ionz_main.c $(LINKLIB)

write_sample: write_sample.o nbody_funcs.o allotarrays.o funcs.o powerspec.o tf_fit.o  
	$(CC) $(CFLAGS) -o write_sample  write_sample.o nbody_funcs.o allotarrays.o funcs.o powerspec.o tf_fit.o  $(LINKLIB)
	rm -rf *.o
	rm -rf *~

write_sample.o:	write_sample.c
	$(CC) -c $(CFLAGS) $(INCLUDE) write_sample.c $(LINKLIB)


nbody_funcs.o:	nbody_funcs.c
	$(CC) -c  $(CFLAGS)   $(INCLUDE) nbody_funcs.c $(LINKLIB)

ionz_funcs.o:	ionz_funcs.c
	$(CC) -c $(CFLAGS) $(INCLUDE) ionz_funcs.c $(LINKLIB) 

allotarrays.o:	allotarrays.c
	$(CC) -c $(CFLAGS) $(INCLUDE) allotarrays.c $(LINKLIB)

powerspec.o:	powerspec.c
	$(CC) -c $(CFLAGS) $(INCLUDE) powerspec.c $(LINKLIB)

funcs.o:	funcs.c
	$(CC) -c  $(CFLAGS)   $(INCLUDE) funcs.c $(LINKLIB)

tf_fit.o:	tf_fit.c
	$(CC) -c $(CFLAGS) $(INCLUDE) tf_fit.c $(LINKLIB)


fof_main:	fof_main.o nbody_funcs.o allotarrays.o funcs.o powerspec.o tf_fit.o 
	$(CC) $(CFLAGS) -o fof_main fof_main.o nbody_funcs.o funcs.o allotarrays.o powerspec.o tf_fit.o $(LINKLIB)
	rm -rf *.o
	rm -rf *~

fof_main.o:	fof_main.c 
	$(CC) -c $(CFLAGS) $(INCLUDE) fof_main.c $(LINKLIB)



mass_func:	mass_func.o nbody_funcs.o allotarrays.o tf_fit.o powerspec.o funcs.o
	$(CC) $(CFLAGS) -o mass_func mass_func.o nbody_funcs.o allotarrays.o powerspec.o tf_fit.o funcs.o $(LINKLIB)
	rm -rf *.o
	rm -rf *~

mass_func.o:	mass_func.c 
	$(CC) -c $(CFLAGS) $(INCLUDE) mass_func.c $(LINKLIB)

funcs_fit.o:	funcs_fit.c 
	$(CC) -c $(CFLAGS) $(INCLUDE) funcs_fit.c $(LINKLIB)


clean:
	rm -rf *.o
	rm -rf *~
