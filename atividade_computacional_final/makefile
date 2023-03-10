CC=g++
CFLAGS=-O3 -lm    
DEPS = aproximacao_splines.h romberg.h 
OBJ = main.o aproximacao_splines.o romberg.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

exe: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean: 
	rm -rf *.o
