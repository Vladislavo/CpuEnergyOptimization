CC = gcc

sources = main.c enrgopt_funcs.c
dependencies = enrgopt.h

enrgopt: $(sources) $(dependencies)
	$(CC) $(sources) -lm -o $@

clean:
	rm enrgopt