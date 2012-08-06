#Makefile for vasp_structure_randomise program

#If you are using mac then use this line
#CFLAGS=-O2 -DMAC

#On systems that have normal getline() use this line
CFLAGS=-O2

#Library dependencies you need kdtree from code.google.com/p/kdtree
LDFLAGS=-lm -lkdtree 

all: vasp_structure_randomise.c
	$(CC) vasp_structure_randomise.c -o vasp_structure_randomise $(CFLAGS) $(LDFLAGS)
clean:
	rm vasp_structure_randomise

