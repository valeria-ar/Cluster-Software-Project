FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main.o matrixB.o listGroup.o listVertex.o spmat.o
	gcc main.o matrixB.o listGroup.o listVertex.o spmat.o -o cluster $(LIBS)
	
clean:
	rm -rf *.o cluster

main.o: main.c	
	gcc $(FLAGS) -c main.c	

matrixB.o: matrixB.c
	gcc $(FLAGS) -c matrixB.c
	
listGroup.o: listGroup.c
	gcc $(FLAGS) -c listGroup.c
	
listVertex.o: listVertex.c
	gcc $(FLAGS) -c listVertex.c
	
spmat.o: spmat.c
	gcc $(FLAGS) -c spmat.c

