prog: main.o premier.o lineaire.o gen.o
	gcc -o prog main.o premier.o lineaire.o gen.o -lgmp

main.o: main.c
	gcc -c main.c -lgmp

premier.o: premier.c
	gcc -c premier.c -lgmp

lineaire.o: lineaire.c
	gcc -c lineaire.c -lgmp

gen.o: gen.c
	gcc -c gen.c -lgmp

clean :
	rm -f prog *.o


