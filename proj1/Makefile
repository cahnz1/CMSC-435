output: main.o
	g++ main.o -o output

main.o: main.cpp 
	g++ -c main.cpp -O3 -I/usr/include/eigen3/ 

clean:
	rm *.o output
