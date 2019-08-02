idr: main.o GenEM.o GenEM_mixture.o controller.o BF_CEF.o
	g++ main.o GenEM.o GenEM_mixture.o controller.o BF_CEF.o -o idr -lgsl -lgslcblas
main.o: main.cc
	g++ -c main.cc
controller.o: controller.cc controller.h
	g++ -c controller.cc
BF_CEF.o: BF_CEF.h BF_CEF.cc
	g++ -c BF_CEF.cc
GenEM_mixture.o: GenEM_mixture.h GenEM_mixture.cc
	g++ -c GenEM_mixture.cc
GenEM.o: GenEM.h GenEM.cc
	g++ -c GenEM.cc
clean: 
	rm *.o idr
