objects = Graph.o Value.o Estimation.o Vertex.o
sources= Graph Value Estimation Vertex
cppfiles= Graph.cpp Value.cpp Estimation.cpp Vertex.cpp 
headers= Graph.hpp Value.hpp Estimation.hpp CommonFile.hpp Vertex.hpp
main= main.cpp

all: parity_game final

parity_game: $(cppfiles) $(headers) 
	for i in $(sources) ; do \
	(g++ -std=c++11 -Wall -c $$i.cpp -o $$i.o); \
	done

final: $(objects) $(main)
	(g++ -std=c++11 -o mpg -Wall $(objects) main.cpp ../lp_dev/liblpsolve55.a -Wl,--no-as-needed -ldl); \
	done
	
clean:
	rm -rf $(objects) mpg
	rm -rf mmpggenerator
	rm -rf tokenringgenerator

mmpggenerator: mmpg.cpp
	g++ mmpg.cpp -o mmpggenerator

tokenringgenerator: tokenring.cpp
	g++ tokenring.cpp -o tokenringgenerator

