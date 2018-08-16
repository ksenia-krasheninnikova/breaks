rootPath=/Users/admin/projects/
halLibPath=${rootPath}/hal/lib/
sonLibPath=${rootPath}/sonLib/lib/
#libSourcesAll = $(wildcard impl/*.cpp)

cpp = h5c++ ${h5prefix}

hal2synteny : main.cpp hal2psl.h psl.h  psl_io.cpp  psl_io.h
	${cpp} -std=c++11 -O0 -g -I${sonLibPath} -I${halLibPath} main.cpp ${halLibPath}/halLiftover.a ${halLibPath}/halLib.a  ${sonLibPath}/sonLib.a -o hal2synteny

#TODO?
#${halLibPath}/halLib.a: 

psl_merger.o: psl_merger.cpp
	g++ -std=c++11 -gdwarf-3 -O0 psl_merger.cpp -c -o psl_merger.o
	
clean:
	rm -rf *.o
	rm -rf hal2synteny*
