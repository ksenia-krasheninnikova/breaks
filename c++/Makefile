rootPath=/Users/admin/projects/
halLibPath=${rootPath}/hal/lib/
sonLibPath=${rootPath}/sonLib/lib/
libSourcesAll = $(wildcard impl/*.cpp)

cpp = h5c++ ${h5prefix}

psl2synteny : main.cpp ${halLibPath}/halLib.a
	${cpp} -std=c++11 -gdwarf-3 -O0 -I${halLibPath} -I${sonLibPath} ${halLibPath}/halLib.a main.cpp -o psl2synteny

#TODO?
#${halLibPath}/halLib.a: 

psl_merger.o: psl_merger.cpp
	g++ -std=c++11 -gdwarf-3 -O0 psl_merger.cpp -c -o psl_merger.o
	
clean:
	rm -rf *.o
	rm -rf psl2synteny*
