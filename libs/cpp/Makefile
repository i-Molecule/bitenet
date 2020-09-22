all: libdescriptor.so libmarchingcubes.so

libdescriptor.so: extern.cpp
	g++ -std=c++17 -fPIC -O3 -shared -I ./include/ ./source/* extern.cpp -o ../libdescriptor.so

libmarchingcubes.so: marching_cubes.cpp
	g++ -std=c++17 -fPIC -O3 -shared -I ./include/ marching_cubes.cpp -o ../libmarchingcubes.so

clean:
	rm -f ../*.so
