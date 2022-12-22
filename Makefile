all: main.cpp
	g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` main.cpp -o main`python3-config --extension-suffix` -lblas -I/usr/include/mkl -I/usr/include/python3.8

test:all
	python3 -m pytest

clean:
	rm -rf *.so __pycache__