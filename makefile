FLAGS = -std=c++11 -Wall -fopenmp -mavx2 -mfma -O3 -march=native

LIBS = -lfftw3f -lm

SRCS = src/dft2d.cpp

run:
	g++ $(FLAGS) $(SRCS) -o run.out $(LIBS)
	./run.out

clean:
	rm -f run.out

plot:
	$(MAKE) run

	python3 plot.py
