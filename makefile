FLAGS = -std=c++11 -Wall -fopenmp -mavx2 -mfma -O3 -march=native

LIBS = -lm

MAIN = main.cpp

SRCS = src/dft2d.cpp src/objf.cpp

run:
	g++ $(FLAGS) $(MAIN) $(SRCS) -o run.out $(LIBS)
	./run.out
	$(MAKE) clean

clean:
	rm -f run.out

plot:
	$(MAKE) run

	python3 plot.py
