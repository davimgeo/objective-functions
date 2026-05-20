FLAGS = -std=c++17 -Wall

LIBS = -lfftw3f -lm

SRCS = src/dft2d.cpp

run:
	g++ $(FLAGS) $(SRCS) -o run.out $(LIBS)
	./run.out

clean:
	rm -f run.out

plot:
	python3 plot.py
