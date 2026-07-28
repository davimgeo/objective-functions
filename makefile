FLAGS = -std=c++11 -Wall -fopenmp -mavx2 -mfma -O3 -march=native
INCLUDE = -Iinclude

LIBS = -lfftw3 -lm

#MAIN = main1d.cpp
MAIN = test.cpp

SRCS = src/1D/*.cpp

run:
	g++ $(FLAGS) $(MAIN) $(INCLUDE) $(SRCS) -o run.out $(LIBS)
	./run.out
	$(MAKE) clean

clean:
	rm -f run.out

plot:
	$(MAKE) run

	python3 plots/plot_1d.py
