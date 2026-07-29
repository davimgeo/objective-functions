FLAGS = -std=c99 -Wall -fopenmp -mavx2 -mfma -O3 -march=native
INCLUDE = -Iinclude

LIBS = -lfftw3f -lm

#MAIN = main1d.c
MAIN = test.c

SRCS = src/1D/*.c

run:
	gcc $(FLAGS) $(MAIN) $(INCLUDE) $(SRCS) -o run.out $(LIBS)
	./run.out
	$(MAKE) clean

clean:
	rm -f run.out

plot:
	$(MAKE) run

	python3 plots/plot_1d.py
