CC = gcc

CFLAGS = -std=gnu99 -O3 -g -march=native -fopenmp -mavx2 -mfma

INCLUDE = \
	-I../../include \
	-I../../src \
	-I../../libs/wave-modeling/include \
	-I../../libs/wave-modeling/src \
	-Iconfig

LIBS = -lfftw3f -lm

OBJF_SRC = $(shell find ../../src -name "*.c")

WAVE_SRC = $(shell find ../../libs/wave-modeling/src -name "*.c")

GRADIENT_SRC = main.c config/config.c

SRC = $(OBJF_SRC) $(WAVE_SRC) $(GRADIENT_SRC)

run:
	@$(CC) $(CFLAGS) $(INCLUDE) $(SRC) $(LIBS) -o gradient_models.out
	@./gradient_models.out

clean:
	rm -f *.out
