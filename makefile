FLAGS = -std=c++17 -Wall -lm

run:
	g++ $(FLAGS) main.cpp -o main.out
	./main.out

clean:
	rm -f main.out

plot:
	python3 plot.py
