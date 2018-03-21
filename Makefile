all:
	gcc sa.c --std=c99 -g -o sa -lm

run:
	./sa
