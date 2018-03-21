all: sa

sa: sa.c
	gcc sa.c --std=c99 -g -o sa -lm

run: sa
	./sa

test.dot: sa
	./sa > test.dot

test.pdf: test.dot
	dot -Tpdf < test.dot > test.pdf
