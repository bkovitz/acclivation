all: sa

sa: sa.c
	gcc sa.c --std=c99 -D_POSIX_C_SOURCE=199309L -g -o sa -lm

run: sa
	./sa

outs: sa
	./sa > out
	./sa >> out
	./sa >> out
	./sa >> out
	./sa >> out
	./sa >> out

test.dot: sa
	./sa > test.dot

test.pdf: test.dot
	dot -Tpdf < test.dot > test.pdf

dot: test.pdf
	evince test.pdf

clean:
	rm test.pdf test.dot

tags:
	ctags *.[ch]

.PHONY: tags run all dot clean
