all: sa

sa: sa.c
	gcc sa.c --std=c99 -Werror -Wall -D_POSIX_C_SOURCE=199309L -g -o sa -lm

run: sa
	./sa

out: sa
	./sa > out

outs: sa
	./sa > out
	./sa >> out
	./sa >> out
	./sa >> out
	./sa >> out
	grep average out

test.dot: sa
	./sa > test.dot

test.pdf: test.dot
	dot -Tpdf < test.dot > test.pdf

dot: test.pdf
	evince test.pdf

phplot:
	./plot_xyz.py phfunc 0 1 2

fitness: sa
	./sa > fitness.out
	./plot_fitness.py fitness.out

OUT=out
with_seed:
	./sa `grep seed $(OUT) | cut -d'=' -f2`

many: sa
	./sa > out
	for i in $$(seq 1 50); do ./sa >> out; done
	grep average out

clean:
	rm test.pdf test.dot

tags:
	ctags *.[ch]

.PHONY: tags run all dot clean plot fitness with_seed out
