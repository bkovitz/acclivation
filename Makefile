DOT = dot -Tpdf

all: sa

sa: sa.c Makefile
	gcc sa.c --std=c99 -Werror -Wall -D_POSIX_C_SOURCE=199309L -g -o sa -lm

run: sa
	./sa

out: sa
	./sa > out
	grep "epoch fitness" out

outs: sa
	./sa > out
	./sa >> out
	./sa >> out
	./sa >> out
	./sa >> out
	grep "epoch fitness" out

%.pdf: %.dot
	$(DOT) < $< > $@

test.dot: sa
	./sa > test.dot

test.pdf: test.dot
	dot -Tpdf < test.dot > test.pdf

dot: test.pdf
	evince test.pdf

phplot:
	./plot_xyz.py phfunc 0 1 2

vfplot:
	./plot_xyz.py vfunc 0 1 4

phrangeplot:
	./plot_xyz.py vfunc 2 3 4 scatter

fitness: sa
	./sa > fitness.out
	./plot_fitness.py fitness.out

OUT=out
with_seed:
	./sa `grep seed $(OUT) | cut -d'=' -f2`

NRUNS = $(shell seq 1 20)
many: sa
	@./sa > out
	@tail -1 out
	@$(foreach i,$(NRUNS),./sa >> out; tail -1 out;)

clean:
	rm sa test.pdf test.dot

tags:
	ctags *.[ch]

.PHONY: tags run all dot clean plot fitness with_seed out
