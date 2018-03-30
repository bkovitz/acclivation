DOT = dot -Tpdf

all: sa

sa: sa.c Makefile
	gcc sa.c --std=c99 -Werror -Wall -D_POSIX_C_SOURCE=199309L -g -o sa -lm

run: sa
	./sa

out: sa
	./sa > out
	#grep "epoch fitness" out
	tail -3 out

outs: sa
	./sa > out
	tail -3 out
	./sa >> out
	tail -3 out
	./sa >> out
	tail -3 out
	./sa >> out
	tail -3 out
	./sa >> out
	tail -3 out

%.pdf: %.dot
	$(DOT) < $< > $@

test.dot: sa
	./sa > test.dot

test.pdf: test.dot
	dot -Tpdf < test.dot > test.pdf

dot: test.pdf
	evince test.pdf

phplot:
	sed -n '/BEGIN PHFUNC/,/END PHFUNC/ {//!p}' out > phfunc
	./plot_xyz.py phfunc 0 1 2

vfplot:
	sed -n '/BEGIN VFUNC/,/END VFUNC/ {//!p}' out > vfunc
	./plot_xyz.py vfunc 0 1 4

phrangeplot:
	sed -n '/BEGIN VFUNC/,/END VFUNC/ {//!p}' out > vfunc
	./plot_xyz.py vfunc 2 3 4 scatter

plots:
	make vfplot &
	make phrangeplot &
	make phplot &

fitness: sa
	./sa > fitness.out
	./plot_fitness.py fitness.out

OUT=out
with_seed:
	./sa `grep seed $(OUT) | cut -d'=' -f2`

NRUNS = $(shell seq 1 20)
many: sa
	@./sa > /tmp/out
	@tail -3 /tmp/out
	@$(foreach i,$(NRUNS),./sa >> /tmp/out; tail -3 /tmp/out;)

clean:
	rm sa test.pdf test.dot

tags:
	ctags *.[ch]

.PHONY: tags run all dot clean plot fitness with_seed out
