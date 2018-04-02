LAT = latex -shell-escape
#TEX = pdflatex --shell-escape
TEXINPUTS=.:./sty:
TEX = TEXINPUTS=.:./sty: latexmk -pdf -xelatex
BIB = bibtex8
TEXFILES = $(wildcard *.tex)
PDFFILES = $(TEXFILES:.tex=.pdf)
BIBFILES = $(wildcard *.bib)
DOT = dot -Tpdf

OK_LINE_ARGS = --ridge_type=0 --bumps=1 --ridge_radius=0.2 --activation_types=1 --mutation_type_ub=10 \
  --knob_type=0 --multi_edges=0 --peak_movement=0 --output_types=0 --c2=1 --c3=0 \
  --spreading_rate=0.2 --edge_weights=0 --c1_lb=0.1 --c1_ub=0.9 --extra_mutation_rate=0.05 \
  --decay=0.8 --allow_move_edge=0

all: sa

$(PDFFILES): $(BIBFILES)
%.pdf: %.tex
	$(TEX) $<

%.pdf: %.dot
	$(DOT) < $< > $@

sa: sa.c Makefile
	gcc sa.c --std=c99 -Werror -Wall -D_POSIX_C_SOURCE=199309L -g -o sa -lm

data: sa run.py add_param_set.py
	./run.py > d.csv
	./add_param_set.py d.csv > data.csv

run: sa
	./sa

out: sa
	./sa > out
	#grep "epoch fitness" out
	tail -3 out

okline: sa
	./sa $(OK_LINE_ARGS) > out
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

outs2: sa
	./sa --ridge_type=1 --bumps=0 --ridge_radius=0.05 --activation_types=1 \
--mutation_type_ub=10 --knob_type=0 --multi_edges=1 --peak_movement=1 \
--output_types=1 --c2=2 --c3=-.45 --spreading_rate=0.1 --edge_weights=1 \
--c1_lb=-1 --c1_ub=1 --extra_mutation_rate=0.01 --decay=0.8 --allow_move_edge=0 \
> out
	./sa --ridge_type=1 --bumps=0 --ridge_radius=0.05 --activation_types=1 \
--mutation_type_ub=10 --knob_type=0 --multi_edges=1 --peak_movement=1 \
--output_types=1 --c2=2 --c3=-.45 --spreading_rate=0.1 --edge_weights=1 \
--c1_lb=-1 --c1_ub=1 --extra_mutation_rate=0.01 --decay=0.8 --allow_move_edge=0 \
>> out
	./sa --ridge_type=1 --bumps=0 --ridge_radius=0.05 --activation_types=1 \
--mutation_type_ub=10 --knob_type=0 --multi_edges=1 --peak_movement=1 \
--output_types=1 --c2=2 --c3=-.45 --spreading_rate=0.1 --edge_weights=1 \
--c1_lb=-1 --c1_ub=1 --extra_mutation_rate=0.01 --decay=0.8 --allow_move_edge=0 \
>> out
	./sa --ridge_type=1 --bumps=0 --ridge_radius=0.05 --activation_types=1 \
--mutation_type_ub=10 --knob_type=0 --multi_edges=1 --peak_movement=1 \
--output_types=1 --c2=2 --c3=-.45 --spreading_rate=0.1 --edge_weights=1 \
--c1_lb=-1 --c1_ub=1 --extra_mutation_rate=0.01 --decay=0.8 --allow_move_edge=0 \
>> out
	./sa --ridge_type=1 --bumps=0 --ridge_radius=0.05 --activation_types=1 \
--mutation_type_ub=10 --knob_type=0 --multi_edges=1 --peak_movement=1 \
--output_types=1 --c2=2 --c3=-.45 --spreading_rate=0.1 --edge_weights=1 \
--c1_lb=-1 --c1_ub=1 --extra_mutation_rate=0.01 --decay=0.8 --allow_move_edge=0 \
>> out

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

tom: sa
	@./sa --seed=520664716 \
			--num_epochs=200 \
  		--ridge_radius=0.200000 \
  		--c2=1.000000 \
			--c3=0.000000 \
			--decay=0.900000 \
  		--spreading_rate=0.200000 \
  		--distance_weight=10.000000 \
  		--bumps=1 \
  		--mutation_type_ub=16 \
  		--extra_mutation_rate=0.100000 \
  		--crossover_freq=0.300000 \
  		--edge_inheritance=5 \
  		--edge_weights=0 \
  		--activation_types=1 \
  		--num_candidates=7 \
  		--generations_per_epoch=20 \
  		--sa_timesteps=20

clean:
	rm sa test.pdf test.dot

tags:
	ctags *.[ch]

.PHONY: tags run all dot clean plot fitness with_seed out tom
