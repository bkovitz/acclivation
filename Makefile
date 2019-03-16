CC=gcc
CFLAGS=--std=gnu99 -Werror -Wall -g
LDFLAGS=--std=gnu99 -g
LDLIBS=-lm

LAT = latex -shell-escape
#TEX = pdflatex --shell-escape
TEXINPUTS=.:./sty:
TEX = TEXINPUTS=.:./sty: latexmk -pdf -xelatex
BIB = bibtex8
TEXFILES = $(wildcard *.tex)
PDFFILES = $(TEXFILES:.tex=.pdf)
BIBFILES = $(wildcard *.bib)
DOT = dot -Tpdf
ifeq ($(shell uname -s), Darwin)
VIEW_PDF = open
else
VIEW_PDF = evince
endif

#Thin ridge along y=x
YXLINE = --ridge_type=0 --c1_lb=-1 --c1_ub=1 

OBLIQUE_LINE = $(YXLINE) --ridge_type=0 --c2=2.0 --c3=-0.45 --c1_lb=-0.275 --c1_ub=0.725

#--c2=1 --c3=0 \
#--c2=2.0 --c3=-0.45 \

CIRCLE = --ridge_type=1 --bumps=0 --ridge_radius=0.1 --peak_movement=1 \
	--c1_lb=-1 --c1_ub=1

RESTRICTED =
#RESTRICTED = --c1_lb=0.2 --c1_ub=0.9

#A parameter set for experimentation. Try the good ideas here, run with
#'make x', and save noteworthy parameter sets under a different name.
X_ARGS = $(YXLINE) $(RESTRICTED) --bumps=0 \
	--reward_coverage=0 \
	--num_epochs=80 --generations_per_epoch=20 \
	--num_organisms=80 --num_candidates=6 \
	--num_nodes=4 --num_edges=1 \
	--input_accs=1 --activation_types=3 --output_types=0 --knob_type=0 \
	--mutation_type_ub=10 --extra_mutation_rate=0.00 --crossover_freq=0.1 \
	--multi_edges=0 --allow_move_edge=0 --edge_weights=1 --edge_inheritance=1 \
	--log=ancestors \
	--spreading_rate=0.05 --decay=1.0 --log=ancestors --seed=968766798 --invu=1 #--seed=1043614093

#QUESTION: What happens to the winner here when you turn a knob? What happens
#when you make single graph edits? Why is this graph stuck?
LOOK_AT_THIS = $(YXLINE) --bumps=1 \
	--num_epochs=40 --generations_per_epoch=20 \
	--num_organisms=80 --num_candidates=6 \
	--num_nodes=4 --num_edges=4 \
	--input_accs=1 --activation_types=3 --output_types=0 --knob_type=0 \
	--mutation_type_ub=10 --extra_mutation_rate=0.00 --crossover_freq=0.05 \
	--multi_edges=0 --allow_move_edge=0 --edge_weights=0 --edge_inheritance=5 \
	--spreading_rate=0.2 --decay=0.6 --control_increment=0.02 --seed=1692985943

OK_LINE_ARGS = --ridge_type=0 --bumps=1 --ridge_radius=0.2 \
	--activation_types=1 --mutation_type_ub=10 --knob_type=0 --multi_edges=0 \
	--peak_movement=0 --output_types=0 --c2=1 --c3=0 --spreading_rate=0.2 \
	--edge_weights=0 --c1_lb=0.1 --c1_ub=0.9 --extra_mutation_rate=0.05 \
	--decay=0.8 --allow_move_edge=0

#CIRCLE_ARGS = --ridge_type=1 --bumps=1 --ridge_radius=0.1 \
#	--activation_types=1 --mutation_type_ub=10 --knob_type=0 --multi_edges=0 \
#	--peak_movement=0 --output_types=1 --c2=1 --c3=0 --spreading_rate=0.2 \
#	--edge_weights=0 --c1_lb=0.1 --c1_ub=0.9 --extra_mutation_rate=0.05 \
#	--decay=0.8 --allow_move_edge=0

#This produced a near-optimal result
GOOD_XLINE_ARGS = --num_epochs=60 --ridge_type=0 --bumps=0 --ridge_radius=0.2 \
	--activation_types=1 --mutation_type_ub=10 --knob_type=0 --multi_edges=0 \
	--peak_movement=1 --output_types=1 --c2=1 --c3=0 --spreading_rate=0.2 \
	--edge_weights=1 --c1_lb=-1 --c1_ub=1 --extra_mutation_rate=0.00 \
	--decay=0.8 --allow_move_edge=0 --crossover_freq=0.1 --edge_inheritance=4 \
	--num_organisms=200

#This produced a very hill-shaped hill for the y=x line
#Until BEN disallowed initial activations that aren't evenly divisible
#by knob_constant. 4-Apr-2018
GOOD_YXLINE_RUN = --num_epochs=120 --ridge_type=0 --bumps=0 --ridge_radius=0.2 \
	--activation_types=3 --mutation_type_ub=10 --knob_type=0 --multi_edges=0 \
	--peak_movement=0 --output_types=0 --c2=1 --c3=0 --spreading_rate=0.1 \
	--edge_weights=1 --c1_lb=-1 --c1_ub=1 --extra_mutation_rate=0.00 \
	--decay=0.8 --allow_move_edge=0 --crossover_freq=0.1 --edge_inheritance=5 \
	--num_organisms=80 --num_nodes=4 --num_edges=10 --seed=853488368

#This produces a respectable line down the middle on y=x with bumps
OK_YXLINE_BUMPS = #Oops, wrong params

all: sa swig

$(PDFFILES): $(BIBFILES)
%.pdf: %.tex
	$(TEX) $<

%.pdf: %.dot
	$(DOT) < $< > $@

#sds.o: sds.c sds.h sdsalloc.h
#	gcc $(CFLAGS) -c sds.c -o sds.o
#
#sa.o: sa.c sds.o coordset.h sds.h
#	gcc $(CFLAGS) -c sa.c -o sa.o
#
#coordset.o: coordset.c coordset.h
#	gcc $(CFLAGS) -c coordset.c -o coordset.o

sa: sa.o sds.o coordset.o

#	gcc $(CFLAGS) $^ -o sa -lm

data: sa run.py add_param_set.py
	./run.py > d.csv
	./add_param_set.py d.csv > data.csv

out: x
#	./sa > out
#	#grep "epoch fitness" out
#	tail -4 out

okline: sa
	./sa $(OK_LINE_ARGS) > out
	tail -4 out

circle: sa
	./sa $(CIRCLE_ARGS) > out
	tail -4 out

# Experimentation target
x: sa
	./sa $(X_ARGS) > out
	tail -4 out
	@echo

lookatthis: sa
	./sa $(LOOK_AT_THIS) > out
	tail -4 out
	@echo

goodxline: sa
	./sa $(GOOD_XLINE_ARGS) > out
	tail -4 out
	@echo

goodyxline: sa
	./sa $(GOOD_YXLINE_RUN) > out
	tail -4 out
	@echo

outs: sa
	./sa > out
	tail -4 out
	./sa >> out
	tail -4 out
	./sa >> out
	tail -4 out
	./sa >> out
	tail -4 out
	./sa >> out
	tail -4 out

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
	$(VIEW_PDF) test.pdf

phplot:
	sed -n '/BEGIN PHFUNC/,/END PHFUNC/ {//!p;}' out > phfunc
	./plot_xyz.py phfunc 0 1 2

vfplot:
	sed -n '/BEGIN VFUNC/,/END VFUNC/ {//!p;}' out > vfunc
	./plot_xyz.py vfunc 0 1 4

phrangeplot:
	sed -n '/BEGIN VFUNC/,/END VFUNC/ {//!p;}' out > vfunc
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

NRUNS = $(shell seq 2 20)
many: sa
	@./sa $(X_ARGS) --run=1 > /tmp/out
	@tail -4 /tmp/out
	@$(foreach i,$(NRUNS),./sa $(X_ARGS) --run=$i >> /tmp/out; tail -4 /tmp/out;)

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

swig: _sa.so

MACH := $(shell uname)
ifeq ($(MACH),Darwin)
  EXTRA_WARN=""
else
  EXTRA_WARN="-Wno-maybe-uninitialized"
endif

_sa.so: sa.i sa.c
	swig -python sa.i
	CFLAGS="-std=gnu99 -g -Wno-strict-prototypes -Wno-return-type $(EXTRA_WARN) -Wno-unused-variable" LDFLAGS="-lm" python setup_sa.py build_ext --inplace

swig_clean:
	rm -rf *.pyc *.so a.out* build sa_wrap.c* sa.py

swig_test:
	python -c "import sa; sa.tom()"

clean: swig_clean
	rm -f sa *.o test.pdf test.dot

tags:
	ctags *.[ch]

.PHONY: tags run all dot clean plot fitness with_seed out tom swig_clean swig_test runs

# NEW EXPERIMENTS 2019-Mar-10 Inspired by Etienne Barnard's suggestion to look
# at feed-forward networks.

GOOD = --sa_timesteps=20 --log=ancestors --bumps=0 --num_organisms=40 --multi_edges=0 --knob_constant=0.1
X1 = --log=ancestors --sa_timesteps=20 --bumps=1 --moat_ub=0.0 --num_organisms=40 --multi_edges=0 --knob_constant=0.1 --activation_types=5 --alpha=0.0  #--num_epochs=1000
#X = $(OBLIQUE_LINE) --log=ancestors --sa_timesteps=20 --bumps=1 --moat_ub=0.0 --num_organisms=40 --multi_edges=0 --knob_constant=0.1 --activation_types=5 --alpha=0.9  #--num_epochs=1000

YX = $(YXLINE) --ridge_radius=0.1 --bumps=0 \
	--knob_type=1 --knob_constant=0.02 --crossover_freq=0.05 --mutation_type_ub=16 \
	--num_organisms=200 --num_candidates=7 \
	--input_accs=1 --activation_types=6 --sa_timesteps=10 --alpha=0.8 \
	--edge_from_phnode=1 --edge_inheritance=1 --multi_edges=0 \
	--num_epochs=40 --dot=1

GOOD_YXLINE1 = $(YX) --ridge_radius=0.2 --num_organisms=800 \
	--dot=1 --log=ancestors --seed=2848818913

GOOD_YXLINE2 = $(YX) --log=ancestors --seed=2436093377

THINYX_WITH_BUMPS = $(YXLINE) --ridge_radius=0.2 --bumps=1 --moat_ub=0.0 \
	--knob_constant=0.02 --crossover_freq=0.05 --mutation_type_ub=16 --num_organisms=100 \
	--input_accs=1 --activation_types=6 --sa_timesteps=10 --alpha=0.8 \
	--edge_from_phnode=1 --edge_inheritance=5 --multi_edges=0 \
	--num_epochs=40

GOOD_THINYX_WITH_BUMPS = $(YXLINE) $(THINYX_WITH_BUMPS)
	--dot=1 --log=ancestors --num_epochs=40  --seed=1583407075

X1 = $(THINYX_WITH_BUMPS) --ridge_radius=0.1 --moat_ub=0.5 --bump_freq=30.0 \
	--flat=1 --flat_multiplier=1.0 --down_bump=1 \
	--num_organisms=80 --edge_inheritance=5 --knob_type=1 --knob_constant=0.01 \
	--num_candidates=8 --num_epochs=60 --viability_lb=0.0 \
  #--seed=577086024 --log=ancestors

# thu night
RUN1 = $(THINYX_WITH_BUMPS) --moat_ub=0.5 --bump_freq=15.0 \
	--num_organisms=800 --edge_inheritance=5 --knob_type=0 \
	--seed=3233094016 --log=ancestors

CLOSE_BUMPS_ACCLIVATION = $(THINYX_WITH_BUMPS) --moat_ub=0 --bump_freq=15.0 --flat=1 \
	--num_organisms=80 --edge_inheritance=5 --knob_type=0 \
	--num_candidates=8 --seed=3859373065 --log=ancestors

TIGHT_FOLDING_FOR_LEAPING = $(THINYX_WITH_BUMPS) \
	--ridge_radius=0.1 --moat_ub=0.5 --bump_freq=30.0 \
	--flat=1 --flat_multiplier=1.0 --down_bump=1 \
	--num_organisms=80 --edge_inheritance=5 --knob_type=1 --knob_constant=0.01 \
	--num_candidates=8 --num_epochs=60 --viability_lb=0.0 \
	--seed=1207166735 --log=ancestors

OBLIQUE_WITH_BUMPS = $(OBLIQUE_LINE) --ridge_radius=1.0 --bumps=1 --moat_ub=0.0 \
	--knob_constant=0.02 --crossover_freq=0.05 --mutation_type_ub=16 --num_organisms=200 \
	--input_accs=1 --activation_types=6 --sa_timesteps=10 --alpha=0.8 \
	--edge_from_phnode=1 --edge_inheritance=5 --multi_edges=0

OK_OBLIQUE_WITH_BUMPS = $(OBLIQUE_WITH_BUMPS) \
	--dot=1 --log=ancestors --num_epochs=40  --seed=2408275062

X2 = $(OBLIQUE_WITH_BUMPS) --ridge_radius=0.2 \
	--dot=0 --num_epochs=40 --sa_timesteps=3 \
	--mutation_type_ub=16 --num_organisms=400 --generations_per_epoch=20

OK = $(YXLINE) --ridge_radius=0.2 --bumps=1 --moat_ub=0.0 \
	--knob_type=1 --knob_constant=0.02 --crossover_freq=0.02 --mutation_type_ub=100 --num_organisms=400 \
	--input_accs=1 --activation_types=6 --sa_timesteps=5 --alpha=0.8 \
	--edge_from_phnode=0 --edge_inheritance=5 --multi_edges=0 \
	--num_epochs=50  --dot=0 #--log=ancestors #--seed=1583407075

C0 =  $(CIRCLE) --ridge_radius=0.1 --bumps=1 --moat_ub=0.0 \
	--knob_type=1 --knob_constant=0.02 --crossover_freq=0.02 --mutation_type_ub=100 --num_organisms=800 --num_candidates=80 \
	--input_accs=1 --activation_types=6 --sa_timesteps=5 --alpha=0.8 \
	--edge_from_phnode=0 --edge_inheritance=5 --multi_edges=0 \
	--num_epochs=50  --dot=0 #--log=ancestors #--seed=1583407075

ARGS = $(GOOD_YXLINE2) # Change this to some other variable to run other parameters
#ARGS = $(OK_OBLIQUE_WITH_BUMPS)
run: all
	./sa $(ARGS) > out
	@grep 'deltas' out

N = $(shell seq 1 20)
runs: all
	@echo "$(ARGS)"
	@rm -f outs/out*
	@$(foreach i,$(N),./sa $(ARGS) --run=$i > outs/out$i; echo -n "$i: "; grep 'fitness deltas' outs/out$i;)

# Targets for generating plots
plot1:
	echo "\
select 4.4.4\n\
plot phfitness show=False filename=e4g4o4phfit.png azim=30.0\n\
exit\n\
" | ./analyze.py

plot_final_fittest:
	echo "\
select 41.20.77\n\
plot phfitness show=False delta=0.005 azim=34.0 elev=64 cmapname=gnuplot_r\n\
plot vfitness show=False delta=0.005 azim=34.0 elev=64 cmapname=gnuplot_r\n\
plot phrange show=False delta=0.01 azim=34.0 elev=64 cmapname=plasma\n\
dot show=False\n\
exit\n\
" | ./analyze.py
