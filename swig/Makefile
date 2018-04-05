all:
	swig -python sa.i
	CFLAGS="-std=gnu99" LDFLAGS="-lm" python setup_sa.py build_ext --inplace

clean:
	rm -rf *.pyc *.so a.out* build sa_wrap.c* sa.py

test:
	python -c "import sa; sa.tom()"
