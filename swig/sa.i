%module sa
%{
#include "../sa.c"
%}

#define WITH_SWIG
%include "../sa.c"
