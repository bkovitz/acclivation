%module sa
%{
#include "sa.c"
%}

#define WITH_SWIG
%include "sa.c"

%extend RUN {
  size_t __len__() { return $self->num_epochs; }

  PyObject *__getitem__(size_t i) {
    PyObject *resultobj = 0;
    EPOCH *result = 0 ;
    if (i < 0) {
      SWIG_exception_fail(SWIG_ArgError(0), "index out of range error");
    } else if (i > ($self->num_epochs-1)) {
      PyErr_SetNone(PyExc_StopIteration);
      return NULL;
    }
    result = (EPOCH *) ($self->epochs[i]);
    resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_EPOCH, 0 |  0 );
    return resultobj;
  fail:
    return NULL;
  }
};

%extend EPOCH {
  size_t __len__() { return $self->num_generations; }

  PyObject *__getitem__(size_t i) {
    PyObject *resultobj = 0;
    GENERATION *result = 0 ;
    if (i < 0) {
      SWIG_exception_fail(SWIG_ArgError(0), "index out of range error");
    } else if (i > ($self->num_generations-1)) {
      PyErr_SetNone(PyExc_StopIteration);
      return NULL;
    }
    result = (GENERATION *) ($self->generations[i]);
    resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_GENERATION, 0 |  0 );
    return resultobj;
  fail:
    return NULL;
  }
};

%extend GENERATION {
  size_t __len__() { return $self->num_organisms; }

  PyObject *__getitem__(size_t i) {
    PyObject *resultobj = 0;
    Organism *result = 0 ;
    if (i < 0) {
      SWIG_exception_fail(SWIG_ArgError(0), "index out of range error");
    } else if (i > ($self->num_organisms-1)) {
      PyErr_SetNone(PyExc_StopIteration);
      return NULL;
    }
    result = (Organism *) ($self->organisms[i]);
    resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_Organism, 0 |  0 );
    return resultobj;
  fail:
    return NULL;
  }
};
