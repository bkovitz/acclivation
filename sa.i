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

%{
//#define _GNU_SOURCE - not needed, Python already does that!
#include <stdio.h>

static ssize_t py_write(void *cookie, const char *buf, size_t size) {
  // Note we might need to acquire the GIL here, depending on what you target exactly
  PyObject *result = PyObject_CallMethodObjArgs(cookie, PyString_FromString("write"),
                                                PyString_FromStringAndSize(buf, size), NULL);

  (void)result; // Should we DECREF?
  return size; // assume OK, should really catch instead though
}

static int py_close(void *cookie) {
  Py_DECREF(cookie);
  return 0;
}

static FILE *fopen_python(PyObject *output) {
  if (PyFile_Check(output)) {
    // See notes at: https://docs.python.org/2/c-api/file.html about GIL
    return PyFile_AsFile(output);
  }

  cookie_io_functions_t funcs = {
    .write = py_write,
    .close = py_close,
  };
  Py_INCREF(output);
  return fopencookie(output, "w", funcs);
}
%}

%typemap(in) FILE * {
  $1 = fopen_python($input);
}

%typemap(freearg) FILE * {
  // Note GIL comment above here also
  // fileno for fopencookie always returns -1
  if (-1 == fileno($1)) fclose($1);
}

%inline %{
  void hello(FILE *out)
  {
    fprintf(out, "Hello How are you\n");
  }
%}
