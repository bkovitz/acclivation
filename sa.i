%module sa
%{
#include "sa.c"
%}

/* capturing FILE output into a python buffer */

%{
#include <stdio.h>

#ifdef __APPLE__
// funopen uses int instead of size_t
static int py_write(void *cookie, const char *buf, int size) {
  // Note we might need to acquire the GIL here, depending on what you target exactly
  PyObject *result = PyObject_CallMethodObjArgs(cookie, PyString_FromString("write"),
                                                PyString_FromStringAndSize(buf, size), NULL);
  (void)result; // Should we DECREF?
  return size; // assume OK, should really catch instead though
}
#else
static ssize_t py_write(void *cookie, const char *buf, size_t size) {
  // Note we might need to acquire the GIL here, depending on what you target exactly
  PyObject *result = PyObject_CallMethodObjArgs(cookie, PyString_FromString("write"),
                                                PyString_FromStringAndSize(buf, size), NULL);
  (void)result; // Should we DECREF?
  return size; // assume OK, should really catch instead though
}
#endif

static int py_close(void *cookie) {
  Py_DECREF(cookie);
  return 0;
}

static FILE *fopen_python(PyObject *output) {
  if (PyFile_Check(output)) {
    // See notes at: https://docs.python.org/2/c-api/file.html about GIL
    return PyFile_AsFile(output);
  }

  Py_INCREF(output);
#ifdef __APPLE__
  return funopen(output, NULL, py_write, NULL, py_close);
#else
  cookie_io_functions_t funcs = {
    .write = py_write,
    .close = py_close,
  };
  return fopencookie(output, "w", funcs);
#endif
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

#define WITH_SWIG
%include "sa.c"

/* access to arrays inside RUN, EPOCH, GENERATION */

%extend RUN {
  size_t __len__() { return $self->num_epochs; }

  PyObject *__getitem__(size_t i) {
    PyObject *resultobj = 0;
    EPOCH *result = 0 ;
    if (i > ($self->num_epochs-1)) {
      PyErr_SetNone(PyExc_StopIteration);
      return NULL;
    }
    result = (EPOCH *) ($self->epochs[i]);
    resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_EPOCH, 0 |  0 );
    return resultobj;
  }
};

%extend EPOCH {
  size_t __len__() { return $self->num_generations; }

  PyObject *__getitem__(size_t i) {
    PyObject *resultobj = 0;
    GENERATION *result = 0 ;
    if (i > ($self->num_generations-1)) {
      PyErr_SetNone(PyExc_StopIteration);
      return NULL;
    }
    result = (GENERATION *) ($self->generations[i]);
    resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_GENERATION, 0 |  0 );
    return resultobj;
  }
};

%extend GENERATION {
  size_t __len__() { return $self->num_organisms; }

  PyObject *__getitem__(size_t i) {
    PyObject *resultobj = 0;
    Organism *result = 0 ;
    if (i > ($self->num_organisms-1)) {
      PyErr_SetNone(PyExc_StopIteration);
      return NULL;
    }
    result = (Organism *) ($self->organisms[i]);
    resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_Organism, 0 |  0 );
    return resultobj;
  }
};
