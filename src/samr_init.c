// Automatically generated, editing not advised.
#ifndef R_samr_H
#define R_samr_H
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("samr", String)
#else
#define _(String) (String)
#endif

#define FDEF(name)  {#name, (DL_FUNC) &F77_SUB(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
void F77_SUB(rankcol)(
	double *x,
	int *n,
	int *ip,
	int *ixr,
	int *iscrat
);


static R_NativePrimitiveArgType rankcol_t[] = {
	REALSXP,
	INTSXP,
	INTSXP,
	INTSXP,
	INTSXP
};
static R_FortranMethodDef fMethods[] = {
	FDEF(rankcol),
	{NULL, NULL, 0}
};

void R_init_samr(DllInfo *dll){
  R_registerRoutines(dll, NULL, NULL, fMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

#endif
