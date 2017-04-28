// This file was automatically generated by odin.
// Do not edit by hand as changes will be lost.
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdbool.h>

// Collect together all the parameters and transient memory
// required to run the model in a struct.
typedef struct larvalMod_pars {
  double initial_Be;
  double initial_Bl;
  double initial_Bp;
  double initial_E;
  double initial_L;
  double initial_P;
  double initial_M;
  double initial_Bm;
  double initial_nt;
  double initial_Reff;
  double dE;
  double dL;
  double dP;
  double uoE;
  double uoL;
  double uP;
  double uM;
  double Y;
  double S;
  double tr;
  double sf;
  double dt;
  double n;
  int dim_rF;
  double *rF;
  double Emax;
  double trx;
} larvalMod_pars;
larvalMod_pars* larvalMod_get_pointer(SEXP larvalMod_ptr, int closed_error);
SEXP larvalMod_set_user(larvalMod_pars *larvalMod_p, SEXP user);

SEXP get_ds_pars();
double get_user_double(SEXP user, const char *name, double default_value);
double* get_user_array_dim(SEXP user, const char *name, bool is_real, int nd, int *dest_dim);
SEXP get_user_array_check_rank(SEXP user, const char *name, int nd);
void get_user_array_copy(SEXP el, const char *name, bool is_real, void *dest);
SEXP get_list_element(SEXP list, const char *name);
double odin_sum1(double *x, int from_i, int to_i);

// Create the pointer; this will establish the struct, allocate
// memory for things that are constant size, and initialize
// constant variables
static void larvalMod_finalize(SEXP larvalMod_ptr);
SEXP larvalMod_create(SEXP user) {
  larvalMod_pars *larvalMod_p = (larvalMod_pars*) Calloc(1, larvalMod_pars);
  larvalMod_p->dE = NA_REAL;
  larvalMod_p->dL = NA_REAL;
  larvalMod_p->dP = NA_REAL;
  larvalMod_p->uoE = NA_REAL;
  larvalMod_p->uoL = NA_REAL;
  larvalMod_p->uP = NA_REAL;
  larvalMod_p->uM = NA_REAL;
  larvalMod_p->Y = NA_REAL;
  larvalMod_p->S = NA_REAL;
  larvalMod_p->tr = NA_REAL;
  larvalMod_p->sf = NA_REAL;
  larvalMod_p->dt = NA_REAL;
  larvalMod_p->n = NA_REAL;
  larvalMod_p->rF = NULL;
  larvalMod_p->Emax = NA_REAL;
  SEXP larvalMod_ptr = PROTECT(R_MakeExternalPtr(larvalMod_p, R_NilValue, R_NilValue));
  R_RegisterCFinalizer(larvalMod_ptr, larvalMod_finalize);
  larvalMod_set_user(larvalMod_p, user);
  UNPROTECT(1);
  return larvalMod_ptr;
}

// Set user-supplied parameter values.
SEXP larvalMod_set_user(larvalMod_pars *larvalMod_p, SEXP user) {
  larvalMod_p->dE = get_user_double(user, "dE", larvalMod_p->dE);
  larvalMod_p->dL = get_user_double(user, "dL", larvalMod_p->dL);
  larvalMod_p->dP = get_user_double(user, "dP", larvalMod_p->dP);
  larvalMod_p->uoE = get_user_double(user, "uoE", larvalMod_p->uoE);
  larvalMod_p->uoL = get_user_double(user, "uoL", larvalMod_p->uoL);
  larvalMod_p->uP = get_user_double(user, "uP", larvalMod_p->uP);
  larvalMod_p->uM = get_user_double(user, "uM", larvalMod_p->uM);
  larvalMod_p->Y = get_user_double(user, "Y", larvalMod_p->Y);
  larvalMod_p->S = get_user_double(user, "S", larvalMod_p->S);
  larvalMod_p->tr = get_user_double(user, "tr", larvalMod_p->tr);
  larvalMod_p->sf = get_user_double(user, "sf", larvalMod_p->sf);
  larvalMod_p->dt = get_user_double(user, "dt", larvalMod_p->dt);
  larvalMod_p->n = get_user_double(user, "n", larvalMod_p->n);
  Free(larvalMod_p->rF);
  {
    double *tmp = (double*)get_user_array_dim(user, "rF", true, 1, &larvalMod_p->dim_rF);
    larvalMod_p->rF = (double*) Calloc(larvalMod_p->dim_rF, double);
    memcpy(larvalMod_p->rF, tmp, larvalMod_p->dim_rF * sizeof(double));
  }
  larvalMod_p->Emax = get_user_double(user, "Emax", larvalMod_p->Emax);
  larvalMod_p->trx = (larvalMod_p->tr / larvalMod_p->dt);
  return R_NilValue;
}
// Wrapper around this for use from R.
SEXP r_larvalMod_set_user(SEXP larvalMod_ptr, SEXP user) {
  larvalMod_pars *larvalMod_p = larvalMod_get_pointer(larvalMod_ptr, 1);
  larvalMod_set_user(larvalMod_p, user);
  return R_NilValue;
}

// Arrange to free all memory we have allocated
// This is called by R automatically when the pointer is
// garbage collected (i.e., when all objects holding the pointer
// go out of scope
void larvalMod_finalize(SEXP larvalMod_ptr) {
  larvalMod_pars *larvalMod_p = larvalMod_get_pointer(larvalMod_ptr, 0);
  if (larvalMod_ptr) {
    Free(larvalMod_p->rF);
    Free(larvalMod_p);
    R_ClearExternalPtr(larvalMod_ptr);
  }
}

SEXP larvalMod_initialise(SEXP larvalMod_ptr, SEXP step_ptr) {
  larvalMod_pars *larvalMod_p = larvalMod_get_pointer(larvalMod_ptr, 1);
  GetRNGstate();
  larvalMod_p->initial_Be = 0;
  larvalMod_p->initial_Bl = 0;
  larvalMod_p->initial_Bp = 0;
  larvalMod_p->initial_E = 31;
  larvalMod_p->initial_L = 2;
  larvalMod_p->initial_P = 0.34000000000000002;
  larvalMod_p->initial_M = 2.5;
  larvalMod_p->initial_Bm = 0;
  larvalMod_p->initial_nt = 0;
  larvalMod_p->initial_Reff = 0;
  SEXP state = PROTECT(allocVector(REALSXP, 10));
  REAL(state)[0] = larvalMod_p->initial_Be;
  REAL(state)[1] = larvalMod_p->initial_Bl;
  REAL(state)[2] = larvalMod_p->initial_Bp;
  REAL(state)[3] = larvalMod_p->initial_Bm;
  REAL(state)[4] = larvalMod_p->initial_nt;
  REAL(state)[5] = larvalMod_p->initial_Reff;
  REAL(state)[6] = larvalMod_p->initial_E;
  REAL(state)[7] = larvalMod_p->initial_L;
  REAL(state)[8] = larvalMod_p->initial_P;
  REAL(state)[9] = larvalMod_p->initial_M;
  PutRNGstate();
  UNPROTECT(1);
  return state;
}

SEXP larvalMod_set_initial(SEXP larvalMod_ptr, SEXP step_ptr, SEXP state_ptr) {
  return R_NilValue;
}

void larvalMod_update(larvalMod_pars *larvalMod_p, size_t step, double *state, double *state_next, double *output) {
  double Be = state[0];
  double Bl = state[1];
  double Bp = state[2];
  double Bm = state[3];
  double nt = state[4];
  double E = state[6];
  double L = state[7];
  double P = state[8];
  double M = state[9];
  double K = (step <= larvalMod_p->trx ? (1 + (larvalMod_p->sf * (((double) 1 / larvalMod_p->trx) * (0)))) : (1 + (larvalMod_p->sf * (((double) 1 / larvalMod_p->trx) * (odin_sum1(larvalMod_p->rF, (step - larvalMod_p->trx) - 1, step - 1))))));
  double uE = larvalMod_p->uoE * (1 + ((E + L) / (K)));
  double uL = larvalMod_p->uoL * (1 + (larvalMod_p->Y * (E + L) / (K)));
  state_next[0] = Rf_rbinom(round(E), (larvalMod_p->dE + uE) * larvalMod_p->dt);
  state_next[1] = Rf_rbinom(round(L), (larvalMod_p->dL + uL) * larvalMod_p->dt);
  state_next[2] = Rf_rbinom(round(P), (larvalMod_p->dP + larvalMod_p->uP) * larvalMod_p->dt);
  state_next[3] = Rf_rbinom(round(M), larvalMod_p->uM * larvalMod_p->dt);
  state_next[4] = Rf_rbinom(round(M), (larvalMod_p->dt / larvalMod_p->S));
  state_next[5] = 0.5 * (larvalMod_p->Emax / (exp(larvalMod_p->uM * larvalMod_p->S) - 1)) * ((double) 1 / (1 + uE * larvalMod_p->dE)) * ((double) 1 / (1 + uL * larvalMod_p->dL)) * ((double) 1 / (1 + larvalMod_p->uP * larvalMod_p->dP));
  state_next[6] = (E - Be > 0 ? E - Be + Rf_rpois(nt * larvalMod_p->n) : Rf_rpois(nt * larvalMod_p->n));
  state_next[7] = (L - Bl > 0 ? L - Bl + Rf_rbinom(round(Be), (larvalMod_p->dE / (uE + larvalMod_p->dE))) : Rf_rbinom(round(Be), (larvalMod_p->dE / (uE + larvalMod_p->dE))));
  state_next[8] = (P - Bp > 0 ? P - Bp + Rf_rbinom(round(Bl), (larvalMod_p->dL / (uL + larvalMod_p->dL))) : Rf_rbinom(round(Bl), (larvalMod_p->dL / (uL + larvalMod_p->dL))));
  state_next[9] = (M - Bm > 0 ? M + (0.5 * (Rf_rbinom(round(Bp), (larvalMod_p->dP / (larvalMod_p->uP + larvalMod_p->dP))))) - Bm : M + (0.5 * (Rf_rbinom(round(Bp), (larvalMod_p->dP / (larvalMod_p->uP + larvalMod_p->dP))))));
}

void larvalMod_update_dde(size_t n, size_t step, double *state, double *state_next,
               size_t n_out, double *output, void *larvalMod_p) {
  larvalMod_update((larvalMod_pars*)larvalMod_p, step, state, state_next, output);
}

SEXP larvalMod_update_r(SEXP larvalMod_ptr, SEXP step, SEXP state) {
  SEXP state_next = PROTECT(allocVector(REALSXP, LENGTH(state)));
  larvalMod_pars *larvalMod_p = larvalMod_get_pointer(larvalMod_ptr, 1);
  double *output = NULL;
  larvalMod_update(larvalMod_p, INTEGER(step)[0], REAL(state), REAL(state_next), output);
  UNPROTECT(1);
  return state_next;
}

// Translate all elements in the struct back to R
// This will mostly be useful for debugging.
SEXP larvalMod_contents(SEXP larvalMod_ptr) {
  larvalMod_pars *larvalMod_p = larvalMod_get_pointer(larvalMod_ptr, 1);
  SEXP state = PROTECT(allocVector(VECSXP, 27));
  SET_VECTOR_ELT(state, 0, ScalarReal(larvalMod_p->initial_Be));
  SET_VECTOR_ELT(state, 1, ScalarReal(larvalMod_p->initial_Bl));
  SET_VECTOR_ELT(state, 2, ScalarReal(larvalMod_p->initial_Bp));
  SET_VECTOR_ELT(state, 3, ScalarReal(larvalMod_p->initial_E));
  SET_VECTOR_ELT(state, 4, ScalarReal(larvalMod_p->initial_L));
  SET_VECTOR_ELT(state, 5, ScalarReal(larvalMod_p->initial_P));
  SET_VECTOR_ELT(state, 6, ScalarReal(larvalMod_p->initial_M));
  SET_VECTOR_ELT(state, 7, ScalarReal(larvalMod_p->initial_Bm));
  SET_VECTOR_ELT(state, 8, ScalarReal(larvalMod_p->initial_nt));
  SET_VECTOR_ELT(state, 9, ScalarReal(larvalMod_p->initial_Reff));
  SET_VECTOR_ELT(state, 10, ScalarReal(larvalMod_p->dE));
  SET_VECTOR_ELT(state, 11, ScalarReal(larvalMod_p->dL));
  SET_VECTOR_ELT(state, 12, ScalarReal(larvalMod_p->dP));
  SET_VECTOR_ELT(state, 13, ScalarReal(larvalMod_p->uoE));
  SET_VECTOR_ELT(state, 14, ScalarReal(larvalMod_p->uoL));
  SET_VECTOR_ELT(state, 15, ScalarReal(larvalMod_p->uP));
  SET_VECTOR_ELT(state, 16, ScalarReal(larvalMod_p->uM));
  SET_VECTOR_ELT(state, 17, ScalarReal(larvalMod_p->Y));
  SET_VECTOR_ELT(state, 18, ScalarReal(larvalMod_p->S));
  SET_VECTOR_ELT(state, 19, ScalarReal(larvalMod_p->tr));
  SET_VECTOR_ELT(state, 20, ScalarReal(larvalMod_p->sf));
  SET_VECTOR_ELT(state, 21, ScalarReal(larvalMod_p->dt));
  SET_VECTOR_ELT(state, 22, ScalarReal(larvalMod_p->n));
  SET_VECTOR_ELT(state, 23, ScalarInteger(larvalMod_p->dim_rF));
  SET_VECTOR_ELT(state, 24, allocVector(REALSXP, larvalMod_p->dim_rF));
  memcpy(REAL(VECTOR_ELT(state, 24)), larvalMod_p->rF, larvalMod_p->dim_rF * sizeof(double));
  SET_VECTOR_ELT(state, 25, ScalarReal(larvalMod_p->Emax));
  SET_VECTOR_ELT(state, 26, ScalarReal(larvalMod_p->trx));
  SEXP state_names = PROTECT(allocVector(STRSXP, 27));
  SET_STRING_ELT(state_names, 0, mkChar("initial_Be"));
  SET_STRING_ELT(state_names, 1, mkChar("initial_Bl"));
  SET_STRING_ELT(state_names, 2, mkChar("initial_Bp"));
  SET_STRING_ELT(state_names, 3, mkChar("initial_E"));
  SET_STRING_ELT(state_names, 4, mkChar("initial_L"));
  SET_STRING_ELT(state_names, 5, mkChar("initial_P"));
  SET_STRING_ELT(state_names, 6, mkChar("initial_M"));
  SET_STRING_ELT(state_names, 7, mkChar("initial_Bm"));
  SET_STRING_ELT(state_names, 8, mkChar("initial_nt"));
  SET_STRING_ELT(state_names, 9, mkChar("initial_Reff"));
  SET_STRING_ELT(state_names, 10, mkChar("dE"));
  SET_STRING_ELT(state_names, 11, mkChar("dL"));
  SET_STRING_ELT(state_names, 12, mkChar("dP"));
  SET_STRING_ELT(state_names, 13, mkChar("uoE"));
  SET_STRING_ELT(state_names, 14, mkChar("uoL"));
  SET_STRING_ELT(state_names, 15, mkChar("uP"));
  SET_STRING_ELT(state_names, 16, mkChar("uM"));
  SET_STRING_ELT(state_names, 17, mkChar("Y"));
  SET_STRING_ELT(state_names, 18, mkChar("S"));
  SET_STRING_ELT(state_names, 19, mkChar("tr"));
  SET_STRING_ELT(state_names, 20, mkChar("sf"));
  SET_STRING_ELT(state_names, 21, mkChar("dt"));
  SET_STRING_ELT(state_names, 22, mkChar("n"));
  SET_STRING_ELT(state_names, 23, mkChar("dim_rF"));
  SET_STRING_ELT(state_names, 24, mkChar("rF"));
  SET_STRING_ELT(state_names, 25, mkChar("Emax"));
  SET_STRING_ELT(state_names, 26, mkChar("trx"));
  setAttrib(state, R_NamesSymbol, state_names);
  UNPROTECT(2);
  return state;
}

// Report back to R information on variable ordering
// The reported information includes position and length of each
// variable, from which offset, etc, can be worked out.
SEXP larvalMod_variable_order(SEXP larvalMod_ptr) {
  SEXP state_len = PROTECT(allocVector(VECSXP, 10));
  SEXP state_names = PROTECT(allocVector(STRSXP, 10));
  SET_VECTOR_ELT(state_len, 0, R_NilValue);
  SET_STRING_ELT(state_names, 0, mkChar("Be"));
  SET_VECTOR_ELT(state_len, 1, R_NilValue);
  SET_STRING_ELT(state_names, 1, mkChar("Bl"));
  SET_VECTOR_ELT(state_len, 2, R_NilValue);
  SET_STRING_ELT(state_names, 2, mkChar("Bp"));
  SET_VECTOR_ELT(state_len, 3, R_NilValue);
  SET_STRING_ELT(state_names, 3, mkChar("Bm"));
  SET_VECTOR_ELT(state_len, 4, R_NilValue);
  SET_STRING_ELT(state_names, 4, mkChar("nt"));
  SET_VECTOR_ELT(state_len, 5, R_NilValue);
  SET_STRING_ELT(state_names, 5, mkChar("Reff"));
  SET_VECTOR_ELT(state_len, 6, R_NilValue);
  SET_STRING_ELT(state_names, 6, mkChar("E"));
  SET_VECTOR_ELT(state_len, 7, R_NilValue);
  SET_STRING_ELT(state_names, 7, mkChar("L"));
  SET_VECTOR_ELT(state_len, 8, R_NilValue);
  SET_STRING_ELT(state_names, 8, mkChar("P"));
  SET_VECTOR_ELT(state_len, 9, R_NilValue);
  SET_STRING_ELT(state_names, 9, mkChar("M"));
  setAttrib(state_len, R_NamesSymbol, state_names);
  UNPROTECT(2);
  return state_len;
}

larvalMod_pars* larvalMod_get_pointer(SEXP larvalMod_ptr, int closed_error) {
  larvalMod_pars *larvalMod_p = NULL;
  if (TYPEOF(larvalMod_ptr) != EXTPTRSXP) {
    Rf_error("Expected an external pointer");
  }
  larvalMod_p = (larvalMod_pars*) R_ExternalPtrAddr(larvalMod_ptr);
  if (!larvalMod_p && closed_error) {
    Rf_error("Pointer has been invalidated");
  }
  return larvalMod_p;
}

SEXP get_ds_pars() {
  static DL_FUNC get_deSolve_gparms = NULL;
  if (get_deSolve_gparms == NULL) {
    get_deSolve_gparms = R_GetCCallable("deSolve", "get_deSolve_gparms");
  }
  return get_deSolve_gparms();
}
double get_user_double(SEXP user, const char *name, double default_value) {
  double ret = default_value;
  SEXP el = get_list_element(user, name);
  if (el != R_NilValue) {
    if (length(el) != 1) {
      Rf_error("Expected scalar numeric for %s", name);
    }
    if (TYPEOF(el) == REALSXP) {
      ret = REAL(el)[0];
    } else if (TYPEOF(el) == INTSXP) {
      ret = INTEGER(el)[0];
    } else {
      Rf_error("Expected a numeric value for %s", name);
    }
  }
  if (ISNA(ret)) {
    Rf_error("Expected value for %s", name);
  }
  return ret;
}
double* get_user_array_dim(SEXP user, const char *name, bool is_real, int nd, int *dest_dim) {
  SEXP el = get_user_array_check_rank(user, name, nd);

  if (nd == 1) {
    dest_dim[0] = LENGTH(el);
  } else {
    SEXP r_dim = PROTECT(coerceVector(getAttrib(el, R_DimSymbol), INTSXP));
    int *dim = INTEGER(r_dim);

    for (size_t i = 0; i < (size_t) nd; ++i) {
      dest_dim[i] = dim[i];
    }

    UNPROTECT(1);
  }

  void *dest;
  if (is_real) {
    dest = (double*) Calloc(length(el), double);
  } else {
    dest = (int*) Calloc(length(el), int);
  }
  get_user_array_copy(el, name, is_real, dest);
  return dest;
}
SEXP get_user_array_check_rank(SEXP user, const char *name, int nd) {
  SEXP el = get_list_element(user, name);
  if (el == R_NilValue) {
    Rf_error("Expected value for %s", name);
  } else {
    if (nd == 1) {
      if (isArray(el)) {
        // this may be too strict as a length-1 dim here will fail
        Rf_error("Expected a vector for %s", name);
      }
    } else {
      SEXP r_dim = getAttrib(el, R_DimSymbol);
      if (r_dim == R_NilValue || LENGTH(r_dim) != nd) {
        if (nd == 2) {
          Rf_error("Expected a matrix for %s", name);
        } else {
          Rf_error("Expected a %dd array for %s", nd, name);
        }
      }
    }
  }
  return el;
}
void get_user_array_copy(SEXP el, const char *name, bool is_real, void *dest) {
  int given_int = TYPEOF(el) == INTSXP;
  size_t n = (size_t) length(el);
  if (is_real) {
    if (given_int) {
      el = PROTECT(coerceVector(el, REALSXP));
    } else if (TYPEOF(el) != REALSXP) {
      Rf_error("Expected a numeric value for %s", name);
    }
    memcpy(dest, REAL(el), n * sizeof(double));
  } else {
    if (TYPEOF(el) == REALSXP) {
      el = PROTECT(coerceVector(el, INTSXP));
    } else if (TYPEOF(el) != INTSXP) {
      Rf_error("Expected a numeric value for %s", name);
    }
    memcpy(dest, INTEGER(el), n * sizeof(int));
  }
  if (given_int == is_real) {
    UNPROTECT(1);
  }
}
SEXP get_list_element(SEXP list, const char *name) {
  SEXP ret = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); ++i) {
    if(strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      ret = VECTOR_ELT(list, i);
      break;
    }
  }
  return ret;
}
double odin_sum1(double *x, int from_i, int to_i) {
  double tot = 0.0;
  for (int i = from_i; i <= to_i; ++i) {
    tot += x[i];
  }
  return tot;
}