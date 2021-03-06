## This file was automatically generated by odin.
## Do not edit by hand as changes will be lost.
.R6_larvalMod <- R6::R6Class(
  "odin_model",
  parent_env = environment(odin::odin),
  public = list(
    name = "larvalMod",
    ptr = NULL,
    ## Cache:
    init = NULL,
    variable_order = NULL,
    names = NULL,
    transform_variables = NULL,
    ## Methods:
    initialize = function(user = NULL) {
      self$ptr <- .Call("larvalMod_create", user, PACKAGE = "odinPackage")
      self$init <- .Call("larvalMod_initialise", self$ptr, NA_real_, PACKAGE = "odinPackage")
      loadNamespace("dde")
      self$update_cache()
    },

    set_user = function(..., user = list(...)) {
      .Call("r_larvalMod_set_user", self$ptr, user, PACKAGE = "odinPackage")
      self$update_cache()
      invisible(self$init)
    },

    update_cache = function() {
      self$variable_order <- .Call("larvalMod_variable_order", self$ptr, PACKAGE = "odinPackage")
      odin_prepare(self, TRUE)
    },

    update = function(step, y) {
      .Call("larvalMod_update_r", self$ptr, step, y, PACKAGE = "odinPackage")
    },

    initial = function(step) {
      self$init
    },

    run = function(step, y = NULL, ..., use_names = TRUE) {
      if (is.null(y)) {
        y <- self$init
      }
      ret <- dde::difeq(y, step, "larvalMod_update_dde", self$ptr,
                        dllname = "odinPackage",
                        parms_are_real = FALSE, ynames = FALSE, ...)
      if (use_names) {
        colnames(ret) <- self$names
      } else {
        colnames(ret) <- NULL
      }
      ret
    },

    contents = function() {
      .Call("larvalMod_contents", self$ptr, PACKAGE = "odinPackage")
    }
  ))
larvalMod <- function(dE, dL, dP, uoE, uoL, uP, uM, Y, S, tr, sf, dt, n, Emax, rF, user = list(dE = dE, dL = dL, dP = dP, uoE = uoE, uoL = uoL, uP = uP, uM = uM, Y = Y, S = S, tr = tr, sf = sf, dt = dt, n = n, Emax = Emax, rF = rF)) {
  .R6_larvalMod$new(user = user)
}
.R6_larvalModDet <- R6::R6Class(
  "odin_model",
  parent_env = environment(odin::odin),
  public = list(
    name = "larvalModDet",
    ptr = NULL,
    use_dde = NULL,
    ## Cache:
    init = NULL,
    variable_order = NULL,
    names = NULL,
    transform_variables = NULL,
    ## Methods:
    initialize = function(user = NULL, use_dde = FALSE) {
      self$ptr <- .Call("larvalModDet_create", user, use_dde, PACKAGE = "odinPackage")
      self$use_dde <- use_dde
      self$init <- .Call("larvalModDet_initialise", self$ptr, NA_real_, PACKAGE = "odinPackage")
      if (use_dde) {
        loadNamespace("dde")
      }
      self$update_cache()
    },

    set_user = function(..., user = list(...)) {
      .Call("r_larvalModDet_set_user", self$ptr, user, PACKAGE = "odinPackage")
      self$update_cache()
      invisible(self$init)
    },

    update_cache = function() {
      self$variable_order <- .Call("larvalModDet_variable_order", self$ptr, PACKAGE = "odinPackage")
      odin_prepare(self, FALSE)
    },

    deriv = function(t, y) {
      .Call("larvalModDet_deriv_r", self$ptr, t, y, PACKAGE = "odinPackage")
    },

    initial = function(t) {
      self$init
    },

    run = function(t, y = NULL, ..., use_names = TRUE) {
      if (is.null(y)) {
        y <- self$init
      }
      if (self$use_dde) {
        ret <- dde::dopri(y, t, "larvalModDet_deriv_dde", self$ptr,
                          dllname = "odinPackage",
                          parms_are_real = FALSE, ynames = FALSE, ...)
      } else {
        ret <- deSolve::ode(y, t, "larvalModDet_deriv_ds", self$ptr,
                            initfunc = "larvalModDet_initmod_ds", dllname = "odinPackage",
                            ...)
      }
      if (use_names) {
        colnames(ret) <- self$names
      } else {
        colnames(ret) <- NULL
      }
      ret
    },

    contents = function() {
      .Call("larvalModDet_contents", self$ptr, PACKAGE = "odinPackage")
    }
  ))
larvalModDet <- function(dE, dL, dP, uoE, uoL, uP, uM, Y, S, tr, sf, dt, B, Emax, rF, user = list(dE = dE, dL = dL, dP = dP, uoE = uoE, uoL = uoL, uP = uP, uM = uM, Y = Y, S = S, tr = tr, sf = sf, dt = dt, B = B, Emax = Emax, rF = rF), use_dde = FALSE) {
  .R6_larvalModDet$new(user = user, use_dde = use_dde)
}
.R6_larvalModP <- R6::R6Class(
  "odin_model",
  parent_env = environment(odin::odin),
  public = list(
    name = "larvalModP",
    ptr = NULL,
    ## Cache:
    init = NULL,
    variable_order = NULL,
    names = NULL,
    transform_variables = NULL,
    ## Methods:
    initialize = function(user = NULL) {
      self$ptr <- .Call("larvalModP_create", user, PACKAGE = "odinPackage")
      self$init <- .Call("larvalModP_initialise", self$ptr, NA_real_, PACKAGE = "odinPackage")
      loadNamespace("dde")
      self$update_cache()
    },

    set_user = function(..., user = list(...)) {
      .Call("r_larvalModP_set_user", self$ptr, user, PACKAGE = "odinPackage")
      self$init <- .Call("larvalModP_initialise", self$ptr, NA_real_, PACKAGE = "odinPackage")
      self$update_cache()
      invisible(self$init)
    },

    update_cache = function() {
      self$variable_order <- .Call("larvalModP_variable_order", self$ptr, PACKAGE = "odinPackage")
      odin_prepare(self, TRUE)
    },

    update = function(step, y) {
      .Call("larvalModP_update_r", self$ptr, step, y, PACKAGE = "odinPackage")
    },

    initial = function(step) {
      self$init
    },

    run = function(step, y = NULL, ..., use_names = TRUE) {
      if (is.null(y)) {
        y <- self$init
      }
      ret <- dde::difeq(y, step, "larvalModP_update_dde", self$ptr,
                        dllname = "odinPackage",
                        parms_are_real = FALSE, ynames = FALSE, ...)
      if (use_names) {
        colnames(ret) <- self$names
      } else {
        colnames(ret) <- NULL
      }
      ret
    },

    contents = function() {
      .Call("larvalModP_contents", self$ptr, PACKAGE = "odinPackage")
    }
  ))
larvalModP <- function(dE, dL, dP, uoE, uoL, uP, uM, Y, S, tr, sf, dt, n, Emax, E0, L0, P0, M0, rF, user = list(dE = dE, dL = dL, dP = dP, uoE = uoE, uoL = uoL, uP = uP, uM = uM, Y = Y, S = S, tr = tr, sf = sf, dt = dt, n = n, Emax = Emax, E0 = E0, L0 = L0, P0 = P0, M0 = M0, rF = rF)) {
  .R6_larvalModP$new(user = user)
}
