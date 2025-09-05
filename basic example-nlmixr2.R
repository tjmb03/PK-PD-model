

install.packages(c("dparser", "lotri", "rxode2ll", "rxode2parse",
                   "rxode2random", "rxode2et", "rxode2",
                   "nlmixr2data", "nlmixr2est", "nlmixr2extra",
                   "nlmixr2plot", "nlmixr2"))


library(nlmixr2)

## The basic model consists of an ini block that has initial estimates
one.compartment <- function() {
  ini({
    tka <- log(1.57); label("Ka")
    tcl <- log(2.72); label("Cl")
    tv <- log(31.5); label("V")
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  # and a model block with the error specification and model specification
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

## The fit is performed by the function nlmixr/nlmixr2 specifying the model, data and estimate
fit <- nlmixr2(one.compartment, theo_sd,  est="saem", saemControl(print=0))
print(fit)
