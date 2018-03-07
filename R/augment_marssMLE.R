###############################################################################################################################################
#  augment method for class marssMLE
#  returns fitted values, residuals, std err of residuals and std residuals
##############################################################################################################################################
augment.marssMLE <- function (x, type.predict = c("observations ", "states"),
                              interval = c("none", "confidence"), level = 0.95, ...) {
  model=x[["model"]]
  resids = residuals(x)
  se.resids = sqrt(apply(resids$var.residuals,3,function(x){takediag(x)}))
  data.dims = attr(model, "model.dims")[["y"]]
  nn = data.dims[1]
  TT = data.dims[2]
  if(type.predict=="observations"){
    data.names = attr(model, "Y.names")
    ret = data.frame(
      .rownames=rep(data.names,each=TT), 
      tid=rep(1:TT,nn), 
      y = vec(t(model$data)),
      .fitted = vec(t(fitted(x))),
      .se.fit = vec(t(se.resids[1:nn,,drop=FALSE])),
      .resids = vec(t(resids$model.residuals)),
      .std.resid = vec(t(resids$std.residuals[1:nn,,drop=FALSE]))
    )
    if(interval=="confidence"){
      ret = cbind(
        .conf.low <- qnorm(level/2)*ret$.se.fit + ret$.fitted,
        .conf.up <- qnorm(1-level/2)*ret$.se.fit + ret$.fitted
      )
    }
  }
  if(type.predict=="states"){
    state.names = attr(model, "X.names")
    state.dims = attr(model, "model.dims")[["x"]]
    mm = state.dims[1]
    ret = data.frame(
      .rownames=rep(state.names,each=TT), 
      tid=rep(1:TT,mm), 
      y = vec(t(model$data)),
      .fitted = vec(t(fitted(x))),
      .se.fit = vec(t(se.resids[1:nn,,drop=FALSE])),
      .resids = vec(t(resids$model.residuals)),
      .std.resid = vec(t(resids$std.residuals[1:nn,,drop=FALSE]))
    )
    if(interval=="confidence"){
      ret = cbind(
        .conf.low <- qnorm(level/2)*ret$.se.fit + ret$.fitted,
        .conf.up <- qnorm(1-level/2)*ret$.se.fit + ret$.fitted
      )
    }
  }
  ret
}