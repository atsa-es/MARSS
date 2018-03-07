MARSSinfo = function(number){
  if(missing(number)){
cat("Pass in a label (in quotes) to get info on a MARSS error or warning message.
     convergence: Non-convergence warnings
     denominv: An error related to denom not invertible
     degenvarcov: Warnings about degenerate variance-covariance matrices or variance going to 0
     diag0blocked: Warning: setting diagonal to 0 blocked at iter=X. logLik was lower in attempt to set 0 diagonals on X. See also R0blocked
     HessianNA: Warning: There are NAs in the Hessian matrix.
     is.marssMLE: Error from is.marssMLE
     kferrors: Stopped at iter=xx in MARSSkem() because numerical errors were generated in MARSSkf
     LLdropped: MARSS warns that log-likelihood dropped
     LLunstable: iter=xx MARSSkf: logLik computation is becoming unstable.  Condition num. of Sigma[t=1] = Inf and of R = Inf.
     modelclass: Your model object is not the right class.
     residvarinv: warning: the variance of the residuals at t = x is not invertible.
     R0blocked: Setting of 0s on the diagonal of R blocked; corresponding x0 should not be estimated.  See also x0R0 error and diag0blocked.
     slowconvergence: MARSS seems to take a long, long, long time to converge
     ts: MARSS complains that you passed in a ts object as data
     varcovstruc: Error: MARSS says the variance-covariance matrix is illegal.
     x0R0: Error concerning setting of x0 in model with R with 0s on diagonal
     V0init: MARSS complains about init values for V0
")
return()
}
if(number=="V0init")
cat( writeLines(strwrap(
"In a variance-covariance matrix, you cannot have 0s on the diagonal. When you pass in a qxq variance-covariance matrix (Q,R or V0) with a 0 on the diagonal, MARSS rewrites this into H%*%V.  V is a pxp variance-covariance matrix made from only the non-zero row/cols in your variance-covariance matrix and H is a q x p matrix with all zero rows corresponding to the 0 diagonals in the variance-covariance matrix that you passed in. By setting, the start  variance to 0, you have forced a 0 on the diagonal of a variance-covariance matrix and that will break the EM algorithm (the BFGS algorithm will probably work since it uses start in a different way). This is probably just an accident in how you specified your start values for your variance-covariance matrices. Check the start values and compare to the model$free values.  If you did not pass in start, then MARSS's function for generating reasonable start values does not work for your model.  So just pass in your own start values for the variances. Note, matrices with a second dim equal to 0 are fine in test$start (and test$par). It just means the parameter is fixed (not estimated).\n"
)))

if(number=="is.marssMLE")
cat( writeLines(strwrap(
"If you got an error from is.marssMLE related to your model, then the first thing to do is look at the list that you passed into the model argument.  Often that will reveal the problem.  If not, then look at your data and make sure it is a nxT matrix (and not a Txn) and doesn't have any weird values in it.  If the problem is still not clearyou need to look at the model that MARSS thinks you are trying to fit.
\n\n
If you used test=MARSS(foo), then test is the MLE object.  If the function exited without giving you the MLE object, try test=MARSS(...,fit=FALSE) to get it.  Type summary(test$model) to see a print out of the model.  If your model is time-varying, this will be very verbose so you'll want to divert the output to a file.  Then try this test$par=test$start, now you have filled in the par element of the MLE object. Try parmat(test,t=1) to see all the parameters at t=1 using the start as the par values.  This might reveal the problem too.  Note, matrices with a second dim equal to 0 are fine in test$start (and test$par). It just means the parameter is fixed (not estimated).\n"
)))

if(number=="denominv")
  cat( writeLines(strwrap(
"This is telling you that you specifified a model that is logically indeterminant. First check your data and covariates (if you have them).  Make sure you didn't make a mistake when entering the data.  For example, a row of data that is all NAs or two rows of c or dthat are the same.  Then look at your model by passing in fit=FALSE to the MARSS() call.  Are you trying to estimate B butyou set Q to zero?  That won't work.  
\n\n
Note if you are estimating D, your error will report problems in A update. If you are estimating C, your error will report problems in U update.  This is because in the MARSS algorithms, the models with D and C are rewritten into a simpler MARSS model with time-varying A and U. If you have set R=0, you might get this error if you are trying to estimate A (or D).\n\n
Did you set a VO (say, diagonal), that is inconsisent with V0T (the covariance matrix implied by the model)?  That can cause problems with the Q update.  Are you estimating C or D, but have rows of c or d that are all zero?  That won't work.  Are you estimating C or D with only one column (one time point) of c or d? Depending on your constraints in C or D that might not work.
\n"
)))

if(number=="convergence")
cat(
  writeLines(
strwrap("MARSS tests both the convergence of the log-likelihood and of the individual 
parameters.  If you just want the log-likelihood, say for model
selection, and don't care too much about the parameter values, then you will be concerned 
mainly that the log-likelihood has converged.  Set abstol to something fairly small
like 0.0001  (in your MARSS call pass in control=list(abstol=0.0001) ).  
Then see if a warning about logLik being converged shows up.  If it doesn't, then you are
probably fine.  The parameters are not at the MLE, but the log-likelihood has 
converged.  This indicates ridges or flat spots in the likelihood.
\n\n
If you are concerned
about getting the MLEs for the parameters and they are showing up as not converged, 
then you'll need to run the algorithm longer (in your MARSS call pass in control=list(maxit=10000) ).
But first think hard about whether you have created a model with ridges and flat 
spots in the likelihood.  Do you have parameters that can create essentially the
same pattern in the data?  Then you may have created a model where the parameters 
are confounded.  Are you trying to fit a model that cannot fit the data?  That
often causes problems.  It's easy to create a MARSS model that is logically inconsistent 
with your data.  Are you trying to estimate both B and U? That is often problematic.  Try demeaning 
your data and setting U to zero.  Are you trying to estimate B and you set tinitx=0? 
tinit=0 is the default, so it is set to this if you did not pass in tinitx in the model list.
You should set tinitx=1 when you are trying to estimate B.
\n"))
)

if(number=="degenvarcov")
  cat( writeLines(strwrap(
"This is not an error but rather an fyi.  Q or R is getting very small.  Because control$allow.degen=TRUE, the code is trying to set Q or R to 0, but in fact, the MLE Q or R is not 0 so setting to 0 is being blocked (correctly).  The code is warning you about this because when Q or R gets really small, you can have numerical problems in the algorithm.  Have you standardized the variances in your data and covariates (explanatory variables)?  In some types of models, that kind of mismatch can drive Q or R towards 0.  This is correct behavior, but you may want to standardize your data so that the variability is on similar scales.\n"
)))

if(number=="x0R0")
cat(writeLines(strwrap(
"Short explanation: This is a constraint imposed by the EM algorithm.  What's happening is 
that x0 cannot be solved 
for because the 0s on the diagonal of R are causing it to disappear from the likelihood.  
Most likely you have set tinitx=1 and your model now has Y_1 = Z x_1. Depending on Z that might not
be solvable for x_1.  If you haven't set R to 0, then pass in allow.degen=FALSE.  
That will stop R being set to 0.  If you did set R to zero, then try setting tinitx=0 
if that makes sense for your model.  You can also try putting a diffuse prior on x0, 
IF you know the implied covariance structure of x0.  However,
if you know that, then setting tinitx=0 is likely ok.
\n\n
Long explanation: If R=0 (or some rows of R) and V10=0 (by setting tinitx=1 and V0=0, 
you have set V10=0), then logically xtT=x10 because when V10=0,
there is no information from the data on x1.  But x1T is the update 
for x10 in the EM algorithm.  This means you cannot estimate x10 since it will 
never change from whatever x10 you start the algorithm with.  Note that the MLE 
value of x10 is not y[1] in this case. In fact, y[1] never enters the
likelihood; you can change y[1] to whatever you want and the likelihood 
doesn't change.  Instead the MLE value of x10 is determined by y[2].  
x11=x10 (for this case) and x21=B x11.  The MLE x10 is then the x10 
that minimizes y[2]-Z x21.
\n\n
If R=0 and you are estimating an initial x (so initial V is set to 0, which is the default), 
you could use method=\"BFGS\".  It will find that MLE x10.  But you probably 
don't want to set tinitx=1.  You would like y[1] to be used not ignored.  
If you use tinitx=0, then your initial x is x00 not x10 and your initial variance
is V00 not V10.  The algorithm will find the x00 that minimized y[1]-Z x10 
where x10=B x00.  This is what you want (most likely).
\n\n
The EM update equations are prone to running into numerical problems when 
R=0 even when tinitx=0 and may exit with an error.  If so, try method=\"BFGS\".
\n\n
If you get this error after the EM algorithm runs for awhile then it probably 
means one of your Q diagonals went to zero and then
the model became illogical.  You can use control=list(alllow.degen=FALSE) to 
stop any Q diagonals being set to zero.  You will likely get a convergence
warning about one of the Qs not having converged (converging to 0 takes an 
infinite number of time steps).
\n\n
If your x0 are independent, you can try a diffuse prior on x0, e.g. V0 diagonal 
equal to say 5, but that is unwise if your model implies that the x0 are 
not independent of each other (or use diffuse=TRUE for a true diffuse prior). 
Don't do that if your model does not imply independent x0 since the prior's 
correlation structure will not match that implied by your model.
\n\n
MARSS does not allow you to specify that x00 or x10 come from the stationary 
distribution of x (assuming it exists), 
i.e. vec(var(x[infinity]))=solve(I-kronecker(B,B))%*%vec(Q) and 
mean(x[infinity])=solve(I-B)%*%u  That's a non-linear constraint. You could get close by 
setting x0, V00 to some value, estimate B and Q, use those to reset x0 and V00, 
re-estimate B and Q, repeat.  But really you should just write a custom LL function
to pass into optim() that specifies that constraint.  So give your custom function 
the B, Q etc, compute x10 and V10 (or x00 and V00) using the stationary 
distribution, and use the kalman filter to give you the log-likelihood.   Note if 
you do something like that, make sure that it is reasonable to assume that x[1]
comes from the stationary distribution (look at your data).  Imposing a initial 
distribution that conflicts violently with your y[1] will lead to highly biased 
parameter estimates
\n\n
Note that these problems occur due to R = 0 (because you set it there or the 
algorithm is going to R = 0 because that is the MLE). Problems with estimating 
x00 (tinitx=0, V0=\"zero\") or x10 (tinitx=1, V0=\"zero\") are much less likely 
to occur when R != 0.\n"
)))
  
if(number=="ts")
  cat( writeLines(strwrap(
"Time series objects have the frequency information embedded (quarter, month, season, etc).  
There are many ways to model quarter, month, season, etc effects and MARSS()
will not guess how you want to model these---there are many, many different places in the 
model that seasonal effects might enter.  You need to pass in the data as a matrix with time going across
the columns (e.g. by using t(as.data.frame.ts(y)) ).  If you want to model quarter, month, etc effects,
you need to use these as covariates in your model.  You can get the frequency information 
as a matrix as with the
command t(as.data.frame.ts(stats:::cycle.ts(y))). Once you have the frequency information (now coded
numeric by the previous command), you can use that in your model to include seasonal effect.  The 
are many different ways in which seasonal effects can modeled (in the x, in the y, in different 
parameters, etc., etc.).
You need to decide how to model them and how to write the model matrices to achieve that.
\n"
)))

if(number=="LLdropped")
cat(writeLines(
    strwrap("The EM algorithm is generally quite robust but it requires inverting the 
variance-covariance matrices and when those inverses
inverses become numerically unstable, the log-likelihood can drop.  The first thing to 
try however is to set 
safe=TRUE in the control list.  This tells MARSS to run the Kalman smoother after each 
parameter update.  This
slows things down, but is a more robust algorithm.  The default is to only run the 
smoother after all parameters
are updated.  If that fails, set maxit (in control) to something smaller than when 
the LL dropping warning starts
and see what is happening to Q and R.  Which one is becoming hard to invert? 
Think long and hard about why this
is happening.  Very likely, something about the way you set up the problem is 
logically forcing this to happen.
\n\n
It may be that you are trying to fit a model that is mathematically inconsistent with 
your data. Are you fitting
a mean-reverting model but the mean implied by the model is different than the mean 
of the data? That won't work.
Are you trying to fit a MARSS model, which w(t) and v(t) errors are a random-walk 
in time and drawn from a multivariate normal, to binned data
where you have multiple time steps at one bin level? Like this, 1,1,1,1,2,2,2,10,10,10,1,1,1,.  
That's not remotely
a random-walk through time.  The binning is not so much the problem.  It's the strings 
of 1's and 2's in a row
that are the problem.  For that kind of binned data, you need some kind of thresholding observation model.
\n\n
If your are fitting models with R=0 or some of the diagonals of R=0, then EM 
can really struggle.  Try BFGS.  If you are fitting
AR-p models with R!=0 and rewritten as a MARSS model, then try using a vague 
prior on x0.  Set x0=matrix(0,1,m)
(or some other appropriate fixed value instead of 0.) and V0=diag(1,m),
where m=number of x's.  That can make it easier to estimate these AR-p with error models.
\n\n
Lastly, try using fun.kf='MARSSkfss' in the MARSS() call.  The tells MARSS() to use the classic Kalman filter/smoother function
rather than Koopman and Durbin's filter/smoother as implemented in the KFAS package.  Normally, 
            the Koopman and Durbin's filter/smoother is more robust but maybe there is something
            about your problem that makes the traditional filter/smoother more robust. Note, they
            normally give identical answers so it would be quite odd to have them different."
            )))

if(number=="kferrors")
  cat( writeLines(strwrap(
"This means the Kalman filter/smoother algorithm became unstable and most likely one of the variances became ill-conditioned.  When that happens the of those matrices are poor, and you will start to get negative values on the diagonals of your variance-covariance matrices.  Once that happens, the inverse of that var-covariance matrix produces an error.  If you get this error, turn on tracing with control$trace=1. This will store the error messages so you can see what is going on.  It may be that you have specified the model in such a way that some of the variances are being forced very close to 0, which makes the var-covariance matrix ill-conditioned. The output from the MARSS call will be the parameter values just before the error occurred.\n"
)))

if(number=="LLunstable")
    cat( writeLines(strwrap(
      "This means, generally, that V0 is very small, say 0, and R is very small and very close to zero.\n"
      )))

  if(number==10)
    cat( writeLines(strwrap(
"Note, this warning is often associated with warnings about the log-likelihood dropping.  The log-likelihood is entering an unstable area, likely a region where it is going steeply to infinity (correctly, probably).
\n\n
The EM algorithm is generally quite robust but it requires inverting the variance-covariance matrices and when those inverses inverses become numerically unstable, the log-likelihood can drop.  The 
first thing to try however is to set safe=TRUE in the control list.  This tells MARSS to run the Kalman smoother after each parameter update.  This slows things down, but is a more robust algorithm.  The default is to only run the smoother after all parameters are updated.  If that fails, set maxit (in control) to something smaller than when the LL dropping warning starts and see what is happening to Q and R.  Which one is becoming hard to invert? Think about why this is happening.  Very likely, something about the way you set up the problem is logically forcing this to happen.
\n\n
It may be that you are trying to fit a model that is mathematically inconsistent with your data. Are you fitting a mean-reverting model but the mean implied by the model is different than the mean of the data? That won't work.  Are you trying to fit a MARSS model, which w(t) and v(t) errors are a random-walk in time and drawn from a multivariate normal, to binned data where you have multiple time steps at one bin level? Like this, 1,1,1,1,2,2,2,10,10,10,1,1,1,.  That's not remotely random-walk through time.  The binning is not so much the problem.  It's the strings of 1's and 2's in a row that are the problem.  For that kind of binned data, you need some kind of thresholding observation model.
\n\n
If your are fitting models with R=0 or some of the diagonals of R=0, then EM can really struggle.  Try BFGS.  If you are fitting AR-p models with R!=0 and rewritten as a MARSS model, then try using a vague prior on x0.  Set x0=matrix(0,1,m) (or some other appropriate fixed value instead of 0.) and V0=diag(1,m), where m=number of x's.  That can make it easier to estimate these AR-p with error models.\n"
)))

if(number=="slowconvergence")
    cat( writeLines(strwrap(
"First thing to do is set silent=2, so you see where MARSS() is taking a long time.  This will
give you an idea of how long each EM iteration is taking so you can estimate how long it will take to get to a
certain number of iterations.  When we get a comment about why the algorithm takes 10,000 iterations to converge, the user is either doing Dynamic
Factor Analysis or they are estimating many variances and they set allow.degen=FALSE.  We'll talk about those two cases.
\n\n
Dynamic Factor Analysis (DFA): Why does this take so long?  By its nature DFA is often a difficult estimation problem because there are two almost equivalent solutions.  The model has a component that looks like this y=z*trend. This is equivalent to y=(z/a)*(a*trend).  That is there exist an an infinite number of trends (a*trend) that will give you the same answer.  However, the likelihood of the (a*trend)'s are not the same since we have a model for the trends---a random walk with variance = 1.  That's pretty flat though for a range of a.  When we have a fairly flat 2D likelihood surface---in this case (z/a)*(a*trend)---EM algorithms take a long time to converge.
\n\n
Variances going to zero: If you set allow.degen=FALSE, and one of your variances is going to zero then it its log is going to negative infinitity and it will take infinite number of iterations to get there (but MARSS() will complain about numerical instability before that).  The log-log convergence test in MARSS is checking for convergence of the log of all the parameters, and clearly the variance going to 0 will not pass this test.  However, you log-likelihood has long converged. So, you want to 'turn off' the convergence test for the parameters and use only the abstol test---which tests if the log-likelihood increased by less than  than some tolerance between iterations.  To do this, pass in a huge value for the slope of the log-log convergence test.  Pass this into your MARSS call: control=list(conv.test.slope.tol=1000)\n"
)))

if(number=="modelobject")
    writeLines(
'Version 3.7 uses model object with attributes while versions 3.4 and earlier did not.  In order, to view 3.4 model fits with MARSS version 3.5+, you need to add on the attributes.  Here is some code to do that.
              
# x is a pre 3.5 marssMLE object from a MARSS call.  x=MARSS(....)
x$marss=x$model
class(x$marss)="marssMODEL"
attr(x$marss,"form")="marss"
attr(x$marss,"par.names")=names(x$marss$fixed)
tmp=x$marss$model.dims
for(i in names(tmp)) if(length(tmp[[i]])==2) tmp[[i]]=cbind(tmp[[i]],1)
tmp$x=c(sqrt(dim(kemz$model$fixed$Q)[1]), dim(kemz$model$data)[2])
tmp$y=c(sqrt(dim(kemz$model$fixed$R)[1]), dim(kemz$model$data)[2])
attr(x$marss,"model.dims")=tmp
attr(x$marss,"X.names")=x$marss$X.names
              
class(x$model)="marssMODEL"
x$model=x$form.info$marxss
attr(x$model,"form")=x$form.info$form
attr(x$model,"par.names")=names(x$model$fixed)
tmp=x$form.info$model.dims
for(i in names(tmp)) if(length(tmp[[i]])==2) tmp[[i]]=cbind(tmp[[i]],1)
tmp$x=c(sqrt(dim(kemz$model$fixed$Q)[1]), dim(kemz$model$data)[2])
tmp$y=c(sqrt(dim(kemz$model$fixed$R)[1]), dim(kemz$model$data)[2])
attr(x$model,"model.dims")=tmp
attr(x$model,"X.names")=x$marss$X.names
              
#now this should work
coef(x, type="matrix")
')

if(number==23)
    cat( writeLines(strwrap(
"This is the same error as number 22 except that the 0s on the diagonal of R are arising because allow.degen=TRUE (this is the default setting in the control list) and R is getting very small.  MARSS attempts to set R to 0, but the constraint that x0 associated with R=0 comes into play.  MARSS then blocks the setting of R to 0 and warns you.  You can set allow.degen=FALSE, but it is just an informational warning.  There is nothing wrong per se.\n"
)))

if(number=="residvarinv")
  cat( writeLines(strwrap(
"The computation of the standardized residuals requires taking the Cholesky decomposition of the joint variance-covariance matrix of the observation and state residuals.  This is matrix is not invertible for some reason.  If you have missing data and a non-diagonal Q or R, try Harvey=FALSE (if you set Harvey=TRUE).  This will compute the exact joint variance-covariance matrix.  However, in some cases, the exact matrix is also not invertible.  This could occur is say Q is non-diagonal and all the data are missing in the last time-step.  When the matrix is not invertible, the standardized residuals for that time-step are set to NAs.\n"
)))

if(number=="varcovstruc")
  cat( writeLines(strwrap(
    "The structure of variance-covariance matrices have many constraints. These constraints arise
     due to the nature of variance-covariance matrices and their estimation, not due to MARSS per se.
     They must be symmetric.  If numeric (meaning no estimated values), they must be positive definite. 
     You cannot fix the covariances and estimate the variances; or vis-a-versa.  There are many 
     constraints on shared values. You cannot have shared estimated values across variances and covariances.
     Within a block with a shared variance, the covariances must be equal (or 0).  Across pairs of 
     different shared variances, the covariances must all be equal.  So if there are 3 variances of 'a' and 
     another of 'b', the covariances between all the 'a' and 'b' must be the same.  The covariances cannot
     be shared across different pairs.  So if there is another variance of 'c', the covariance between it and
     the 'a' and 'b' must be different.  If there are blocks within the matrix (so it is a block diagonal
     matrix), there can be no shared values across blocks unless the blocks are identical.
     If you set method=BFGS, however, there are extra constraints.  In this case, the code
     requires that the matrices be diagonal, unconstrained, or equalvarcov.  This is due to the fact
     that the code use a Cholesky decomposition to ensure that the matrices stay positive definite
     during the estimation iterations.\n"
  )))

  if(number=="HessianNA")
    cat( writeLines(strwrap(
      "The variance-covariance matrix can be estimated (large sample estimator) from
      the inverse of the Hessian of the log-likelihood function at the MLE parameter
      values.  The Hessian is the second partial derivative of a matrix function.
      The Hessian
      of the log-likelihood function at the MLEs is the observed Fisher information.
      The observed Fisher information is an estimator of large-sample 
      variance-covariance
      matrix of the estimated parameters.  The MARSS package provides 3 ways to 
      compute the Hessian: the recursive algorithm by Harvey (1989), a numerical 
      estimate using the dfHess() function from the nmle package, and a numerical
      estimate from the optim() function.  The calculation of the Hessian 
      associated with the variance terms (Q & R) is prone to numerical errors.
      When this happens, an NA is put on the diagonal of the Hessian for that
      parameter value. No standard errors or CIs can be computed for that value.
      A Hessian with many NAs is probably a sign that you have a poor model
      (meaning your model is not a good description of the data) or you do not
      have enough data given the complexity of your model.\n"
  )))

}
