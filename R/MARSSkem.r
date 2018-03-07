#######################################################################################################
#   KEM function
#   Minimal error checking is done.  You should run is.marssMLE(MLEobj) before calling this.
#   Maximization using an EM algorithm with Kalman filter
#######################################################################################################
MARSSkem = function( MLEobj ) {
  MODELobj=MLEobj[["marss"]]
  # This is a core function and does not check if user specified a legal or solveable model. 
  # y is MLEobj$marss$data with the missing values replaced by 0
  kf.x0 = ifelse(MODELobj[["tinitx"]]==1,"x10","x00")  #the initial conditions treatment "x00" x0 is at t=0 or "x01" x0 is at t=1
  #kf.x0=x00 prior is defined as being E[x(t=0)|y(t=0)]; xtt[0]=x0; Vtt[0]=V0
  #kf.x1=x10 prior is defined as being E[x(t=0)|y(t=0)]; xtt1[1]=x0; Vtt1[1]=V0
  
  #The model will be form = marss, so use base function for that form here
  constr.type = describe_marss( MODELobj )
  #Check that model is allowed given the EM algorithm constaints; returns some info on the model structure
  if(MLEobj[["control"]][["trace"]] != -1){
    errhead = "\nErrors were caught in MARSSkemcheck \n"
    errmsg = " Try using foo=MARSS(..., fit=FALSE), then summary(foo$model) to see what model you are trying to fit.\n" 
    tmp=MARSSkemcheck( MLEobj )
    if(!tmp$ok){ cat(c(errhead, tmp$msg, errmsg)); stop("Stopped in MARSSkemcheck() due to specification problem(s).\n", call.=FALSE) }
  }
  
  #set up holders for warning messages
  msg=NULL; stop.msg=NULL; msg.kem=NULL; msg.kf=NULL; msg.conv=NULL #error messages
  
  ## attach would be risky here since user might have one of these variables in their workspace    
  y = MODELobj[["data"]]#must have time going across columns
  d = MODELobj[["free"]]      # D or free matrix
  f = MODELobj[["fixed"]]   # f matrix
  inits = MLEobj[["start"]]
  model.el = attr(MODELobj, "par.names")
  model.dims = attr(MODELobj, "model.dims")
  n =  model.dims[["data"]][1]; TT = model.dims[["data"]][2]; m = model.dims[["x"]][1]
  Id = list(m = diag(1,m), n = diag(1,n)); IIm=diag(1,m) # identity matrices
  
  control = MLEobj$control
  stopped.with.errors=FALSE; kf=NULL; condition.limit=1E10
  
  ## Set up MLE object for the iterations
  MLEobj.iter = MLEobj
  MLEobj.iter$constr.type=constr.type
  MLEobj.iter$par=list()
  ## assign the starting parameter values; use fixed values where fixed otherwise use inits
  for(elem in model.el) MLEobj.iter$par[[elem]]=inits[[elem]]
  
  ## make a list of time-varying and fixed parameters
  time.varying = fixed = list()
  for(elem in model.el) {
    if( is.fixed(d[[elem]]) ){ 
      MLEobj.iter$par[[elem]]=matrix(0,0,1)
      fixed[[elem]] = TRUE
    }else{ fixed[[elem]] = FALSE }
    if( model.dims[[elem]][3] == 1){
      time.varying[[elem]] = FALSE
    }else{ time.varying[[elem]] = TRUE }  #time-varying
  }
  #flags for whether a 0 was set on R or Q diagonals; determines whether various is.zero diagonal matrices are recomputed
  set.degen=list(Q=FALSE, R=FALSE, V0=FALSE)
  #define a couple constants that come up a lot
  #this is L%*V0%*%t(L)
  tmpV=tcrossprod(parmat(MLEobj.iter, "L", t=1)$L %*% parmat(MLEobj.iter, "V0", t=1)$V0, parmat(MLEobj.iter, "L", t=1)$L)
  IIzV0=diag(as.numeric(diag(tmpV)==0),m)
  IImIIzV0 = (IIm-IIzV0)
  
  ## zero out the missing values
  y[is.na(y)]=0
  
  ## Set up variable for debuging and diagnostics
  iter.record=list(par=NULL,logLik=NULL)  
  
  ################# The main EM loop which will run until tol reached or max.iter reached
  #######################################################################################   
  
  #set up the convergence flags
  cvg = 1 + control$abstol
  MLEobj.iter$logLik = NA #start with no value
  # 72 means no info yet; 0 means converged
  MLEobj.iter$conv.test=list(convergence=72, messages="No convergence testing performed.\n", not.converged.params=names(coef(MLEobj.iter,type="vector")), converged.params=c() )
  
  if(control$silent==2) cat("EM iteration: ")
  for(iter in 1:(control$maxit+1)) { #+1 so that the iter.record and kf are run for the last EM iteration
    if(control$silent==2) cat(iter," ")
    ################# E STEP Estimate states given U,Q,A,R,B,X0 via Kalman filter
    #####################################################################################
    kf.last = kf
    kf = MARSSkf( MLEobj.iter )  #kf selects the function based on MLEobj$fun.kf
    if(!kf$ok) { 
      if(control$trace>0){ msg.kf=c(msg.kf,paste("iter=",iter," ",kf$errors) )
      }else msg.kf=kf$errors
      stop.msg = paste("Stopped at iter=",iter," in MARSSkem() because numerical errors were generated in the Kalman filter.\n",sep="")
      stopped.with.errors=TRUE; break
    }
    MLEobj.iter$kf=kf
    MLEobj.iter$logLik=kf$logLik
    if(control$demean.states) {
      xbar = try(apply(cbind(kf$x0T,kf$xtT),1,mean) )
      MLEobj.iter$kf$xtT = kf$xtT-xbar
      MLEobj.iter$kf$x0T = kf$x0T-xbar
    }
    
    Ey = MARSShatyt( MLEobj.iter )
    if(!Ey$ok) { 
      if(control$trace>0){ msg.kf=c(msg.kf,paste("iter=",iter," ",Ey$errors) )
      }else msg.kf=Ey$errors
      stop.msg = paste("Stopped at iter=",iter," in MARSSkem() because numerical errors were generated in MARSShatyt\n",sep="")
      stopped.with.errors=TRUE; break
    }
    MLEobj.iter$Ey = Ey
    
    # This is a diagnostic line that checks if the solution is becoming unstable; likelike.old is set first at iter=1
    if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(MLEobj.iter$logLik) == TRUE ) cvg = MLEobj.iter$logLik - loglike.old  
    if(iter > 2 & cvg < -sqrt(.Machine$double.eps)) {
      if(control$trace>0){ 
        msg.kem=c(msg.kem,paste("iter=",iter," LogLike DROPPED.  old=", loglike.old, " new=", MLEobj.iter$logLik, "\n", sep=""))
      }else msg.kem = "MARSSkem: The soln became unstable and logLik DROPPED.\n"
    }
    
    ################
    # Keep a record of the iterations for debugging and convergence diagnostics
    ################################################################
    if(control$trace>0){ # if trace is on, keep the full record over all iterations
      iter.record$par=rbind(iter.record$par,coef(MLEobj.iter,type="vector"))
      iter.record$logLik=c(iter.record$logLik,MLEobj.iter$logLik)
      if(!is.null(MLEobj.iter[["kf"]][["errors"]])) {
        msg.kf=c(msg.kf, paste("iter=",iter," ", kf$errors, sep=""))
      }
      MLEobj.iter$iter.record=iter.record
    }else{ #Otherwise keep just last (control$conv.test.deltaT+1) iterations for diagnostics
      iter.record$par=rbind(iter.record$par,coef(MLEobj.iter,type="vector"))
      iter.record$logLik=c(iter.record$logLik,MLEobj.iter$logLik)
      tmp.len=dim(iter.record$par)[1]
      if(tmp.len>(control$conv.test.deltaT+1)) {
        iter.record$par=as.matrix(iter.record$par[(tmp.len-control$conv.test.deltaT):tmp.len,,drop=FALSE])
        iter.record$logLik = iter.record$logLik[(tmp.len-control$conv.test.deltaT):tmp.len]
      }
      MLEobj.iter$iter.record=iter.record
    }
    
    ################
    # Convergence Test
    ################################################################
    if(iter >= control$minit){  # then do convergence testing
      if( cvg >= 0 && cvg < control$abstol ){
        if(iter>=control$min.iter.conv.test){ 
          MLEobj.iter$conv.test = loglog.conv.test(iter.record, iter, deltaT=control$conv.test.deltaT, tol=control$conv.test.slope.tol)
          if(MLEobj.iter$conv.test$convergence!=1) break  #1=not converged, keep going; 0=converged; anything else=problem
        }else MLEobj.iter$conv.test$convergence=3  #abstol reached log-log hasn't run yet because min.iter.cov.test not reached
      }else MLEobj.iter$conv.test$convergence=4 #minit reached but not abstol
    }
    if(iter>control$maxit){ iter=control$maxit; break } #reset iter to maxit since needed to determine if stopped due to reaching maxit
    
    # Store loglike for comparison to new one after parameters are updated
    loglike.old = MLEobj.iter$logLik
    
    #set the parameters at t=1
    par1=parmat(MLEobj.iter,t=1)
    
    ################# M STEP update U,Q,A,R,B,X0 via ML given x(t) estimate
    # Update Q and R
    # Run Kalman smoother again to update the hidden states expectations
    # Update the other parameters
    
    ################################################################
    # Get new R subject to its constraints
    ################################################################
    # For the degen test, I require that d is a design matrix;
    if(control[["allow.degen"]]){
      tmp=degen.test("R", MLEobj.iter, iter) #this will test degeneracy and replace diags with 0s if needed
      MLEobj.iter=tmp$MLEobj; msg.kem=c(msg.kem, tmp$msg)
      #update d and f because some of the R diagonals may have been set to 0
      if(tmp$set.degen){ #then some diagonals set to 0 so need to update values
        d$R=MLEobj.iter$marss$free$R
        f$R=MLEobj.iter$marss$fixed$R
        kf=MLEobj.iter$kf
        Ey=MLEobj.iter$Ey
        fixed$R = is.fixed(d$R)
        set.degen$R=TRUE  #flag so that the identity matrices are redone
      }
    }
    
    #Now run the standard EM update equations
    if(!fixed[["R"]]){
      sum1 = t.dR.dR = 0
      Z=par1$Z
      A=par1$A
      for (i in 1:TT) {
        if(time.varying[["Z"]] & i>1){ Z = parmat(MLEobj.iter, "Z", t=i)$Z }
        if(time.varying[["A"]] & i>1){ A = parmat(MLEobj.iter, "A", t=i)$A }
        if(time.varying[["R"]] | i==1){ 
          dR = sub3D(d[["R"]],t=i) #by def, i goes to TT if time-varying
          t.dR.dR = t.dR.dR + crossprod(dR)
        }
        hatyt = Ey[["ytT"]][,i,drop=FALSE]; hatyxt=sub3D(Ey[["yxtT"]],t=i); hatOt = sub3D(Ey[["OtT"]],t=i)
        hatPt = kf[["VtT"]][,,i]+tcrossprod(kf[["xtT"]][,i,drop=FALSE])
        hatxt = kf[["xtT"]][,i,drop=FALSE]
        sum1a = (hatOt - tcrossprod(hatyxt, Z) - tcrossprod(Z, hatyxt)- tcrossprod(hatyt, A) - tcrossprod(A, hatyt) + 
                   tcrossprod(Z%*%hatPt, Z) + tcrossprod(Z%*%hatxt, A) + tcrossprod(A, Z%*%hatxt) + tcrossprod(A)) #A%*%t.A
        sum1a = symm(sum1a) #enforce symmetry function from MARSSkf
        sum1 = sum1 + crossprod(dR, vec(sum1a))
      }
      if(time.varying[["R"]]){ 
        if(length(t.dR.dR)==1){ inv.dR=1/t.dR.dR }else{ inv.dR = pcholinv(t.dR.dR) }
      }else{ 
        if(length(t.dR.dR)==1){ inv.dR=(1/t.dR.dR)/TT }else{ inv.dR = pcholinv(t.dR.dR)/TT }
      }       
      
      MLEobj.iter[["par"]][["R"]] = inv.dR%*%sum1 #.03
      par1[["R"]]=parmat(MLEobj.iter,"R",t=1)$R
      
      #Start~~~~~~~~Error checking
      R=par1$R  #reset
      for(i in 1:model.dims[["R"]][3]){
        if(time.varying$R & i>1) R=parmat(MLEobj.iter,"R",t=i)$R  
        if(any(eigen(R,symmetric=TRUE,only.values=TRUE)$values<0)) {
          stop.msg=paste("Stopped at iter=",iter," in MARSSkem: solution became unstable. R update is not positive definite.\n",sep="")
          stopped.with.errors=TRUE;  
          break }
      }
      if(stopped.with.errors) break      
      if( control$safe && !fixed[["R"]] ){  
        new.kf=rerun.kf("R", MLEobj.iter, iter)
        if(!new.kf$ok){
          stopped.with.errors=TRUE
          msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
          break
        }else{
          kf=new.kf$kf
          MLEobj.iter$kf = kf
          MLEobj.iter$logLik=kf$logLik
          Ey = MARSShatyt( MLEobj.iter )
          MLEobj.iter$Ey = Ey
          msg.kem=c(msg.kem, new.kf$msg.kem)
        }
      }
      
    } #R not fixed
    
    ################
    # Get new Q subject to its constraints
    ################################################################
    #Start the testing for 0s along the diagonal of Q
    # For the degen test, I require that d is a design matrix;
    if(control[["allow.degen"]]){
      tmp=degen.test("Q", MLEobj.iter, iter) #this will test degeneracy and replace diags with 0s if needed
      MLEobj.iter=tmp$MLEobj; msg.kem=c(msg.kem, tmp$msg)
      if(tmp$set.degen){
        d$Q=MLEobj.iter$marss$free$Q
        f$Q=MLEobj.iter$marss$fixed$Q
        kf=MLEobj.iter$kf
        Ey=MLEobj.iter$Ey
        fixed$Q = is.fixed(d$Q)
        set.degen$Q=TRUE  #flag so that the identity matrices are redone
      }
    }
    
    #Then do the regular EM update
    if( !fixed[["Q"]] ){ #dim d$Q =0 or d$Q all zeros
      # If you treat x0 as at t=1 then 
      # S00 = 0; S11 = 0; S10 = 0; X1 = 0; X0 = 0; TT.numer = TT-1
      # Otherwise if x0 is at t=0 follow Shumway and Stoffer (S&S2006 Eqn 6.67-69)
      
      dQ = sub3D(d$Q,t=1) #won't be all zeros due to !is.fixed
      B = par1$B
      U = par1$U
      if(kf.x0=="x00"){
        TT.numer = TT; 
        # IImIIzV0 = (IIm-IIz$V0[,,1]); IIzV0 = IIz$V0[,,1]      
        X0 =  IImIIzV0%*%kf$x0T + IIzV0%*%par1$x0
        S00 = kf$V0T + tcrossprod(X0)
        S10 = kf$Vtt1T[,,1] + tcrossprod(kf$xtT[,1,drop=FALSE], X0);  #where diag.V0=0 kf$Vtt1T[,,1]=0 since x at t-1 is fixed
        S11 = kf$VtT[,,1] + tcrossprod(kf$xtT[,1,drop=FALSE])
        X1 = kf$xtT[,1,drop=FALSE]
        sum1a = (S11 - tcrossprod(B, S10) - tcrossprod(S10, B) + tcrossprod(B%*%S00, B)
                 - tcrossprod(U, X1) - tcrossprod(X1, U) + tcrossprod(U, B%*%X0) + tcrossprod(B%*%X0, U) + tcrossprod(U)) #U%*%t.U
        sum1a = symm(sum1a) #symmetry function from MARSSkf
        sum1 = crossprod(dQ, vec(sum1a))
        t.dQ.dQ = crossprod(dQ)
      }
      if(kf.x0=="x10"){
        sum1 = 0; t.dQ.dQ=0; TT.numer = TT-1
        #Gharamani treatment of initial condition; the initial condition specifies x at t=1
      }
      
      for (i in 2:TT) {
        S00 = kf[["VtT"]][,,i-1] + tcrossprod(kf[["xtT"]][,i-1,drop=FALSE]) #sum 2:T E(xt1T%*%t(xt1T))
        S10 = kf[["Vtt1T"]][,,i] + tcrossprod(kf[["xtT"]][,i,drop=FALSE],kf[["xtT"]][,i-1,drop=FALSE])  #sum 2:T E(xtT%*%t(xt1T))
        S11 = kf[["VtT"]][,,i] + tcrossprod(kf[["xtT"]][,i,drop=FALSE])   #sum 2:T E(xtT%*%t(xt1T))
        X0 = kf[["xtT"]][,i-1,drop=FALSE]      #sum 2:T E(xt1T)
        X1 = kf[["xtT"]][,i,drop=FALSE]        #sum 2:T E(xtT)
        if(time.varying[["B"]] ){ B = parmat(MLEobj.iter, "B", t=i)$B }
        if(time.varying[["U"]] ){ U = parmat(MLEobj.iter, "U", t=i)$U }
        if(time.varying[["Q"]] ){ 
          dQ = sub3D(d[["Q"]],t=i)
          t.dQ.dQ = t.dQ.dQ + crossprod(dQ)
        }
        sum1a = (S11 - tcrossprod(B,S10) - tcrossprod(S10, B) + tcrossprod(B%*%S00,B)
                 - tcrossprod(U, X1) - tcrossprod(X1, U) + tcrossprod(U,B%*%X0) + tcrossprod(B%*%X0, U) + tcrossprod(U))
        sum1a = symm(sum1a) #enforce symmetry function from MARSSkf
        sum1 = sum1 + crossprod(dQ, vec(sum1a)) 
      }
      
      #pcholinv because there might be all zero cols in dQ; won't equal matrix(0,1,1) since !is.fixed
      if(time.varying[["Q"]]){ 
        if(length(t.dQ.dQ)==1){ inv.dQ=1/t.dQ.dQ }else{ inv.dQ = pcholinv(t.dQ.dQ) }
      }else{
        t.dQ.dQ = crossprod(dQ) #t.dQ%*%dQ
        if(length(t.dQ.dQ)==1){ inv.dQ=(1/t.dQ.dQ)/TT.numer }else{ inv.dQ = pcholinv(t.dQ.dQ)/TT.numer }
      }
      #0 will appear in par where there are all 0 cols in d since inv.dQ will be 0 row/col there      
      MLEobj.iter$par$Q=inv.dQ%*%sum1
      par1$Q=parmat(MLEobj.iter,"Q",t=1)$Q
      
      #Start~~~~~~~~~~~~Error checking
      Q=par1$Q #reset
      for(i in 1:max(dim(d$Q)[3],dim(f$Q)[3])){
        if(time.varying$Q & i>1) Q=parmat(MLEobj.iter, "Q", t=i)$Q  
        if(any(eigen(Q,symmetric=TRUE,only.values=TRUE)$values<0)) {
          stop.msg=paste("Stopped at iter=",iter," in MARSSkem: solution became unstable. Q update is not positive definite.\n",sep="")
          stopped.with.errors=TRUE; break 
        }
      }
      if(stopped.with.errors) break      
      if( control$safe &&  !fixed[["Q"]] ){
        new.kf=rerun.kf("Q", MLEobj.iter, iter)
        if(!new.kf$ok){
          stopped.with.errors=TRUE
          msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
          break
        }else{
          kf=new.kf$kf
          MLEobj.iter$kf=kf
          MLEobj.iter$logLik=kf$logLik
          Ey = MARSShatyt( MLEobj.iter )
          MLEobj.iter$Ey=Ey
          msg.kem=c(msg.kem, new.kf$msg.kem)
        }
      }
    } #if Q not fixed
    
    #This code sets up the IIz and star (inverse) matrices needed
    #This only needs to be done at iter=1 or if Q, R or V0 might have changed
    if( !fixed[["Q"]] | !fixed[["R"]] | !fixed[["V0"]] | set.degen[["Q"]] | set.degen[["R"]] | iter==1 ){
      #Set up the variance matrices needed for the degenerate case
      if(iter==1){ 
        IIz=star=list()  #IIz location of 0s on diagonal
        elems=c("Q", "R", "V0") 
        for(elem in elems){ #set up the arrays; only needed at iter=1
          #var-cov = G Q t(G) and H R t(H)
          elem1="L"
          if(elem=="Q") elem1 = "G"
          if(elem=="R") elem1 = "H"
          star[[elem]]=IIz[[elem]]=array(0,dim=c(model.dims[[elem1]][1],model.dims[[elem1]][1],max(model.dims[[elem]][3],model.dims[[elem1]][3]) ))
        }
      }else{ 
        elems=c("Q", "R", "V0")[c((!fixed[["Q"]] | set.degen[["Q"]]), (!fixed[["R"]] | set.degen[["R"]]), !fixed[["V0"]])] 
      }
      for( elem in elems ){
        elem1="L"
        if(elem=="Q") elem1 = "G"
        if(elem=="R") elem1 = "H"
        thedim=model.dims[[elem1]][1] 
        maxT=model.dims[[elem]][3] #Q, R, or V0
        maxT1=model.dims[[elem1]][3] #G, H, or L
        #star$elem is mathbb{elem}; the bolded Q, R, and L in section 4.4 and 5
        #move up so only done once; star[[elem]]=IIz[[elem]]=array(0,dim=c(thedim,thedim,maxT))
        if(maxT==1)  QRV = parmat(MLEobj.iter,elem,t=1)[[elem]]
        if(maxT1==1) GHL = parmat(MLEobj.iter,elem1,t=1)[[elem1]]
        for(i in 1:max(maxT,maxT1)){
          if(maxT!=1)   QRV = parmat(MLEobj.iter,elem,t=i)[[elem]]
          if(maxT1!=1)  GHL = parmat(MLEobj.iter,elem1,t=i)[[elem1]]
          #this is being done to find the zeros on the diagonal
          pari = tcrossprod(GHL %*% QRV, GHL)
          #These are the identity matrices used to identify the location of deterministic rows of x;
          #section 7.2 in EM Derivations "Idntifying the fully deterministic x rows"
          #IIz means location of 0s on diagonal of var-cov matrix GHL%*%QRV%*%t(GHL)
          if(set.degen[[elem]] | iter==1){
            IIz[[elem]][,,i]  = diag(as.numeric(diag(pari)==0),thedim)
            #I the locations of 0s on diagonal of Q are time-constant; see section 7.2
            if(elem=="Q" & max(maxT, maxT1) !=1){
              if(!all.equal(IIz[[elem]][,,1],IIz[[elem]][,,i])){ 
                stop.msg = paste("Stopped at iter=",iter," in MARSSkem. IIz$Q (location of 0s on diagonal) must be time invariant.\nYou probably want to set allow.degen=FALSE if it is true.\n",sep="")
                stopped.with.errors=TRUE; break }
            }
          }
          #this is defining the bolded Q, R and K in section 4.4 and 5 of the EM Derivation
          starmultiplier = tcrossprod(pcholinv(crossprod(GHL)),GHL)
          star[[elem]][,,i] = crossprod(starmultiplier, pcholinv(QRV)%*%starmultiplier)
        }
        set.degen[[elem]]=FALSE #reset so this code is not run again
        if(elem=="V0" & (iter==1 | set.degen$V0)){ #then 0 location in V0 has potentially changed (via degen.test(V0))
          #set up the diag matrices needed often
          IIzV0=sub3D(IIz$V0,t=1)
          IImIIzV0 = (IIm-IIzV0)
        }      
      }  #for over elems
      if(stopped.with.errors) break
    } #if Q, R or V0 is not fixed or a value was set to 0 via allow.degen or iter=1
    
    ################################################################
    # Get new x0 subject to its constraints
    ################################################################
    if( !fixed[["x0"]] ){  # some element needs estimating
      f.x0=sub3D(f$x0,t=1) 
      d.x0=sub3D(d$x0,t=1) 
      A=par1$A; Z=par1$Z; B=par1$B; U=par1$U
      Qinv = sub3D(star$Q,t=1); diag.Q=1-takediag(IIz$Q[,,1]) 
      Rinv = sub3D(star$R,t=1); diag.R1=1-takediag(IIz$R[,,1]); 
      IIz.R = sub3D(IIz$R,t=1) 
      x0.degen.update = FALSE
      diag.V0=1-takediag(IIzV0)
      nQ0=sum(diag.Q==0)
      if(any(diag.Q==0)) x0.degen.update=!is.fixed(d.x0[diag.Q==0,,drop=FALSE])       
      numer=denom=0
      if(kf.x0=="x00") hatxt0=kf$x0T else hatxt0=kf$xtT[,1,drop=FALSE]
      if(any(diag.V0==1)){ #meaning some V0 positive
        denom = crossprod( d.x0, star$V0[,,1]%*%d.x0) 
        numer = crossprod( d.x0, hatxt0-f.x0) 
      }
      
      if(any(diag.V0==0)){
        AdjM=B; AdjM[AdjM!=0]=1
        if(kf.x0=="x00"){
          #set up values for t=1
          Bstar=B; Bstar.tm=IIm
          Ustar=U; Ustar.tm=0*U
          Mt=AdjM
          IId.tm=IIm #t=0; no w's
          IId=makediag(1-diag.Q) #only can be w free if Q==0
          if(any(diag.Q==0)&any(diag.Q!=0)) #which did not pick up a zero from B
            IId[diag.Q==0,diag.Q==0][1 + 0:(nQ0 - 1)*(nQ0 + 1)]=apply(Mt[diag.Q==0,diag.Q!=0,drop=FALSE]==0,1,all)
          Delta7 = kf$xtT[,1,drop=FALSE]-B%*%(IImIIzV0%*%hatxt0 + IIzV0%*%f.x0)-U
          Delta8 = B%*%IIzV0%*%d.x0 #since IId.tm and Bstar.tm are IIm
          numer = numer+crossprod(Delta8, Qinv%*%Delta7)
          denom = denom+crossprod(Delta8, Qinv%*%Delta8)
        }else{  #x10
          Bstar=IIm; Ustar=0*U
          IId.tm=0; IId=IIm
          Mt=IIm
        }
        if(any(IId==1)){ #again t=1
          Delta5=Ey$ytT[,1,drop=FALSE]-Z%*%((IIm-IId)%*%kf$xtT[,1,drop=FALSE])-Z%*%IId%*%(Bstar%*%(IImIIzV0%*%hatxt0 + IIzV0%*%f.x0)+Ustar)-A
          Delta6=Z%*%IId%*%Bstar%*%(IIzV0%*%d.x0)         
          if(any(diag.R1==0)){
            if(any(crossprod(Delta6,IIz.R)%*%Delta6 !=0)){
              stop.msg = paste("Stopped at iter=",iter," in MARSSkem at x0 update.\n There are 0s on R diagonal. x0 assoc with these must be fixed (not estimated)\n when using the EM algorithm. Try method=\"BFGS\".  Type MARSSinfo(\"x0R0\") for help.\n", sep="")
              stopped.with.errors=TRUE
              break
            }
          }
          numer = numer + crossprod(Delta6, Rinv%*%Delta5) 
          denom = denom + crossprod(Delta6, Rinv%*%Delta6) 
        }
        #t>1 will break out as soon as no IId=1
        for(t in 2:TT){ #start at 2 if x00; at 3 if x10
          if(!any(IId==1)) break
          if(time.varying$A) A = parmat(MLEobj.iter,c("A"),t=t)$A
          if(time.varying$B) B = parmat(MLEobj.iter,c("B"),t=t)$B
          if(time.varying$U) U = parmat(MLEobj.iter,c("U"),t=t)$U
          if(time.varying$Z) Z = parmat(MLEobj.iter,c("Z"),t=t)$Z
          if(time.varying$R){ Rinv = star$R[,,t]; IIz.R=sub3D(IIz$R,t=t) }
          if(time.varying$Q) Qinv = star$Q[,,t]
          Ustar.tm=Ustar
          Ustar=B%*%Ustar + U
          Bstar.tm=Bstar
          Bstar=B%*%Bstar
          if(t<=(m+1)){
            IId.tm=IId
            Mt=AdjM%*%Mt
            IId=makediag(1-diag.Q) #only can be w free if Q==0
            if(any(diag.Q==0)&any(diag.Q!=0)) 
              IId[diag.Q==0,diag.Q==0][1 + 0:(nQ0 - 1)*(nQ0 + 1)]=apply(Mt[diag.Q==0,diag.Q!=0,drop=FALSE]==0,1,all)
          }
          if(any(IId==1)){
            Delta5=Ey$ytT[,t,drop=FALSE]-Z%*%((IIm-IId)%*%kf$xtT[,t,drop=FALSE])-Z%*%IId%*%(Bstar%*%(IImIIzV0%*%hatxt0 + IIzV0%*%f.x0)+Ustar)-A
            Delta6=Z%*%IId%*%Bstar%*%(IIzV0%*%d.x0)
            #Deal with Delta6=0 and Rinv=Inf, so 0*Inf
            if(any(diag.R1==0)){
              if(any(crossprod(Delta6, IIz.R)%*%Delta6 !=0)){
                stop.msg = paste("Stopped at iter=",iter," in MARSSkem at x0 update.\n There are 0s on R diagonal. x0 assoc with these must be fixed (not estimated)\n when using the EM algorithm. Try method=\"BFGS\".  Type MARSSinfo(\"x0R0\") for help.\n", sep="")
                stopped.with.errors=TRUE
                break
              }
            }
            numer = numer + crossprod(Delta6, Rinv%*%Delta5)
            denom = denom + crossprod(Delta6, Rinv%*%Delta6)
          }
          if(any(IId.tm==1)){
            Delta7 = kf$xtT[,t,drop=FALSE]-B%*%(IIm-IId.tm)%*%kf$xtT[,t-1,drop=FALSE]-B%*%IId.tm%*%(Bstar.tm%*%(IImIIzV0%*%hatxt0 + IIzV0%*%f.x0)+Ustar.tm)-U
            Delta8 = B%*%IId.tm%*%Bstar.tm%*%(IIzV0%*%d.x0)
            numer = numer + crossprod(Delta8, Qinv%*%Delta7) 
            denom = denom + crossprod(Delta8, Qinv%*%Delta8) 
          }
        } #for t
      } #any diag.LAM=0
      if(length(denom)==1){ denom = 1/denom }else{ denom=try(pcholinv( denom ) ) }
      if(inherits(denom, "try-error") | (length(denom)==1 && denom==Inf)){
        stop.msg = paste("Stopped at iter=",iter," in MARSSkem at x0 update. denom is not invertible. \n This means that some of the x0 cannot be estimated. Type MARSSinfo('denominv') for more info. \n", sep="")
        stopped.with.errors=TRUE
        break
      }
      MLEobj.iter$par$x0 = denom%*%numer
      if( !is.matrix(MLEobj.iter$par$x0) ) MLEobj.iter$par$x0=matrix(MLEobj.iter$par$x0,dim(d$x0)[2],1)
      par1$x0=parmat(MLEobj.iter,"x0",t=1)$x0
      
      #~~~~~~~~Error checking  
      if( control$safe &&  !fixed[["x0"]] ){
        new.kf=rerun.kf("x0", MLEobj.iter, iter)
        if(!new.kf$ok){
          stopped.with.errors=TRUE
          msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
          break
        }else{
          kf=new.kf$kf
          MLEobj.iter$kf=kf
          MLEobj.iter$logLik=kf$logLik
          Ey = MARSShatyt( MLEobj.iter )
          MLEobj.iter$Ey=Ey
          msg.kem=c(msg.kem, new.kf$msg.kem)
        }
      } 
    } #x0 is not fixed
    
    ################
    # Get new V0 subject to its constraints
    ################################################################
    if( !fixed[["V0"]] ){  # some element needs estimating (obviously V0!=0)
      dV0=sub3D(d$V0,t=1)
      V0.update = chol2inv(chol(crossprod(dV0)))%*%crossprod(dV0, vec(kf$V0T))
      MLEobj.iter$par$V0 = V0.update
      if(!is.matrix(MLEobj.iter$par$V0) ) MLEobj.iter$par$V0=matrix(MLEobj.iter$par$V0,dim(d$V0)[2],1)
      par1$V0=parmat(MLEobj.iter,"V0",t=1)$V0
      
      #~~~~~~~~Error checking
      if(any(eigen(par1$V0,symmetric=TRUE,only.values=TRUE)$values<0)) {
        tmp=""
        if(any(abs(Re(eigen(par1$B,only.values=TRUE)$values)>1))) tmp="Your B matrix is outside the unit circle.  This is likely the problem.\n"
        stop.msg=paste("Stopped at iter=",iter," in MARSSkem: solution became unstable. V0 update is not positive definite.\n",tmp,sep="")
        stopped.with.errors=TRUE;  
        break } 
      
      if( control$safe ){
        new.kf=rerun.kf("V0", MLEobj.iter, iter)
        if(!new.kf$ok){
          stopped.with.errors=TRUE
          msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
          break
        }else{
          kf=new.kf$kf
          MLEobj.iter$kf=kf
          MLEobj.iter$logLik=kf$logLik
          Ey = MARSShatyt( MLEobj.iter )
          MLEobj.iter$Ey=Ey
          msg.kem=c(msg.kem, new.kf$msg.kem)
        }
      } 
    } #if not fixed V0
    
    
    ################
    # Get new A subject to its constraints (update of R will use this)
    ##############################################################
    if( !fixed[["A"]] ){ #if there is anything to update
      #note if Z and f.a are constant then we can write this as numer = (Ey$ytT - Z %*% kf$xtT)%*%matrix(1,dim(kf$xtT)[2],1)-TT*f$A
      numer=denom=0
      Z=par1$Z #reset
      starR=star[["R"]][,,1]
      for(i in 1:TT){
        if(time.varying[["Z"]] & i>1) Z = parmat(MLEobj.iter,"Z",t=i)$Z
        if(time.varying[["A"]] | i==1){
          dA=sub3D(d[["A"]],t=i)
          fA=sub3D(f$A,t=i)
        }
        if(time.varying[["R"]]) starR=star[["R"]][,,i]
        numer = numer + crossprod(dA, starR%*%(Ey[["ytT"]][,i,drop=FALSE]-Z%*%kf[["xtT"]][,i,drop=FALSE]-fA))
        denom = denom + crossprod(dA, starR%*%dA)
      }
      if(length(denom)==1){ denom = try(1/denom) }else{ denom=try(chol2inv(chol( denom ))) }
      if(inherits(denom, "try-error")){
        stop.msg = paste("Stopped at iter=",iter," in MARSSkem at A update. denom is not invertible. \n If some of your R diagonals equal 0, then A elements corresponding to R==0 cannot be estimated.\n The problem may be with your D matrix (if you have one) also. Type MARSSinfo('denominv') for more info. \n", sep="")
        stopped.with.errors=TRUE;  break }
      MLEobj.iter$par$A = denom%*%numer
      
      #not sure why denom%*%numer would ever not be a matrix
      if(!is.matrix(MLEobj.iter$par$A) ) MLEobj.iter$par$A=matrix(MLEobj.iter$par$A,dim(d$A)[2],1)
      par1$A = parmat(MLEobj.iter,"A",t=1)$A
      
      #~~~~~~~~Error checking  
      if( control$safe ){
        new.kf=rerun.kf("A", MLEobj.iter, iter)
        if(!new.kf$ok){
          stopped.with.errors=TRUE
          msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
          break
        }else{
          kf=new.kf$kf
          MLEobj.iter$kf=kf
          MLEobj.iter$logLik=kf$logLik
          Ey = MARSShatyt( MLEobj.iter )
          MLEobj.iter$Ey=Ey
          msg.kem=c(msg.kem, new.kf$msg.kem)
        }
      }
    } #A not fixed
    
    ################
    # Get new U subject to its constraints (update of Q and B will use this)
    ########################################################################
    if( !fixed[["U"]] ){ #if there is anything to update
      numer = matrix(0,m,1); denom = matrix(0,m,m) #this is the start if kf.x0="x10"
      fU=sub3D(f$U,t=1); dU=sub3D(d$U,t=1)
      B=par1$B; Z=par1$Z; A=par1$A   #reset
      Qinv = star$Q[,,1]; Rinv = star$R[,,1]
      #U.degen.update = FALSE #CUT?
      diag.Q=1-takediag(IIz$Q[,,1]); diag.V0=1-takediag(IIzV0)
      nQ0=sum(diag.Q==0)
      #if(any(diag.Q==0)) U.degen.update=!all(d$U[diag.Q==0,,]==0)    #CUT?   
      numer=denom=0
      if(kf.x0=="x00") hatxt0=kf$x0T else hatxt0=kf$xtT[,1,drop=FALSE]
      E.x0 = IImIIzV0%*%hatxt0+IIzV0%*%par1$x0
      AdjM=B; AdjM[AdjM!=0]=1
      if(kf.x0=="x00"){
        #set up values for t=1
        Bstar=B; Bstar.tm=IIm
        fstar=fU; fstar.tm=0*fU
        Dstar=dU; Dstar.tm=0*dU
        Mt=AdjM
        IId.tm=IIm #t=0; no w's
        IId=makediag(1-diag.Q) #only can be w free if Q==0
        if(any(diag.Q==0)&any(diag.Q!=0)) #which did not pick up a zero from B
          IId[diag.Q==0,diag.Q==0][1 + 0:(nQ0 - 1)*(nQ0 + 1)]=apply(Mt[diag.Q==0,diag.Q!=0,drop=FALSE]==0,1,all)
        #x_1-B(I-Id)xtm-B Id (B* E.x0 + f*)-fU; f*=0, B*=I; xtm=E.x0 so reduces to the following
        Delta3 = kf$xtT[,1,drop=FALSE]-B%*%E.x0-fU
        Delta4 = dU #since IId.tm and Bstar.tm are IIm
        numer = numer + crossprod(Delta4, Qinv%*%Delta3) 
        denom = denom + crossprod(Delta4, Qinv%*%Delta4) 
      }else{  #x10
        Bstar=IIm; fstar=0*fU; Dstar=0*dU
        IId.tm=0*IIm; IId=IIm
        Mt=IIm
      }
      if(any(IId==1)){ #again t=1
        Delta1=Ey$ytT[,1,drop=FALSE]-Z%*%((IIm-IId)%*%kf$xtT[,1,drop=FALSE])-Z%*%(IId%*%(Bstar%*%E.x0+fstar))-A
        Delta2=Z%*%IId%*%Dstar
        numer = numer + crossprod(Delta2, Rinv%*%Delta1) 
        denom = denom + crossprod(Delta2, Rinv%*%Delta2) 
      }
      for(t in 2:TT){ 
        if(time.varying$U){ fU=sub3D(f$U,t=t); dU=sub3D(d$U,t=t) }
        if(time.varying$B) B = parmat(MLEobj.iter,c("B"),t=t)$B
        if(time.varying$A) A = parmat(MLEobj.iter,c("A"),t=t)$A
        if(time.varying$Z) Z = parmat(MLEobj.iter,c("Z"),t=t)$Z
        if(time.varying$R) Rinv = star$R[,,t]
        if(time.varying$Q) Qinv = star$Q[,,t]
        fstar.tm=fstar
        fstar=B%*%fstar + fU
        Dstar.tm=Dstar
        Dstar=B%*%Dstar + dU
        Bstar.tm=Bstar
        Bstar=B%*%Bstar
        if(t<=(m+1)){
          IId.tm=IId
          IId=makediag(1-diag.Q) #only can be w free if Q==0
          if(any(diag.Q==0)&any(diag.Q!=0)){
            Mt=AdjM%*%Mt 
            IId[diag.Q==0,diag.Q==0][1 + 0:(nQ0 - 1)*(nQ0 + 1)]=apply(Mt[diag.Q==0,diag.Q!=0,drop=FALSE]==0,1,all)
          }
        }
        if(any(IId==1)){
          Delta1=Ey$ytT[,t,drop=FALSE]-Z%*%((IIm-IId)%*%kf$xtT[,t,drop=FALSE])-Z%*%(IId%*%(Bstar%*%E.x0+fstar))-A
          Delta2=Z%*%IId%*%Dstar
          numer = numer + crossprod(Delta2, Rinv%*%Delta1) 
          denom = denom + crossprod(Delta2, Rinv%*%Delta2) 
        }
        Delta3 = kf$xtT[,t,drop=FALSE]-B%*%((IIm-IId.tm)%*%kf$xtT[,t-1,drop=FALSE])-B%*%(IId.tm%*%(Bstar.tm%*%E.x0+fstar.tm))-fU
        Delta4 = dU + B%*%IId.tm%*%Dstar.tm
        numer = numer + crossprod(Delta4, Qinv%*%Delta3) 
        denom = denom + crossprod(Delta4, Qinv%*%Delta4) 
      } #for i
      if(length(denom)==1){ denom = 1/denom }else{ denom=try(chol2inv(chol( denom ) ),silent=TRUE) }
      if(inherits(denom, "try-error") | (length(denom)==1 && denom==0)){
        stop.msg = paste("Stopped at iter=",iter," in MARSSkem at U update. denom is not invertible.\n This means some of the U (+ C) terms cannot be estimated.\n Type MARSSinfo('denominv') for more info. \n", sep="")
        stopped.with.errors=TRUE
        break
      }
      MLEobj.iter$par$U = denom%*%numer
      if( !is.matrix(MLEobj.iter$par$U) ) MLEobj.iter$par$U=matrix(MLEobj.iter$par$U,dim(d$U)[2],1)
      par1$U = parmat(MLEobj.iter,"U",t=1)$U  
      
      #~~~~~~~~Error checking  
      if( control$safe &&  !fixed[["U"]] ){
        new.kf=rerun.kf("U", MLEobj.iter, iter)
        if(!new.kf$ok){
          stopped.with.errors=TRUE
          msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
          break
        }else{
          kf=new.kf$kf
          MLEobj.iter$kf=kf
          MLEobj.iter$logLik=kf$logLik
          Ey = MARSShatyt( MLEobj.iter )
          MLEobj.iter$Ey=Ey
          msg.kem=c(msg.kem, new.kf$msg.kem)
        }
      }  
    } #any U not fixed
    
    ################
    # Get new B subject to its constraints
    ################################################################
    if( !fixed[["B"]] ) {
      #t.kf.xtT = t(kf$xtT) #move t() out of for loop
      #need these for t=1 whether kf.x0=x00 or not
      U=par1$U    #reset
      starQ=sub3D(star$Q,t=1)
      dB=sub3D(d$B,t=1); fB=sub3D(f$B,t=1)
      if(kf.x0=="x00"){  #prior is defined as being E[x(t=0)|y(t=0)]; xtt[0]=x0; Vtt[0]=V0
        hatxtm =  IImIIzV0%*%kf$x0T + IIzV0%*%par1$x0
        hatVtm = IImIIzV0%*%kf$V0T%*%IImIIzV0 + IIzV0%*%par1$V0%*%IIzV0
        hatxt = kf$xtT[,1,drop=FALSE]
        Ptm = hatVtm + tcrossprod(hatxtm,kf$x0T) #note def of t.hatxtm = kf$x0T not t(hatxtm)
        Pttm = kf$Vtt1T[,,1] + tcrossprod(hatxt,kf$x0T) #note def of t.hatstm for t=0
        kronPtmQ = kronecker(Ptm,starQ)
        denom = crossprod(dB, kronPtmQ%*%dB)
        numer = crossprod(dB, (vec(starQ%*%(Pttm - tcrossprod(U,kf$x0T))) - kronPtmQ%*%fB))
      }else{ #prior is defined as being E[x(t=1)|y(t=0)]; xtt1[1]=x0; Vtt1[1]=V0
        denom = numer = 0 #see Ghahramani and Hinton treatment.  Summation starts at t=2
      }        
      
      for (i in 2:TT) { 
        hatxtm = kf$xtT[,i-1,drop=FALSE]
        hatVtm = kf$VtT[,,i-1] 
        hatxt = kf$xtT[,i,drop=FALSE]
        Ptm = hatVtm + tcrossprod(hatxtm) 
        Pttm = kf$Vtt1T[,,i] + tcrossprod(hatxt,hatxtm) 
        if(time.varying$Q) starQ=sub3D(star$Q,t=i)
        kronPtmQ = kronecker(Ptm,starQ)
        if(time.varying$B){
          dB=sub3D(d$B,t=i)
          fB=sub3D(f$B,t=i)
        }
        if(time.varying$U) U = parmat(MLEobj.iter,"U",t=i)$U
        denom = denom + crossprod(dB, kronPtmQ%*%dB)
        numer = numer + crossprod(dB, vec(starQ%*%(Pttm - tcrossprod(U, hatxtm) )) - kronPtmQ%*%fB )
      } #for i
      if(length(denom)==1){ denom = try(1/denom) }else{ denom=try(chol2inv(chol( denom ) )) }
      if(inherits(denom, "try-error")){
        stop.msg = paste("Stopped at iter=",iter," in MARSSkem at B update. denom is not invertible.\n Type MARSSinfo('denominv') for more info. \n", sep="")
        stopped.with.errors=TRUE;  break }
      MLEobj.iter$par$B = denom%*%numer
      if( !is.matrix(MLEobj.iter$par$B) ) MLEobj.iter$par$B=matrix(MLEobj.iter$par$B,dim(d$B)[2],1)
      par1$B = parmat(MLEobj.iter,"B",t=1)$B
      
      #~~~~~~~~Error checking      
      if( control$safe ){
        new.kf=rerun.kf("B", MLEobj.iter, iter)
        if(!new.kf$ok){
          stopped.with.errors=TRUE
          msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
          break
        }else{
          kf=new.kf$kf
          MLEobj.iter$kf=kf
          MLEobj.iter$logLik=kf$logLik
          Ey = MARSShatyt( MLEobj.iter )
          MLEobj.iter$Ey=Ey
          msg.kem=c(msg.kem, new.kf$msg.kem)
        }
      } 
      if( control$trace>0 ) {
        Ck = kappa(denom)
        if(Ck>condition.limit) msg.kem=c(msg.kem,paste("iter=",iter," Unstable B estimate because P_{t-1,t-1} is ill-conditioned. C =",round(Ck), "\n", sep=""))
        for(i in 1:max(dim(f$B)[3],dim(d$B)[3])){
          parB=parmat(MLEobj.iter,"B",t=i)$B
          if(any(abs(Re(eigen(parB,only.values=TRUE)$values))>1)) msg.kem=c(msg.kem,paste("iter=",iter,",t=",i," B update is outside the unit circle.", "\n", sep=""))
        }
      }
    } #if !is.fixed B
    
    
    ################
    # Get new Z subject to its constraints
    ################################################################
    if( !fixed[["Z"]] ){
      numer = denom = 0
      A = par1$A #reset
      starR=star$R[,,1]
      for (i in 1:TT) {
        hatxt = kf$xtT[,i,drop=FALSE]
        Pt = kf$VtT[,,i] + tcrossprod(hatxt) 
        hatyxt = Ey$yxtT[,,i]
        if(time.varying$A & i>1) A=parmat(MLEobj.iter,"A",t=i)$A
        if(time.varying$R) starR=star$R[,,i]
        if(time.varying$Z | i==1){
          fZ=sub3D(f$Z,t=i)
          dZ=sub3D(d$Z,t=i)
        }
        PkronR = kronecker( Pt, starR )
        numer = numer + crossprod(dZ, vec(starR%*%(hatyxt-tcrossprod(A, hatxt))) - PkronR%*%fZ)
        denom = denom + crossprod(dZ, PkronR%*%dZ)
      } #for i
      if(length(denom)==1){ denom = try(1/denom) }else{ denom=try(chol2inv(chol( denom ) )) }
      if(inherits(denom, "try-error")){ 
        stop.msg = paste("Stopped at iter=",iter," in MARSSkem in Z update.  denom is not invertible.\n Type MARSSinfo('denominv') for more info. \n", sep="")
        stopped.with.errors=TRUE;  break }
      MLEobj.iter$par$Z = denom%*%numer
      
      #not sure this is needed; is guarding against R returning a vector
      if( !is.matrix(MLEobj.iter$par$Z) ) MLEobj.iter$par$Z=matrix(MLEobj.iter$par$Z,dim(d$Z)[2],1)
      par1$Z = parmat(MLEobj.iter,"Z",t=1)$Z            
      
      #Start~~~~~~~~~~~~Error checking
      if( control$safe ){
        new.kf=rerun.kf("Z", MLEobj.iter, iter)
        if(!new.kf$ok){
          stopped.with.errors=TRUE
          msg.kf=c(msg.kf, new.kf$msg.kf); stop.msg=new.kf$stop.msg
          break
        }else{
          kf=new.kf$kf
          MLEobj.iter$kf=kf
          MLEobj.iter$logLik=kf$logLik
          Ey = MARSShatyt( MLEobj.iter )
          MLEobj.iter$Ey=Ey
          msg.kem=c(msg.kem, new.kf$msg.kem)
        }
      }  
      if( control$trace>0 ) {
        Ck = kappa(denom)
        if(Ck>condition.limit) msg.kem=c(msg.kem,paste("iter=",iter," Unstable Z estimate because P_{t,t} is ill-conditioned. C =",round(Ck), sep=""))
      }
    }#if !is.fixed Z 
  }  # end inner iter loop
  if(control$silent==2) cat("\n")
  
  #prepare the MLEobj to return which has the elements set here
  MLEobj.return = MLEobj
  MLEobj.return$iter.record = iter.record
  MLEobj.return$numIter = iter
  
  if(stopped.with.errors){
    if( control$silent==2 ) cat("Stopped due to numerical instability or errors. Print $errors from output for info or set silent=FALSE.\n")      
    #print brief msg.  Full msg printed if silent=F
    msg=c(stop.msg,"par, kf, states, iter, loglike are the last values before the error.\n")
    if(!control$safe) {
      msg=c(msg,"Try control$safe=TRUE which uses a slower but slightly more robust algorithm.\n")
    }
    if(!control$trace>0) {
      msg=c(msg,"Use control$trace=1 to generate a more detailed error report. See user guide for insight.\n")
    }
    ## Attach any algorithm errors to the MLEobj
    if(control$trace>0 && !is.null(msg.kem)) msg=c(msg,"\nMARSSkem errors. Type MARSSinfo() for help.\n",msg.kem)
    if(control$trace>0 && !is.null(msg.kf)) msg=c(msg,"\nMARSSkf errors. Type MARSSinfo() for help.\n",msg.kf,"\n")    
    MLEobj.return$errors=msg
    
    MLEobj.return$par=MLEobj.iter$par
    MLEobj.return$kf = kf.last
    MLEobj.return$states = kf.last$xtT
    MLEobj.return$convergence = 52
    MLEobj.return$logLik = MLEobj.iter$logLik
    return(MLEobj.return)
  } #stopped with errors
  
  ########### Did not stop with errors 
  ## Set the convergence information
  ## Output depends on how it converged and how iterations were determined.  Possible conv.test$convergence values are
  #  MLEobj.iter$conv.test$convergence==72  means no convergense testing done at all because exited with errors BEFORE minit
  #  MLEobj.iter$conv.test$convergence==4 means stopped by hitting maxit; abstol not reached so no log-log testing done
  #  MLEobj.iter$conv.test$convergence==3 means abstol reached but no log-log test since control$min.iter.conv.test not reached
  #  MLEobj.iter$conv.test$convergence==0  means abstol reached and log-log test passed (=CONVERGED)
  #  MLEobj.iter$conv.test$convergence==1  means stopped by hitting maxit; abstol reached but NOT log-log
  #  MLEobj.iter$conv.test$convergence==-1 or -2  means log-log test returned errors
  catinfo=!control$silent || control$silent==2
  MLEobj.return$convergence = 72  #debugging should be changed below
  if(MLEobj.iter$conv.test$convergence==72){
    MLEobj.return$convergence = 52
    msg.conv=MLEobj.iter$conv.test$messages
    if( catinfo ) cat(paste("Error! EM algorithm exited at iter=",iter," before minit reached.\n Minit was ",control$minit,".\n",sep=""))
  }
  if(MLEobj.iter$conv.test$convergence<0){
    MLEobj.return$convergence = 62
    msg.conv=MLEobj.iter$conv.test$messages
    if( catinfo ) cat(paste("Error! EM algorithm exited due to errors reported by log-log test function.\n" ,sep=""))
  }
  #if it returns 4, then abstol never reached before maxit hit.  Abstol must be hit before log-log
  #test run.  Thus we need to run the log-log test here IF we have enough iterations
  if(MLEobj.iter$conv.test$convergence==4){
    tmp.msg=paste("Warning! Reached maxit before parameters converged. Maxit was ",control$maxit,".\n",sep="")
    if( iter>=control$min.iter.conv.test ){  
      #loglog test has not been run because abstol never reached
      loglog.test = loglog.conv.test(iter.record, iter, deltaT=control$conv.test.deltaT, tol=control$conv.test.slope.tol)
      MLEobj.iter$conv.test$messages=loglog.test$messages
      if(loglog.test$convergence==0){
        MLEobj.return$convergence=11 #log-log passed but not abstol
        if( catinfo ) cat(paste("Warning! log-log convergence only. Maxit (=",control$maxit,") reached before abstol convergence.\n",sep="")) 
      }    
      if(loglog.test$convergence<0){
        MLEobj.return$convergence=63
        if( catinfo ) cat(paste(tmp.msg, " abstol not reached and log-log convergence returned errors.\n",sep=""))
      }      
      if(loglog.test$convergence==1){
        MLEobj.return$convergence=1
        if( catinfo ) cat(paste(tmp.msg, " neither abstol nor log-log convergence tests were passed.\n",sep=""))
      }      
      if(loglog.test$convergence>1){ #loglog test should never return > 1
        MLEobj.return$convergence=72
        if( catinfo ) cat(paste(tmp.msg, " abstol not reached and log-log convergence returned errors.\n",sep=""))
      }     
    }else{ #can'd do log-log test
      MLEobj.return$convergence=12  #can't do the log-log test
      if( catinfo ) cat(paste(tmp.msg, " abstol not reached and no log-log test info since maxit less than min.iter.conv.test.\n",sep=""))
    }
  }
  #At this point $conv.test$convergence==4 still and $convergence is 52,62,72 errors;
  #1 neither abstol nor loglog; 11 loglog only; 12 no abstol and no info on loglog;
  
  #if conv.test$convergence == 3, abstol reached but no info on loglog since maxit < min.iter.conv.test
  if(MLEobj.iter$conv.test$convergence==3){
    tmp.msg=paste("Warning! Abstol convergence only. no info on log-log convergence.\n",sep="")
    MLEobj.return$convergence = 3
    if( catinfo ) cat(paste(tmp.msg, " Maxit (=",control$maxit,") < min.iter.conv.test (=",control$min.iter.conv.test,") so not log-log test.\n", sep=""))
  }
  #if conv.tst$convergence==1, then abstol reached but no loglog convergence
  if(MLEobj.iter$conv.test$convergence==1){
    MLEobj.return$convergence = 10
    msg.conv=MLEobj.iter$conv.test$messages
    if( catinfo ) cat(paste("Warning! Abstol convergence only. Maxit (=",control$maxit,") reached before log-log convergence.\n",sep=""))
  }
  if(MLEobj.iter$conv.test$convergence==0){
    MLEobj.return$convergence = 0
    if( catinfo ) {
      if(iter==control$minit){ 
        cat(paste("Success! algorithm run for ",iter," iterations. abstol and log-log tests passed.\n",sep=""))
      }else{ cat(paste("Success! abstol and log-log tests passed at ",iter," iterations.\n",sep="")) }
      if(control$conv.test.slope.tol>0.1) cat(paste("Alert: conv.test.slope.tol is ",control$conv.test.slope.tol,".\nTest with smaller values (<0.1) to ensure convergence.\n",sep="")) 
    } #!silent
  }  
  msg.conv=MLEobj.iter$conv.test$messages
  if(!is.null(msg.conv)) msg=c(msg, "\nConvergence warnings\n", msg.conv)
  ##############################################################
  
  ## Other misc output
  MLEobj.return$par=MLEobj.iter$par
  MLEobj.return$states = MLEobj.iter$kf$xtT
  MLEobj.return$logLik = MLEobj.iter$logLik
  
  if(!is.null(msg.kem)){ msg.kem=c("\nMARSSkem warnings. Type MARSSinfo() for help.\n", msg.kem); msg=c(msg, msg.kem) }
  if(!is.null(msg.kf)) { msg.kf=c("\nMARSSkf warnings. Type MARSSinfo() for help.\n", msg.kf); msg=c(msg, msg.kf) }
  if((!is.null(msg.kem) || !is.null(msg.kf)) && control$trace<1){  msg = c(msg,  "\nUse control$trace=1 to generate a more detailed error report.\n") }
  if((!is.null(msg.kem) || !is.null(msg.kf)) && (!control$silent || control$silent==2) ){
    cat("Alert: Numerical warnings were generated. Print the $errors element of output to see the warnings.\n")
  }
  
  ## Attach any algorithm errors to the MLEobj
  MLEobj.return$errors=msg
  
  return(MLEobj.return)
}

## Run log-log convergence diagnostics
# 0 converged; 1 not converged; negative problem
loglog.conv.test = function(iter.record, iter, params.to.test=c("Z","U","x0","R","Q","A","logLik"), deltaT=9, tol=0.5){
  if( !is.list(iter.record) || !all(c("par","logLik") %in% names(iter.record)) || 
        !any(params.to.test %in% c(names(iter.record$par),names(iter.record))) ||
        length(dim(iter.record$par))!=2 || dim(iter.record$par)[1]<=1 || is.null(colnames(iter.record$par)) ){ 
    msg="par list not a proper list (with par and logLik) or too short for conv test or has no column names.\n"
    return( list(convergence=-1, messages=msg) )
  }else {
    if("logLik" %in% params.to.test){
      #exp because we don't want the log of the log and subtract mean so exp(LL) doesn't = Inf
      #we are looking at slope so doesn't matter if we sub off the mean
      iter.record.par = cbind(iter.record$par,logLik=exp(iter.record$logLik-mean(iter.record$logLik))) 
      #iter.record.par = cbind(iter.record$par,logLik=iter.record$logLik) 
    }else iter.record.par=iter.record$par 
    names.iter=colnames(iter.record.par)
    names.sub=strsplit(names.iter,"\\.")
    num.names = length(names.sub)
    p.elems=NULL
    for(j in 1:num.names)p.elems=c(p.elems,names.sub[[j]][1])
    num.varcov = sum( p.elems %in% params.to.test )
    test.conv=rep(0,num.names)
    for( j in 1:num.names ){
      if( p.elems[j] %in% params.to.test ) {
        test.len2=dim(iter.record.par)[1]
        test.len1=max(1,test.len2-deltaT)
        test.len=(iter-min(test.len2-1, deltaT)):iter 
        test.par = abs(iter.record.par[test.len1:test.len2,j])
        if(any(test.par==0)) test.par = test.par+1   
        #test.loglog=lm(log(test.par)~log(test.len))
        test.loglog=(log(test.par[length(test.par)])-log(test.par[1]))/(log(test.len[length(test.len)])-log(test.len[1]))
        #test.conv[j]=test.loglog$coef[2]
        test.conv[j]=test.loglog
      }
    }   
  }
  if(any(is.na(test.conv))) {
    msg="The log-log degeneracy test produced NAs.\n"
    return( list(convergence=-2, messages=msg) )
  } 
  if(!is.null(test.conv) && !any(is.na(test.conv)) && any(abs(test.conv)>tol)){
    msg=paste("Warning: the ",names.iter[abs(test.conv)>tol]," parameter value has not converged.\n")
    msg=c(msg,"Type MARSSinfo(\"convergence\") for more info on this warning.\n")
    return( list(convergence=1, messages=msg, not.converged.params=names.iter[abs(test.conv)>tol], converged.params=names.iter[abs(test.conv)<=tol]) ) 
  }else { return( list(convergence=0, messages=NULL, not.converged.params=names.iter[abs(test.conv)>tol], converged.params=names.iter[abs(test.conv)<=tol] ) ) }  #0 means converged successfully
}

rerun.kf = function(elem, MLEobj, iter){    #Start~~~~~~~~Error checking
  if(iter==1) cvg2 = 1 + MLEobj$control$abstol
  msg.kem=NULL; msg.kf=NULL
  loglike.old = MLEobj$logLik 
  kf = MARSSkf( MLEobj )
  if(MLEobj$control$demean.states) {
    xbar = apply(cbind(kf$x0T,kf$xtT),1,mean)
    kf$xtT = kf$xtT-xbar
    kf$x0T = kf$x0T-xbar
  }
  if(!kf$ok){ 
    msg.kf=paste("iter=",iter," ", elem," update ",kf$errors,sep=""); 
    stop.msg = paste("Stopped at iter=",iter," in MARSSkem after ", elem," update: numerical errors in ",MLEobj$fun.kf,".\n",sep="")
    return(list(ok=FALSE, msg.kf=msg.kf, stop.msg=stop.msg)) }
  loglike.new = kf$logLik
  if(iter>1 && is.finite(loglike.old) == TRUE && is.finite(loglike.new) == TRUE ) cvg2 = loglike.new - loglike.old  
  if(iter > 2 & cvg2 < -sqrt(.Machine$double.eps)) {
    if(MLEobj$control$trace>0){
      msg.kem=paste("iter=",iter," LogLike DROPPED in ",elem," update. logLik old=", loglike.old, " new=", loglike.new,"\n", sep="")
    }else msg.kem = paste("MARSSkem: The soln became unstable and logLik DROPPED in the",elem, "updates.\n")
  }
  return(list(kf=kf, msg.kem=msg.kem, msg.kf=msg.kf, ok=TRUE))
}

degen.test = function(elem, MLEobj, iter){
  if( is.fixed(MLEobj$marss$free[[elem]]) ) return(list(MLEobj=MLEobj, msg=NULL, set.degen=FALSE))
  if( MLEobj$constr.type[[elem]]=="time-varying" ) return(list(MLEobj=MLEobj, msg=NULL, set.degen=FALSE))
  if( !MLEobj$control$allow.degen ) return(list(MLEobj=MLEobj, msg=NULL, set.degen=FALSE))
  if( iter<=MLEobj$control$min.degen.iter ) return(list(MLEobj=MLEobj, msg=NULL, set.degen=FALSE))
  if( !is.design(MLEobj$marss$free[[elem]], zero.rows.ok=TRUE) )
    return(list(MLEobj=MLEobj, msg=NULL, set.degen=FALSE))  #strict, i.e. only 0 and 1
  isDiag = function(x){ all(x[!diag(nrow(x))] == 0) }
  if( !isDiag(coef(MLEobj, type="matrix", form="marss")[[elem]] ) )
    return(list(MLEobj=MLEobj, msg=NULL, set.degen=FALSE))  #only allow setting of 0s if diagonal
  
  #diagonal, not fixed, not time-varying, allow.degen set, iter>min iter and free is a design matrix
  #So can proceed
  #if here then not time-varying
  degen.par= abs(MLEobj$par[[elem]]) < MLEobj$control$degen.lim
  msg.degen=NULL
  set.degen=FALSE
  dim.elem = attr(MLEobj$marss, "model.dims")[[elem]][1]
  if(any(degen.par)){
    for(i in which(degen.par)){
      MLEobj.tmp=MLEobj
      #set corresponding par to 0 and corresponding col of free to 0
      MLEobj.tmp$par[[elem]][i,1]=0
      MLEobj.tmp$marss$free[[elem]][,i,1]=0 #req not time-varying so set t=1
      #need to check that setting a R or Q diag to 0 doesn't lead to a improper model
      kemcheck=MARSSkemcheck(MLEobj.tmp)
      if(kemcheck$ok){
        new.kf = MARSSkf( MLEobj.tmp )
        loglike.old=MLEobj$logLik
        if(!new.kf$ok) msg.degen=c(msg.degen,paste("iter=",iter," MARSSkf returned error in attempt to set 0 diagonals for ", elem,"\n  ", new.kf$errors,"Perhaps Q and R are both going to 0?\n", sep="") ); 
        if(new.kf$ok && is.finite(loglike.old) && is.finite(new.kf$logLik) ) tmp.cvg2 = new.kf$logLik - loglike.old  else tmp.cvg2=Inf
        if(new.kf$ok && tmp.cvg2 < -sqrt(.Machine$double.eps)) {
          msg.degen=c(msg.degen,paste("iter=",iter," Setting diagonal to 0 blocked. logLik was lower in attempt to set 0 diagonals on ",elem," logLik old=", loglike.old, " new=", new.kf$logLik,  ", See MARSSinfo(\"",elem,"0blocked\").\n", sep=""))
        }
        if(new.kf$ok && tmp.cvg2 > -sqrt(.Machine$double.eps)) { #this means degenerate elem has lower LL, so accept it
          MLEobj=MLEobj.tmp
          MLEobj$kf=new.kf
          MLEobj$Ey=MARSShatyt( MLEobj ) #needs the updated kf
          if(MLEobj$control$demean.states) {
            xbar = apply(cbind(new.kf$x0T,new.kf$xtT),1,mean)
            MLEobj$kf$xtT = new.kf$xtT-xbar
            MLEobj$kf$x0T = new.kf$x0T-xbar
          }
          MLEobj$logLik=new.kf$logLik
          set.degen=TRUE
        } 
      }else{ msg.degen=c( msg.degen, paste("iter=",iter," Setting element of ",elem," to 0, blocked.  See MARSSinfo(\"",elem,"0blocked\"). The error is due to the following MARSSkemcheck errors.\n  MARSSkemcheck error: ", kemcheck$msg,sep="")) }
    } #for degen.par; do one by one
  } #update MLEobj
  #set.degen is a flag to say if any 0 elements were set
  return(list(MLEobj=MLEobj, msg=msg.degen, set.degen=set.degen))
}

