#######################################################################################################
#   MARSSkfss function
#   Kalman filter and smoother as presented in Shumway & Stoffer (2006)
#   ** All eqn refs are to 2nd ed of Shumway & Stoffer (2006): Time Series Analysis and Its Applications
# 5-17-12  I removed the tests re inversion of Z when R=0 and put that in MARSSkem
#######################################################################################################
MARSSkfss = function( MLEobj ) {
    condition.limit=1E10
    condition.limit.Ft=1E5 #because the Ft is used to compute LL and LL drop limit is about 2E-8

    MODELobj = MLEobj$marss
    debugkf = MLEobj$control$trace
    n=dim(MODELobj$data)[1]; TT=dim(MODELobj$data)[2]; m=dim(MODELobj$fixed$x0)[1]

    #create the YM matrix
    YM=matrix(as.numeric(!is.na(MODELobj$data)),n,TT)
    #Make sure the missing vals in y are zeroed out if there are any
    y=MODELobj$data
    y[YM==0]=0
    
    if(MODELobj$tinitx==1){ init.state="x10" }else{ init.state="x00" }
	  msg=NULL
    #Construct needed identity matrices
    I.m = diag(1,m); I.n = diag(1,n)
    
    #initialize matrices
    # for notation purposes, 't' represents current point in time, 'T' represents the length of the series
    Vtt = Vtt1 = VtT = Vtt1T =J = array(0,dim=c(m,m,TT))     
    # Vtt is analogous to S&S Ptt, var[xt,xt|y(1:t)];  Vtt1 is analogous to S&S Ptt1, cov[xt,xt1|y(1:t)]
    # VtT and Vtt1T are for the backwards smoother; VtT is var[xt,xt|y(1:T)]  and Vtt1T is cov[xt,xt1|y(1:T)]
    # J see eqn 6.49
    xtt = xtt1 = xtT = matrix(0,m,TT)       
    # xtt is E[x(t) | Y(t)]; xtt1 is E[x(t) | Y(t-1)]; xtT is E[x | y(1:T)]
    vt = matrix(0,n,TT)     # these are innovations, parentheses in 6.21; vt equivalent to epsilon, eqn 6.62 
    Ft = array(0,dim=c(n,n,TT))      # used for likelihood, Ft equivalent to sigma matrix eqn 6.62
    Kt = array(0, dim=c(m,n,TT))     # 3D matrix of Kalman gain, EW added 11/14/08

    #Note diff in param names from S&S;B=Phi, Z=A, A not in S&S
    model.elem = names(MODELobj$fixed)
    time.varying = c()
    for(elem in model.elem){
      if( (dim(MODELobj$free[[elem]])[3] != 1) | (dim(MODELobj$fixed[[elem]])[3] != 1))  #not time-varying
           time.varying = c(time.varying, elem)
    }
    pari=parmat(MLEobj,t=1)
    Z=pari$Z; A=pari$A; B=pari$B; U=pari$U; x0=pari$x0; 
    R=tcrossprod(pari$H %*% pari$R, pari$H) 
    Q=tcrossprod(pari$G %*% pari$Q, pari$G)
    V0=tcrossprod(pari$L%*%pari$V0, pari$L)
    if(n==1){ diag.R=unname(R) }else{ diag.R = unname(R)[1 + 0:(n - 1)*(n + 1)] }
    #Check that if any R are 0 then model is solveable
    OmgRVtt = I.m; diag.OmgRVtt = rep(1,m)
    if( any(diag.R==0) ){
      Z.R0=Z[diag.R==0,,drop=FALSE]; Z.R0.n0=Z.R0[,colSums(Z.R0)!=0,drop=FALSE] # The Z cols where there is a val
      if(dim(Z.R0.n0)[1]==dim(Z.R0.n0)[2]){ #then it must be invertible
         Ck=try(solve(Z.R0.n0))
         if(class(Ck)=="try-error"){
           return(list(ok=FALSE,errors="Some R diagonal elements are 0, but Z is such that model is indeterminate in this case."))
         }else{ OmgRVtt = diag(ifelse(colSums(Z.R0)!=0,0,1),m) }#mxm; sets the p elem of Vtt to 0 because these are determined by data
      }else{ #then num rows must be greater than num cols
         if(dim(Z.R0.n0)[1]>dim(Z.R0.n0)[2]){
           return(list(ok=FALSE,errors="Some R diagonal elements are 0, and Z is such that model is over-determined."))
         }else{
           OmgRVtt = diag(ifelse(colSums(Z.R0)!=0,0,1),m)
         }
      }
      diag.OmgRVtt = takediag(OmgRVtt)
    }
    if(m==1){ diag.Q=unname(Q) }else{ diag.Q = unname(Q)[1 + 0:(m - 1)*(m + 1)] }
    t.B = matrix(B,m,m,byrow=TRUE)   #faster transpose
        
    ##############################################
    #FORWARD PASS (K filter) gets you E[x(t) given y(1:t)]
    ##############################################
    # In the following, vt and Ft are needed for the likelihood calculation
		# the missing values will contribute 0.0 for the LL calc
    # R_mod is needed for the corrected likelihood calculation when there are missing values
		# See section 12.3 in Time Series: Theory and Methods (1991) Peter J. Brockwell, Richard A. Davis
		# put 1's on the diagonal where there are missing values and zero out the rows and columns

    for (t in 1:TT) {
      if(length(time.varying)!=0){
        #update the time.varying ones
        pari[time.varying]=parmat(MLEobj,time.varying,t=t)
        Z=pari$Z; A=pari$A; B=pari$B; U=pari$U; #update
        if("R" %in% time.varying | "H" %in% time.varying){
          R=tcrossprod(pari$H %*% pari$R, pari$H) 
          if(n==1){ diag.R=unname(R) }else{ diag.R = unname(R)[1 + 0:(n - 1)*(n + 1)] }
          #Check that if any R are 0 then model is solveable
          OmgRVtt = I.m; diag.OmgRVtt = rep(1,m)
          if(any(diag.R==0)){
            Z.R0=Z[diag.R==0,,drop=FALSE]; Z.R0.n0=Z.R0[,colSums(Z.R0)!=0,drop=FALSE] # The Z cols where there is a val
            if(dim(Z.R0.n0)[1]==dim(Z.R0.n0)[2]){ #then it must be invertible
              Ck=try(solve(Z.R0.n0))
              if(class(Ck)=="try-error"){
               return(list(ok=FALSE,errors="Some R diagonal elements are 0, but Z is such that model is indeterminate in this case."))
              }else{ OmgRVtt = diag(ifelse(colSums(Z.R0)!=0,0,1),m) }#mxm; sets the p elem of Vtt to 0 because these are determined by data
            }else{ 
              if(dim(Z.R0.n0)[1]>dim(Z.R0.n0)[2]){
                return(list(ok=FALSE,errors="Some R diagonal elements are 0, and Z is such that model is over-determined."))
              } 
            }
          diag.OmgRVtt = takediag(OmgRVtt)
          }
        } #if R in time.varying
        if("Q" %in% time.varying | "G" %in% time.varying){
          Q=tcrossprod(pari$G %*% pari$Q, pari$G)
          if(m==1){ diag.Q=unname(Q) }else{ diag.Q = unname(Q)[1 + 0:(m - 1)*(m + 1)] }
        }
        if("B" %in% time.varying) t.B = matrix(B,m,m,byrow=TRUE)   #faster transpose
      }
    #missing value modifications per S&S2006 eq 6.78
    if(any(YM[,t]==0)){
      Mt = I.n; Mt[YM[,t]==0,]=0  #much faster than makediag(YM)
      I.2 = I.n-Mt  
      Zt = Mt%*%Z #If Y missing, that row is 0 in Zt
      At = Mt%*%A
      Omg1=I.n[YM[,t]==1,,drop=FALSE]
      t.Omg1 = I.n[,YM[,t]==1,drop=FALSE] 
      #per 6.78 in Shumway and Stoffer
      #Need to 0 out the covariance between R assoc with non-missing and R assoc with missing values
      Rt = Mt%*%R%*%Mt + I.2%*%R%*%I.2    
    }else { Zt=Z; Rt=R; At=A; Omg1=I.n; t.Omg1=I.n; Mt=I.n }
    if(m*n==1) t.Zt = Zt else t.Zt = matrix(Zt,m,n,byrow=TRUE) #faster transpose
    
    #t=1 treatment depends on how you define the initial condition.  Either as x at t=1 or x at t=0
    if(t==1) {
        if(init.state=="x00") {
          xtt1[,1] = B%*%x0 + U   #Shumway and Stoffer treatment of initial states # eqn 6.19   (pi is defined as t=0)
          Vtt1[,,1] = B%*%V0%*%t.B + Q          # eqn 6.20
        }
        if(init.state=="x10") {    #Ghahramani treatment of initial states uses x10 and has no x00 (pi is defined as t=1)
         xtt1[,1] = x0         
         Vtt1[,,1] = V0
        }
    }else {   #t!=1
       xtt1[,t] = B%*%xtt[,t-1,drop=FALSE] + U  #xtt1 denotes x_t^(t-1), eqn 6.19; B here is B[t]
       Vtt1[,,t] = B%*%Vtt[,,t-1]%*%t.B + Q     #eqn 6.20; B here is B[t]
    }
    if(m!=1) Vtt1[,,t] = symm(Vtt1[,,t])   #in general Vtt1 is not symmetric but here it is since Vtt and Q are

    #Set up the inverse needed in Kt (part corresponding to no missing values)
    #This is used in Kt, if Vtt1=0, then no info from y on those xt and corrs Kt rows =0 since Kt=Vtt1*t.Z*siginv
    siginv1 = Zt%*%Vtt1[,,t]%*%t.Zt + Rt
    # bracketed piece of eqn 6.23 modified per 6.78
    # Because R diag might be 0, the bracketed bit might have 0 diagonals.  Inv by pcholinv deals with this
    # by putting 0 row/cols where 0s appear on diagonal

    if(debugkf == -1){ siginv = pcholinv(siginv1) #skip error-checking
    }else{
        siginv=try(pcholinv(siginv1), silent = TRUE)      
        if(n==1) diag.siginv1 = unname(siginv1) else diag.siginv1=unname(siginv1)[1 + 0:(n - 1)*(n + 1)]   #much faster way to get the diagonal
        #Catch errors before entering chol2inv
        if(class(siginv)=="try-error") {
          if(any(diag.siginv1!=0)){
            Ck1 = try(kappa(siginv1[diag.siginv1!=0,diag.siginv1!=0,drop=FALSE]))
            Ck1 = ifelse(class(Ck1)=="try-error","Inf",round(Ck1))
            msg1=paste("Condition num. of siginv1[t=",t,"] = ",Ck1," ",sep="")
          }
          msg2=""
          if(any(diag.R!=0) & (t==1 || "R"%in%time.varying) ){
            Ck4 = try(kappa((diag(1,n)[diag.R!=0,])%*%R%*%t(diag(1,n)[diag.R!=0,]))) 
            Ck4 = ifelse(class(Ck4)=="try-error","Inf",round(Ck4))
            msg2=paste("Condition num. of R = ",Ck4," ",sep="")
          }
          return(list(ok=FALSE, errors=paste("Stopped in MARSSkfss: chol(Z%*%Vtt1[,,",t,"]%*%t(Z)+R) error. ",msg1," ", msg2,"\n",sep="") ) )       
        }
     }  #don't error check if debugkf=-1
    ####### End of Error-checking for this section
    
    if(n!=1) siginv =symm(siginv)
    Kt[,,t] =  Vtt1[,,t]%*%t.Zt%*%siginv
    Kt.tmp = sub3D(Kt,t=t) # stop R from changing matrix dim; drop=FALSE won't work here

    vt[,t] = y[,t,drop=FALSE] - (Zt%*%xtt1[,t,drop=FALSE]+At) #need to hold on to this for loglike calc; will be 0 when y is missing
    # eqn 6.21
    xtt[,t]=xtt1[,t,drop=FALSE] + Kt.tmp%*%vt[,t,drop=FALSE]     
    Vtt[,,t] = Vtt1[,,t]-Kt.tmp%*%Zt%*%Vtt1[,,t]  # eqn 6.22, detail after 6.28, modified Z per 6.78
    if(m!=1) Vtt[,,t] = symm(Vtt[,,t]) #to ensure its symetric
    #zero out rows cols as needed when R diag = 0
    OmgRVtt.t = OmgRVtt
    if(any(diag.OmgRVtt==0)) diag(OmgRVtt.t) = diag.OmgRVtt + t(!(Z==0))%*%(diag.R==0 & YM[,t]==0)
    Vtt[,,t] = OmgRVtt.t%*%Vtt[,,t]%*%OmgRVtt.t  
    
    # Variables needed for the likelihood calculation; see comments above
    R_mod = (I.n-Mt) + Mt%*%R%*%Mt #not in S&S; see MARSS documention per LL calc when missing values; R here is R[t]
    Ft[,,t] = Zt%*%Vtt1[,,t]%*%t.Zt+R_mod #need to hold on to this for loglike calc ; 1 on diagonal when y is missing
    if(n!=1) Ft[,,t] = symm(Ft[,,t]) #to ensure its symetric
       
    ####### Attach warnings to output if filter is becoming numerically unstable
    if(debugkf>0) {
          Ck1=Ck2=Ck3=Ck4=1
          if(!all(diag.siginv1==0)) Ck1 = kappa(siginv1[diag.siginv1!=0,diag.siginv1!=0,drop=FALSE])
          if(!all(diag.Q==0)) Ck2 = kappa(Vtt1[diag.Q!=0,diag.Q!=0,t])
          if(!all(diag.R==0) & t>1) Ck3 = kappa(Ft[diag.R!=0,diag.R!=0,t])
          if(Ck1>condition.limit && !all(Kt[,,t]==0) ) 
             msg=rbind(msg,paste("MARSSkfss: solution is becoming unstable.  Condition num. of siginv1[t=",t,"] = ",round(Ck1),"\n",sep=""))
          if(Ck2>condition.limit && t>1) 
             msg=rbind(msg,paste("MARSSkfss: solution is becoming unstable.  Condition num. of Vtt1[t=",t,"] = ",round(Ck2),"\n",sep=""))
          if( Ck3>condition.limit.Ft ){
             if(!all(diag.R==0)) Ck4 = kappa( R[diag.R!=0,diag.R!=0] )
             msg=rbind(msg,paste("MARSSkfss: logLik computation is becoming unstable.  Condition num. of Sigma[t=",t,"] = ",round(Ck3)," and of R = ",round(Ck4),".\n",sep=""))
             }
     }
     #Abandon if solution is so unstable that Vtt diagonal became negative
     if(m==1) diag.Vtt = unname(Vtt[,,t]) else diag.Vtt=unname(Vtt[,,t])[1 + 0:(m - 1)*(m + 1)]   #much faster way to get the diagonal
     if( any(diag.Vtt<0) ){
          return(list(ok=FALSE, 
          errors=paste("Stopped in MARSSkfss: soln became unstable and negative values appeared on the diagonal of Vtt at t=",t,".\n",sep="") ) )
     }
    ####### End Error-checking

    } #End of the Kalman filter recursion (for i to 1:TT)

    ######################################################
    #BACKWARD PASS (Kalman smoother) gets you E[x(t)|y(1:T)] from E[x(t)|y(1:t)]
    ######################################################
    xtT[,TT] = xtt[,TT,drop=FALSE]
    VtT[,,TT] = Vtt[,,TT]
    #indexing is 0 to T for the backwards smoother recursions
    s = seq(TT,2)
    for(i in 1:(TT-1)) {
      t=s[i]
      # Zt = Z; Zt[YM[,t]==0,]=0   #MUCH faster than defining Mt using diag(YM); commented out since doesn't seem to be used
      if("B" %in% time.varying){
        B = parmat(MLEobj, "B", t=t)$B  #t since in 6.49, B[t] appears
        if(m==1) t.B=B else t.B = matrix(B,m,m,byrow=TRUE) 
      }
      if("Q" %in% time.varying){
        Q = parmat(MLEobj, "Q", t=t)$Q
        if(m==1){ diag.Q=unname(Q) }else{ diag.Q = unname(Q)[1 + 0:(m - 1)*(m + 1)] } 
      }
      #deal with any 0s on diagonal of Vtt1; these can arise due to 0s in V0, B, + Q
      #0s on diag of Vtt1 will break the Kalman smoother if t>1
      if(m==1) diag.Vtt1 = unname(Vtt1[,,t]) else diag.Vtt1=unname(Vtt1[,,t])[1 + 0:(m - 1)*(m + 1)]   #much faster way to get the diagonal
      if( any(diag.Vtt1<0) ) #abandon if problems like this
          return(list(ok=FALSE, 
          errors=paste("Stopped in MARSSkfss: soln became unstable and negative values appeared on the diagonal of Vtt1.\n") ) )
      
      #Error-checking for 0s on diagonal of Vtt1 that they are allowed 
      if(debugkf!=-1 & any(diag.Vtt1==0) ) {  
        #deal with 0s that are ok if there are corresponding 0s on Q diagonal
        Q0s=all(which(diag.Vtt1==0)%in%which(diag.Q==0))
        #Q0s=identical(which(diag.Q==0),which(diag.Vtt1==0))
        if(!Q0s && (init.state=="x00" || (init.state=="x10" && t>1)) ){
          return(list(ok=FALSE, errors=paste("Stopped in MARSSkfss: soln became unstable when zeros appeared on the diagonal of Vtt1 at t=",t,".\n") ) )
        }
      }
      if(m==1){ Vinv=pcholinv(matrix(Vtt1[,,t],1,1)) 
      }else{
        Vinv=pcholinv(Vtt1[,,t])
        Vinv = symm(Vinv)  #to enforce symmetry after chol2inv call
      }
      J[,,t-1] = Vtt[,,t-1]%*%t.B%*%Vinv  # eqn 6.49 and 1s on diag when Q=0; Here it is t.B[t]

      xtT[,t-1] = xtt[,t-1,drop=FALSE] + J[,,t-1]%*%(xtT[,t,drop=FALSE]-xtt1[,t,drop=FALSE])     # eqn 6.47
      if(m==1) t.J = J[,,t-1] else t.J = matrix(J[,,t-1],m,m,byrow=TRUE) #faster transpose
      VtT[,,t-1] = Vtt[,,t-1] + J[,,t-1]%*%(VtT[,,t]-Vtt1[,,t])%*%t.J  # eqn 6.48
      #VtT[,,t-1] = (VtT[,,t-1]+matrix(VtT[,,t-1],m,m,byrow=TRUE))/2     #should not be necessary here
    } #end of the smoother
 
    #define J0
    if(init.state=="x00") { #Shumway and Stoffer treatment of initial conditions; LAM and pi defined for x_0
      if("B" %in% time.varying){
        B = parmat(MLEobj, "B", t=1)$B
        if(m==1) t.B=B else t.B = matrix(B,m,m,byrow=TRUE) 
      }    
      if("Q" %in% time.varying){
        Q = parmat(MLEobj, "Q", t=1)$Q 
        if(m==1){ diag.Q=unname(Q) }else{ diag.Q = takediag(unname(Q)) } 
      }      
      #deal with any 0s on diagonal of Vtt1; these can arise due to 0s in V0, B, + Q
      #0s on diag of Vtt1 will break the Kalman smoother if t>1
      diag.Vtt1 = unname(Vtt1[,,1]); diag.Vtt1=diag.Vtt1[1 + 0:(m - 1)*(m + 1)]   #much faster way to get the diagonal
      if(debugkf!=-1 & any(diag.Vtt1==0) ) {  
        #deal with 0s that are ok if there are corresponding 0s on Q diagonal
        Q0s=identical(which(diag.Q==0),which(diag.Vtt1==0))
        if(!Q0s && (init.state=="x00" || (init.state=="x10" && t>1)) ){
          return(list(ok=FALSE, errors=paste("Stopped in MARSSkfss: soln became unstable when zeros appeared on the diagonal of Vtt1 at t=1.\n") ) )
        }
      }
      if(m==1){ Vinv=pcholinv(matrix(Vtt1[,,1],1,1))  #pcholinv doesn't like vectors
      }else{
        Vinv=pcholinv(Vtt1[,,1])
        Vinv = symm(Vinv)  #to enforce symmetry after chol2inv call
      }
      J0 = V0%*%t.B%*%Vinv  # eqn 6.49 and 1s on diag when Q=0; Here it is t.B[1]
      x0T = x0 + J0%*%(xtT[,1,drop=FALSE]-xtt1[,1,drop=FALSE]);          # eqn 6.47
      V0T = V0 + J0%*%(VtT[,,1]-Vtt1[,,1])%*%t(J0)   # eqn 6.48
      V0T = symm(V0T) #enforce symmetry
    }
    if(init.state=="x10") { #Ghahramani treatment of initial states; LAM and pi defined for x_1
      if(m==1) J0=matrix(J[,,1],1,1) else J0 = J[,,1]
      x0T = xtT[,1,drop=FALSE]
      if(m==1) V0T=matrix(VtT[,,1],1,1) else V0T = VtT[,,1]
    }
    
    #LAG 1 Covariance smoother
    #run another backward recursion to get E[x(t)x(t-1)|y(T)]
    if("Z" %in% time.varying){ Z = parmat(MLEobj, "Z", t=TT)$Z }  #in 6.55, Z[TT] appears
    Zt = Z; Zt[YM[,TT]==0,]=0     #much faster than Mt%*%Z
    if("B" %in% time.varying){ B = parmat(MLEobj, "B", t=TT)$B }  #in 6.55, B[TT] appears
    KT = matrix(Kt[,,TT], m, n); #funny array call to prevent R from restructuring dims
    Vtt1T[,,TT] = (I.m - KT%*%Zt)%*%B%*%Vtt[,,TT-1] #eqn. 6.55 this is Var(x(T)x(T-1)|y(T)); not symmetric
    s = seq(TT,3)
    for (i in 1:(TT-2)) {
       t = s[i]
       if("B" %in% time.varying){ B = parmat(MLEobj, "B", t=t)$B } #in 6.56, B[t] appears
       if(m==1) t.J = J[,,t-2] else t.J = matrix(J[,,t-2],m,m,byrow=TRUE) #faster transpose
       Vtt1T[,,t-1] = Vtt[,,t-1]%*%t.J + J[,,t-1]%*%(Vtt1T[,,t]-B%*%Vtt[,,t-1])%*%t.J   #eqn 6.56
    }
    if(init.state=="x00"){
      if("B" %in% time.varying){ B = parmat(MLEobj, "B", t=2)$B }
      Vtt1T[,,1] = Vtt[,,1]%*%t(J0) + J[,,1]%*%(Vtt1T[,,2]-B%*%Vtt[,,1])%*%t(J0)
    }
    if(init.state=="x10") Vtt1T[,,1] = NA

    ###########################################################
    #Calculate log likelihood, see eqn 6.62
    #Innovations form of the likelihood
    rtn.list = list(xtT = xtT, VtT = VtT, Vtt1T = Vtt1T, x0T = x0T, V0T = V0T, Vtt = Vtt,
            Vtt1 = Vtt1, J=J, J0=J0, Kt=Kt, xtt1 = xtt1, xtt=xtt, Innov=vt, Sigma=Ft)
    loglike = -sum(YM)/2*log(2*base::pi)    #sum(YM) is the number of data points
    for (t in 1:TT) {
      if(n==1) diag.Ft=Ft[,,t] else diag.Ft=unname(Ft[,,t])[1 + 0:(n - 1)*(n + 1)]
      if( any(diag.Ft==0)){
        if(t>1 || (t==1 & init.state=="x00")){
         return(c(rtn.list,list(ok=FALSE, logLik = NaN,
         errors = paste("One of the diagonal elements of Sigma[,,",t,"]=0. That should never happen when t>1 or t=1 and tinitx=0.  \n Are both Q[i,i] and R[i,i] being set to 0?\n",sep=""))))
        }else{ #t=1 so ok. get the det of Ft and deal with 0s that might appear on diag of Ft when t=1 and V0=0 and R=0 and tinitx=1
#2-5-15 This isn't an error.  x10 can be != y[1] when R=0; x11 cannot be.
#          if(any(abs(vt[diag.Ft==0,1])>1E-16)){
#              return(c(rtn.list,list(ok=FALSE, logLik = -Inf, errors = "V0=0, tinitx=1, R=0 and y[1] does not match x0. You shouldn't estimate x0 when R=0.\n")))
#          }
          OmgF1=makediag(1,n)[diag.Ft!=0,,drop=FALSE] #permutation matrix
          #need to remove those y[1] associated with Ft[,,1]==0 that were non-missing
          loglike = loglike + sum(diag.Ft==0 & YM[,1]!=0)/2*log(2*base::pi)   
          
          if(dim(OmgF1)[1]==0){ #no non-zero Ft[,,1]
            detFt=1 #means R and diag(Ft[,,1] all 0; will become 0 when logged
          }else{
            #when R(i,i) is 0 then vt_t(i) will be zero and Sigma[i,i,1] will be 0 if V0=0.
            #OmgF1 makes sure we don't try to take 1/0 
            if(length(OmgF1%*%Ft[,,t]%*%t(OmgF1))==1)
              detFt = OmgF1%*%Ft[,,t]%*%t(OmgF1)
            else detFt = det(OmgF1%*%Ft[,,t]%*%t(OmgF1))
          }
          #get the inv of Ft
          if(n==1){ Ftinv=pcholinv(matrix(Ft[,,t],1,1)) 
          }else{
            Ftinv = pcholinv(Ft[,,t])  #pcholinv deals with 0s on diagonal
            Ftinv = symm(Ftinv)  #to enforce symmetry after chol2inv call
          }
        }
      }else{ #not any diag of Ft==0
        if(n==1){ 
          detFt=Ft[,,t]
          Ftinv=1/Ft[,,t]
        }else{
          detFt=det(Ft[,,t])
          Ftinv = chol2inv(chol(Ft[,,t])) #don't use pcholinv since it is slower
          Ftinv = symm(Ftinv)  #to enforce symmetry after chol2inv call
        }
      }
      if( detFt<0 || !is.finite(log(detFt)) )
          return(c(rtn.list,list(ok=FALSE, logLik=NaN, Sigma=Ft, errors=paste("Stopped in MARSSkfss: log(det(Ft[,,",t,"]))=NA.\n",sep="") ) ) )

      #matrix call here is a transpose
      loglike = loglike - (1/2) %*% matrix(vt[,t],1,n) %*% Ftinv %*% vt[,t,drop=FALSE] - (1/2)*log(detFt)
      loglike = as.vector(loglike)
    }
    if( !is.finite(loglike) ) return(c(rtn.list,list(ok=FALSE, errors=paste("Stopped in MARSSkfss: loglike computed to NA.\n")) ) )

    return(c(rtn.list,list(logLik = loglike, ok=TRUE, errors = msg)))
}

symm = function(x){
t.x = matrix(x,dim(x)[2],dim(x)[1],byrow=TRUE)
x=(x+t.x)/2
x
}
