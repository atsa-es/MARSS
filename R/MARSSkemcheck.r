MARSSkemcheck = function( MLEobj ){
#This checks that the model can be handled by the MARSSkem algorithm
# Most of this is implementing the restrictions in Summary of Requirements for Degenerate Models in derivation
MODELobj=MLEobj$marss
fixed = MODELobj$fixed
free = MODELobj$free
par.dims=attr(MODELobj,"model.dims")
m=par.dims[["x"]][1]; n=par.dims[["y"]][1]
TT=par.dims[["data"]][2]
pseudolim=1E-8
ok=TRUE
msg=NULL

# If TT=2 then kf will break
  if(TT<=2){
    msg=c(msg, "The number of time steps is <=2.\nMore than 2 data points are needed to estimate parameters.\n")
    ok=FALSE
    }

  ############ Check that if fixed, B is within the unit circle
   for(t in 1:par.dims[["B"]][3]){
    ifixed=min(t, dim3(fixed$B))
    if( is.fixed(subFree(free$B, t=ifixed)) ){  
    #parmat needs MODELobj and par list
      if(is.null(MLEobj$par$B)) tmpparB=MLEobj$start$B else tmpparB=MLEobj$par$B
    tmp.MLEobj=list( marss=MODELobj,par=list(B=tmpparB) ) #B is fixed but par might have cols from other times
    parB = parmat(tmp.MLEobj,"B",t=t)$B
    if( !all(abs(Re(eigen(parB,only.values=TRUE)$values))<=1)){ 
      msg=c(msg, " All the eigenvalues of B must be within the unit circle: all(abs(Re(eigen(fixed$B)$values))<=1)\n")
      ok=FALSE
    } }
   } #end for over time to check B
   
   ############ Check that if R has 0s, then the corresponding row of A and Z are fixed
   # See Summary of Requirements for Degenenate Models in EMDerivations.pdf
   el="R"
   #extracts the r=0 rows (or cols)
   diag.rows = 1 + 0:(par.dims[[el]][1] - 1)*(par.dims[[el]][1] + 1)
   Tmax=0
   for(par.test in c("R","Z","A")){
     Tmax = max(Tmax, par.dims[[par.test]])
   }

   if(is.null(MLEobj[["par"]])) MLEobj$par=MLEobj$start
   #Run with MARSSkfss because it has more internal checks and error messages
   MLEobj$kf=MARSSkfss(MLEobj)
   if(!MLEobj$kf$ok){ return(list(ok=FALSE, msg=c(MLEobj$kf$errors, msg))) }
   MLEobj$Ey=MARSShatyt(MLEobj)
   msg.tmp=NULL
   for(i in 1:TT){
      ifixed=min(i,dim3(fixed[[el]])); ifree=min(i,dim3(free[[el]]))
      free.elem.i = subFree(free[[el]], t=ifree)
      fixed.elem.i = subFixed(fixed[[el]], t=ifixed)
      # zeros on the diagonals
      zero.diags = zeroDiags(fixed.elem.i, free.elem.i, par.dims[[el]][1]) 
      if(any(zero.diags)){
        #II0 = makediag(as.numeric(zero.diags))
        II0 = Diagonal(x=as.numeric(zero.diags)) # x is 0s and 1s
        if(i<=Tmax){
          for( par.test in c("Z","A")){
            #II = diag(1,par.dims[[par.test]][2])
            II = Diagonal(1,par.dims[[par.test]][2])
            ifree.par=min(i,dim(free[[par.test]])[3])
            dpart = subFree2D(free[[par.test]],t=ifree.par)
            par.not.fixed = any( ((II %x% II0)%*%dpart)!=0 )
      
            if(par.not.fixed){
              msg=c(msg, paste("t=",i,": For method=kem (EM), if an element of the diagonal of R is 0, the corresponding row of ", par.test, " must be fixed.\n",sep=""))
            }
           } #for par.test
        } #if par needs to be tested; don't need to do over all TT just up to Tmax

        Z=parmat(MLEobj,"Z",t=i)$Z
        A=parmat(MLEobj,"A",t=i)$A
        Z.R0=Z[zero.diags,,drop=FALSE]
        A.R0=A[zero.diags,,drop=FALSE]
        y.R0=MLEobj$Ey$ytT[zero.diags,i,drop=FALSE]
        y.resid = y.R0 - Z.R0%*%MLEobj$kf$xtT[,i,drop=FALSE] - A.R0
        if(any(y.resid>pseudolim)){
          msg.tmp=c(msg.tmp, paste(" Z, A, and y for the R=0 rows at time=",i," do not agree. For the R=0 rows, E(y) must equal Z * E(x) + A.\n",sep=""))
          }
       }#any zero.diags
  } #end for loop over time
  if(!is.null(msg)){
    msg=c(msg, msg.tmp[1:min(10,length(msg.tmp))]) #error msgs could be really long, only print first 10
    ok=FALSE
  }
      

   ############ Check that II_q^{0} is time constant
   # See Summary of Requirements for Degenenate Models in EMDerivations.pdf
   Tmax=0
   for(par.test in c("Q")){
     Tmax = max(Tmax, par.dims[[par.test]])
   }
   II0=list()
   for(i in 1:Tmax){
      for(el in c("Q")){
        ifixed=min(i,dim3(fixed[[el]])); ifree=min(i,dim3(free[[el]])) 
        free.elem.i = subFree(free[[el]], t=ifree)
        fixed.elem.i = subFixed(fixed[[el]], t=ifixed)
        # zeros on the diagonals
        zero.diags = zeroDiags(fixed.elem.i, free.elem.i, par.dims[[el]][1]) 
        #II0Q = makediag(as.numeric(zero.diags))
        #II0Q = unname(II0[[el]])
        II0Q = Diagonal(x=as.numeric(zero.diags))
      }
      if(i==1){
         II0Q.1=II0Q
      }else{
         if(!isTRUE(all.equal(II0Q, II0Q.1))){
            msg=c(msg, paste("t=",i,": The placement of 0 variances in Q must be time constant.\n",sep=""))
            ok=FALSE
         }
      }
   } #end for loop over time


   ############ Check that if R and Q both have 0s, then the U and Bs are appropriately fixed
   # See Summary of Requirements for Degenenate Models in EMDerivations.pdf
   Tmax=0
   for(par.test in c("R","Q","U","B","Z")){
     Tmax = max(Tmax, par.dims[[par.test]][3])
   }
   II0=list()
   for(i in 1:Tmax){
      for(el in c("R","Q")){
        ifixed=min(i,dim3(fixed[[el]])); ifree=min(i,dim3(free[[el]])) 
        free.elem.i = subFree(free[[el]], t=ifree)
        fixed.elem.i = subFixed(fixed[[el]], t=ifixed)
        # zeros on the diagonals
        zero.diags = zeroDiags(fixed.elem.i, free.elem.i, par.dims[[el]][1]) 
        #II0[[el]] = makediag(as.numeric(zero.diags))
        II0[[el]] = Diagonal(x=as.numeric(zero.diags))
      }
      for( el in c("B","U")){
          ifree=min(i,dim3(free[[el]]))
          #I don't know what par$Z will be.  I want to create a par$Z where there is a non-zero value for any potentially non-zero Z's
          tmp.MLEobj=list(marss=MODELobj, par=list(Z=matrix(1,dim2(free$Z),1)))
          tmp.MLEobj$marss$fixed$Z[tmp.MLEobj$marss$fixed$Z!=0]=1
          tmp.MLEobj$marss$free$Z[tmp.MLEobj$marss$free$Z!=0]=1
          if(isMatrix) attr(tmp.MLEobj$marss$free$Z, "free.dims")=attr(MODELobj$free$Z, "free.dims")
          parZ=parmat(tmp.MLEobj,"Z",t=i)$Z
          #II = diag(1,par.dims[[el]][2])
          II = Diagonal(par.dims[[el]][2])
          dpart = subFree2D(free[[el]],t=ifree) #subFree2D ensures that if 2D vec form, then the matrix ret is 2D
          par.not.fixed = any( ((II %x% (II0$R%*%parZ%*%II0$Q))%*%dpart)!=0 )
      
          if(par.not.fixed){
            msg=c(msg, paste("t=",i,": For method=kem (EM), if an element of the diagonal of R & Q is 0, the corresponding row of ", par.test, " must be fixed.\n",sep=""))
            ok=FALSE
          }
        } #for par.test
   } #end for loop over time   

   ############ Check that B^{0} is fixed
   # See Summary of Requirements for Degenenate Models in EMDerivations.pdf
   Tmax=0
   for(par.test in c("B")){
     Tmax = max(Tmax, par.dims[[par.test]][3])
   }
   II0=list()
   for(i in 1:Tmax){
      for(el in c("Q")){
        ifixed=min(i,dim3(fixed[[el]])); ifree=min(i,dim3(free[[el]])) 
        free.elem.i = subFree(free[[el]], t=ifree)
        fixed.elem.i = subFixed(fixed[[el]], t=ifixed)
        # zeros on the diagonals
        zero.diags = zeroDiags(fixed.elem.i, free.elem.i, par.dims[[el]][1]) 
        #II0[[el]] = makediag(as.numeric(zero.diags))
        II0[[el]] = Diagonal(x=as.numeric(zero.diags))
      }
      for( el in c("B")){
          ifree=min(i,dim3(free[[el]]))
          II = Diagonal(par.dims[[el]][2])
          dpart = subFree2D(free[[el]],t=ifree)
          par.not.fixed = any( ((II %x% II0$Q)%*%dpart)!=0 )
      
          if(par.not.fixed){
            msg=c(msg, paste("t=",i,": If an element of the diagonal of Q is 0, the corresponding row and col of ", par.test, " must be fixed.\n",sep=""))
            ok=FALSE
          }
        } #for par.test
   } #end for loop over time   


   ############ Check that if u^{0} or xi^{0} are estimated, B adjacency matrix is time invariant
   # See Summary of Requirements for Degenenate Models in EMDerivations.pdf
   Tmax=0
   for(par.test in c("U")){
     Tmax = max(Tmax, par.dims[[par.test]][3])
   }
   # location of zeros on Q matrix must be time constant so doesn't matter what t is used
   t=1
   II0=list()
   el="Q"
   ifixed=1; ifree=1
   free.elem.i = subFree(free[[el]], t=ifree)
   fixed.elem.i = subFixed(fixed[[el]], t=ifixed)
   # zeros on the diagonals
   zero.diags = zeroDiags(fixed.elem.i, free.elem.i, par.dims[[el]][1]) 
   #II0[[el]] = makediag(as.numeric(zero.diags))
   II0[[el]] = Diagonal(x=as.numeric(zero.diags))
   
   dpart=subFree2D(free$x0,t=1)
   test.adj = any( (II0$Q%*%dpart)!=0 )
   for(i in 1:Tmax){
    dpart=subFree2D(free$U,t=i)
    test.adj = test.adj & any( (II0$Q%*%dpart)!=0 )  #II0$Q required to be time constant above
   }
   
   if(test.adj){ #means x0^{0} or u^{0} being estimated
   Tmax=0
   for(par.test in c("U", "B")){
     Tmax = max(Tmax, par.dims[[par.test]][3])
   }
   for(i in 1:Tmax){
      el="B"
      ifree=min(i,dim3(free[[el]]))
      #I don't know what par$B will be.  I want to create a par$B where there is a non-zero value for any potentially non-zero B's
      tmp.MLEobj=list(marss=MODELobj, par=list(B=matrix(1,dim2(free$B),1)))
      tmp.MLEobj$marss$fixed[[el]][tmp.MLEobj$marss$fixed[[el]]!=0]=1
      tmp.MLEobj$marss$free[[el]][tmp.MLEobj$marss$free[[el]]!=0]=1
      if(isMatrix) attr(tmp.MLEobj$marss$free[[el]], "free.dims")=attr(MODELobj$free[[el]], "free.dims")
      adjB=parmat(tmp.MLEobj,el,t=i)[[el]]
      adjB[adjB!=0]=1; adjB=unname(adjB)
      if(i==1){
         adjB.1=adjB
      }else{
         if(!isTRUE(all.equal(adjB, adjB.1))){
            msg=c(msg, paste("t=",i,": If u^{0} or xi^{0} are estimated, the adjacency matrix specified by B must be time constant.\n",sep=""))
            ok=FALSE
         }
      }
   } #end for loop over time
   } #if the adj test needs to be done
   
############ Check that II_q^{d} and II_q^{is} are time constant
   # See Summary of Requirements for Degenenate Models in EMDerivations.pdf
   Tmax=0
   for(par.test in c("Q")){
     Tmax = max(Tmax, par.dims[[par.test]][3])
   }
   IIz=IIp=OMGz=OMGp=IId=IIis=list()
   for(i in 1:Tmax){
      el="Q"
      ifixed=min(i,dim3(fixed[[el]])); ifree=min(i,dim3(free[[el]])) 
      free.elem.i = subFree(free[[el]], t=ifree)
      fixed.elem.i = subFixed(fixed[[el]], t=ifixed)
      # zeros on the diagonals
      zero.diags = zeroDiags(fixed.elem.i, free.elem.i, par.dims[[el]][1]) 
      # IIz[[el]] = makediag(as.numeric(zero.diags))
      # IIp[[el]] = diag(1,par.dims[[el]][1])-IIz[[el]]
      # OMGz[[el]] = diag(1,par.dims[[el]][1])[diag(IIz[[el]])==1,,drop=FALSE]
      # OMGp[[el]] = diag(1,par.dims[[el]][1])[diag(IIp[[el]])==1,,drop=FALSE]
      IIz[[el]] = Diagonal(x=as.numeric(zero.diags))
      IIp[[el]] = Diagonal(par.dims[[el]][1])-IIz[[el]]
      OMGz[[el]] = Diagonal(par.dims[[el]][1])[diag(IIz[[el]])==1,,drop=FALSE]
      OMGp[[el]] = Diagonal(par.dims[[el]][1])[diag(IIp[[el]])==1,,drop=FALSE]
        
        #I don't know what par$B will be.  I want to create a par$B where there is a non-zero value for any potentially non-zero B's
        tmp.MLEobj=list(marss=MODELobj, par=list(B=matrix(1,dim2(free$B),1)))
        tmp.MLEobj$marss$fixed$B[tmp.MLEobj$marss$fixed$B!=0]=1
        tmp.MLEobj$marss$free$B[tmp.MLEobj$marss$free$B!=0]=1
        if(isMatrix) attr(tmp.MLEobj$marss$free$B, "free.dims")=attr(MODELobj$free$B, "free.dims")
        Adj.mat=parmat(tmp.MLEobj,"B",t=i)$B
        Adj.mat[Adj.mat!=0]=1; Adj.mat=unname(Adj.mat)
        Adj.mat.pow.m = matrix.power(Adj.mat, m) #to find all the linkages

        Q.0.rows.of.Adj.mat = OMGz[[el]]%*%Adj.mat.pow.m
        #These are the columns corresponding to the directly stochastic bit: Q.0.rows.of.Adj.mat%*%t(OMGp$Q)
        if(dim(OMGp[[el]])[1]!=0){
        tmp=Q.0.rows.of.Adj.mat
        #which rows of the + columns are all zero
        if(dim(Q.0.rows.of.Adj.mat)[1]!=0) tmp=apply(Q.0.rows.of.Adj.mat%*%t(OMGp[[el]])==0,1,all) 
        tmp=t(OMGz[[el]])%*%matrix(as.numeric(tmp),ncol=1) #expand back outl 1s where the deterministic x's are
        #IId[[el]]=makediag(tmp); IId[[el]]=unname(IId[[el]])
        #IIis[[el]]=diag(1,m)-IId[[el]]-IIp[[el]]; IIis[[el]]=unname(IIis[[el]])
        IId[[el]]=Diagonal(x=tmp[,1]) # [,1] to vectorize
        IIis[[el]]=Diagonal(m)-IId[[el]]-IIp[[el]]
        }else{ #Q all 0
         # IId[[el]]=diag(1,par.dims[[el]][1])
         # IIis[[el]]=diag(0,par.dims[[el]][1])
         IId[[el]]=Diagonal(par.dims[[el]][1])
         IIis[[el]]=Diagonal(par.dims[[el]][1],0)
        }
        if(i==1){
          IId.1=IId[[el]]
          IIis.1=IIis[[el]]
        }else{
          if(!isTRUE(all.equal(IId[[el]], IId.1))){
            msg=c(msg, paste("t=",i,": The location of the deterministic x's must be time constant.\n",sep=""))
            ok=FALSE
          }
          if(!isTRUE(all.equal(IIis[[el]], IIis.1))){
            msg=c(msg, paste("t=",i,": The location of the indirectly stochastic x's must be time constant.\n",sep=""))
            ok=FALSE
          }
       }
     } #end for loop over time

return(list(ok=ok, msg=msg))
}
