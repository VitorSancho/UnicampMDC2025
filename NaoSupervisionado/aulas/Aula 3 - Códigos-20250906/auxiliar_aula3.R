####################
# AUX functions rNMF
####################

Nnls <- function(M = M, C = C){
  fun1 = function(x){ return(nnls(M, x)$x)}
  return(apply(C, 2, fun1))
}



is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

checkargs <- function(A, k, alpha, beta, maxit, tol, gamma, ini.W, ini.zeta, my.seed, variation, quiet, nreg, p, n){
  if(!is.matrix(A)) {stop("The input A is not a matrix. Consider as.matrix(A)?")}
  if(!is.numeric(A)) {stop("The input A is not a numeric matrix.")}
  if(!all(A >= 0)) {stop("Not all entries are non-negative.")}
  if(length(k) != 1 | !is.wholenumber(k) | k <= 0) {stop("Something is wrong about k. It should be a positive integer.")}
  if(length(alpha) != 1 | !is.numeric(alpha) | alpha < 0) {stop("Argument 'alpha' must be non-negative.")}
  if(length(beta) != 1 | !is.numeric(beta) | beta < 0) {stop("Argument 'beta' must be non-negative.")}
  if(length(maxit) != 1 | !is.wholenumber(maxit) | maxit <= 0) {stop("Something is wrong about maxit. It should be a positive integer.")}
  if(length(tol) != 1 | !is.numeric(tol) | tol < 0 | tol > 1) {stop("tol must be a number in [0,1)")}
  if(gamma != FALSE){
    if(length(gamma) != 1 | !is.numeric(gamma) | gamma < 0 | gamma > 1) {stop("Argument 'gamma' must be a number in [0,1), or 'FALSE'.")}
  }
  if(!is.null(ini.W)){
    if(!is.matrix(ini.W) | !is.numeric(ini.W) | !all(ini.W >= 0) | nrow(ini.W) != nrow(A) | ncol(ini.W) != k){stop("ini.W must be a p by k non-negative numeric matrix.")}
  }
  if(!is.null(ini.zeta)){
    if(!is.matrix(ini.zeta) | !is.logical(ini.zeta) | nrow(ini.zeta) != nrow(A) | ncol(ini.zeta) != ncol(A)) {stop("ini.zeta must be a logical matrix of the same size as A.")}
    if(sum(c(!ini.zeta)) > round(gamma * p * n)) {stop("ini.zeta contains too many FALSES (outliers); Increase trimming percentage?")}
  }
  if(!is.null(my.seed)){
    if(length(my.seed) != 1 | !is.numeric(my.seed) | !is.wholenumber(my.seed) | my.seed <= 0) {stop("Something is wrong about my.seed. It should be a positive integer.")}
  }
  if(length(variation) != 1 | !is.character(variation) | !(variation %in% c('col', 'row', 'cell', 'smooth', 'row', 'rowsmooth'))) {stop("Argument variation must be one of the following strings 'col', 'row', 'cell', 'smooth', 'row' or 'rowsmooth'")}
  if(length(quiet) != 1 | !is.logical(quiet)) {stop("The argument 'quiet' must be logical.")}
  if(length(nreg) != 1 | !is.numeric(nreg) | !is.wholenumber(nreg) | nreg <= 0) {stop("The argument 'nreg' must be a positive ineger.")}
}

rnmf <- function(X, k = 5, alpha = 0, beta = 0, maxit = 20, tol = 0.001,
                 gamma = FALSE, ini.W = NULL, ini.zeta = NULL, my.seed = NULL, 
                 variation = "cell", quiet = FALSE, nreg = 1, showprogress = TRUE)
{
  tic <- proc.time() # Start a clock.
  p <- nrow(X)
  n <- ncol(X)
  ## Checks if the arguments are valid. Display an error message if not.
  checkargs(X, k, alpha, beta, maxit, tol, gamma, ini.W, ini.zeta, my.seed, 
            variation, quiet, nreg, p, n)
  
  ## Create a data frame of 4 columns: value = entries of X; 
  ## (x,y) = coordinates; outs = is it an outlier?
  X.f <- data.frame(value = as.vector(X), x = rep(1 : p, n),
                    y = rep(1 : n, each = p), outs = FALSE) 
  
  ## to.trim (list) stores entries trimmed in each iteration.
  if(gamma > 0) to.trim <- vector("list")
  else to.trim <- NULL
  ## 'obj' (vector) stores values of the objective function in each iteration.
  obj <- 1
  
  ## Initialize W
  if(missing(ini.W)){
    if(missing(my.seed)){
      W <- initM(large = max(X), nrow = p, ncol = k, small = 0)
    }else{ 
      W <- initM(large = max(X), nrow = p, ncol = k, small = 0,
                 my.seed = my.seed)
    }
  }else{
    W <- ini.W
  }
  ## Creates a progress bar.
  if(showprogress) pb <- txtProgressBar(min = 0, max = maxit, style = 3)
  
  if(variation == "cell" & gamma > 0){
    ##--------------------------
    ## Cell trimming
    ##--------------------------
    ## Initialize zeta (logic matrix the same size as X)
    if(!is.null(ini.zeta)){
      if(sum(ini.zeta) < round(gamma * p * n)) {
        warning("The percentage of FALSES (outliers) in ini.zeta 
                        is smaller than 'gamma'; Other outliers are randomly picked.")
        more.outliers <- sample(which(!ini.zeta), 
                                round(gamma * p * n) - sum(ini.zeta))
        ini.zeta[more.outliers] <- TRUE
      }else{}
      zeta <- ini.zeta
    }else{
      zeta <- matrix(FALSE, nrow = p, ncol = n)
      zeta[sample(1 : (p * n), round(gamma * p * n))] <- TRUE ## randomize zeta
    }
    
    ## Start iterations.
    for(i in 1 : maxit){
      if(showprogress) setTxtProgressBar(pb, i) ## update the progress bar
      ##------Stage 1------##
      ## Given W, fit H
      H <- Nnls.trimH(W, X, !zeta, beta, k, n)
      for(j in 1:nreg){
        ## Find residuals
        R <- abs(X - W %*% H)
        to.trim[[i]] <- order(R, decreasing = TRUE)[1 : round(gamma * n * p)]
        ## to.trim[[i]] = which(rank(R) > round((1 - gamma) * n * p))
        ## Update zeta
        zeta <- matrix(FALSE, nrow = p, ncol = n)
        zeta[to.trim[[i]]] <- TRUE
        ## Refit H
        H <- Nnls.trimH(W, X, !zeta, beta, k, n)
      }
      
      ##------Stage 2------##
      ## Given H, Fit W
      W <- Nnls.trimW(H, X, !zeta, alpha, k, p)
      for(j in 1:nreg) {
        ## Find residuals
        R <- abs(X - W %*% H)
        ##to.trim[[i]] = which(rank(R) > round((1 - gamma) * n * p))
        to.trim[[i]] <- order(R, decreasing=TRUE)[1 : round(gamma * n * p)]
        ## Update zeta
        zeta <- matrix(FALSE, nrow = p, ncol = n)
        zeta[to.trim[[i]]] <- TRUE
        ## Refit W
        W <- Nnls.trimW(H, X, !zeta, alpha, k, p)
        J <- nmlz(W) # Find the normalizing matrix of W
        W <- W %*% J; H = solve(J) %*% H
      }
      obj[i] <- l2((X - W %*% H)[!zeta]) + l2(W) + sum(colSums(abs(H))^2)
      
      ## Check convergence
      #             if(i > 1){
      #                 if(all(to.trim[[i]] == to.trim[[i - 1]]) &
      #                    sum((W - W.prev)^2) / sum(W.prev^2) < tol) break
      #             }else{}
      if(i > 1){
        if(sum((W - W.prev)^2) / sum(W.prev^2) < tol) break
      }else{}
      W.prev <- W
    }
  }else{
    ##--------------------------
    ## Row, Col and Smooth trimming
    ##--------------------------  
    for(i in 1 : maxit){
      flag <- FALSE
      if(showprogress) setTxtProgressBar(pb, i)  # update the progress bar
      ##------Stage 1------##
      #Given W, find H in ||MH - C||
      C <- rbind(X, matrix(0,1,n))
      M <- rbind(W, sqrt(beta) * matrix(1, 1, k))
      H <- Nnls(M,C)
      if(gamma > 0){ ## Trim
        R <- (X - W %*% H)^2
        if(variation == "col"){
          # No change here
        }else if(variation == "smooth"){
          to.trim[[i]] <- which(rank(abs(R)) > round((1 - gamma) * n * p))
          X.s <- matrix(smoothing2(X.f, to.trim[[i]], p,n, frame = TRUE)$value,p,n)
          ## Fit H with trimmed X and W. 
          C.s <- rbind(X.s, matrix(0, nrow = 1, ncol = n))
          H.prev <- H
          H <- Nnls(M,C.s)
        }else if(variation == "row"){
          rsqr <- apply(R, 1, l2) ## Find SS of residuals of rows
          to.trim[[i]] <- which(rank(rsqr) > round((1 - gamma) * p))
          X.trim <- X[-to.trim[[i]],]
          W.trim <- W[-to.trim[[i]],]
          ## Fit H with trimmed X and W. 
          M.trim <- rbind(W.trim, sqrt(beta) * matrix(1, nrow = 1, ncol = k))
          C.trim <- rbind(X.trim, matrix(0, nrow = 1, ncol = n))
          H.prev <- H
          H <- Nnls(M.trim, C.trim)
        }else{
          stop("Wrong mode. Try one of the following: 'cell', 'col', 'row', 'smooth'")
        }
      }else{}
      
      ##------Stage 2------##
      ##Given H, find W in min||MW^T - C||.
      C <- rbind(t(X), matrix(0, nrow = k, ncol = p))
      M <- rbind(t(H), sqrt(alpha) * diag(k))
      W <- t(Nnls(M, C))
      if(gamma > 0){
        R <- (X - W %*% H)^2
        if(variation == "col"){
          rsqr <- apply(R, 2, l2) ## Find SS of residuals of columns
          to.trim[[i]] <- which(rank(rsqr) > round((1 - gamma) * n))
          X.trim <- X[,-to.trim[[i]]]
          H.trim <- H[,-to.trim[[i]]]
          ## Fit W with trimmed X and H. 
          M.trim <- rbind(t(H.trim), sqrt(alpha) * diag(k))
          C.trim <- rbind(t(X.trim), matrix(0, nrow = k, ncol = p))
          W <- t(Nnls(M.trim,C.trim))
        }else if(variation == "smooth"){
          R <- (X - W %*% H)^2
          to.trim[[i]] <- which(rank(abs(R)) > round((1 - gamma) * n * p))
          X.s <- matrix(smoothing2(X.f, to.trim[[i]], p, n, frame = TRUE)$value,p,n)
          ## Fit W with trimmed X and H. 
          C.s <- rbind(t(X.s), matrix(0, nrow = k, ncol = p))
          W <- t(Nnls(M,C.s))
        }else if(variation == "row"){
          ## No change
        }else{
          stop("Wrong mode. Try one of the following: 'cell', 'col', 'row', 'smooth'")
        }
      }
      ##J = nmlz(W) # Find the normalizing matrix of W
      ##W = W %*% J; H = solve(J) %*% H
      
      ## Convergence?
      if(i > 1){
        if(setequal(to.trim[[i]], to.trim[[i-1]])){
          if(sum((W - W.prev) ^ 2) / sum(W.prev ^ 2) < tol) break
        }
      }else{}
      W.prev <- W
    }
  }
  ## Close the link to the progress bar. 
  if(showprogress) close(pb)
  fit <- W %*% H
  
  ## Print report
  if(!quiet){
    if(!gamma){
      cat("Done. Time used:","\n")
      print(proc.time() - tic)
      cat("No trimming.\n",
          "Input matrix dimension: ",p," by ",n, "\n",
          "Left matrix: ", p," by ", k,". Right matrix: ", k," by ", n,"\n",
          "alpha = ", alpha,". beta = ",beta,"\n",
          "Number of max iterations = ",maxit,"\n",
          "Number of iterations = ", i, "\n", sep = ""
      )
    }else{
      cat("Done. Time used: ","\n")
      print(proc.time() - tic)
      cat("\n Trimming mode = \"", variation, "\". Proportion trimmed: ",gamma, "\n",
          "Input matrix dimension: ",p," by ",n, "\n",
          "Left matrix: ",p," by ",k,". Right matrix: ",k," by ",n,"\n",
          "alpha = ",alpha,". beta = ",beta,"\n",
          "Number of max iterations = ",maxit,"\n",
          "Number of iterations = ",i,"\n", sep = ""
      )
    }
  }
  return(invisible(list(W = W, H = H, fit = fit,
                        trimmed = to.trim, niter = i)))
}

initM = function(large, nrow, ncol, small = 0, my.seed = NULL){
  if(!missing(my.seed)) set.seed(my.seed)
  M = matrix(runif(nrow * ncol, small, large), nrow = nrow, ncol = ncol)
  return(M)
}



####################
# AUX functions lle
####################
lle <-
  function(X,m,k,reg=2,ss=FALSE,p=0.5,id=FALSE,nnk=TRUE,eps=1,iLLE=FALSE,v=0.99) {
    
    #vector of chosen data (subset selection)
    choise <- c() #remains empty if there's no subset selection
    
    #find neighbours
    if( iLLE ) cat("finding neighbours using iLLE\n") else cat("finding neighbours\n") 
    if( nnk ) nns <- find_nn_k(X,k,iLLE)	else nns <- find_nn_eps(X,eps)
    
    #calculate weights
    cat("calculating weights\n")
    res_wgts <- find_weights(nns,X,m,reg,ss,p,id,v)
    
    #if subset selection, the neighbours and weights have to be 
    #calculated again since the dataset changed
    if( ss ){
      
      #use reduced dataset
      X <- res_wgts$X
      choise <- res_wgts$choise
      
      cat("finding neighbours again\n")
      if( nnk ) nns <- find_nn_k(X,k,iLLE)	else nns <- find_nn_eps(X,eps)
      
      cat("calculating weights again\n")
      res_wgts <- find_weights(nns,X,m,reg,FALSE,p,id,v)
    }
    
    #extracting data (weights, intrnsic dim)
    wgts <- res_wgts$wgts
    id <- res_wgts$id
    
    #compute coordinates
    cat("computing coordinates\n")
    Y <- find_coords(wgts,nns,N=dim(X)[1],n=dim(X)[2],m)
    
    return( list(Y=Y,X=X,choise=choise,id=id) )
  }

find_nn_k <-
  function(X,k,iLLE=FALSE){
    
    #calculate distance-matrix
    nns <- as.matrix(dist(X))
    
    #get ranks of all entries
    nns <- t(apply(nns,1,rank))
    
    #choose the k+1 largest entries without the first (the data point itself)
    nns <- ( nns<=k+1 & nns>1 )
    
    #optional: improved LLE
    if( iLLE ){
      N <- dim(X)[1]
      n <- dim(X)[2]
      nns2 <- nns
      nns <- as.matrix(dist(X))
      for( i in 1:N){
        if( i%%100==0 ) cat(i,"von",N,"\n")
        for( j in 1:N){
          Mi <- sqrt( sum( rowSums( matrix( c(t(X[i,])) - c(t(X[nns2[i,],])), nrow=k, ncol=n, byrow=TRUE)^2 )^2 ) )
          Mj <- sqrt( sum( rowSums( matrix( c(t(X[j,])) - c(t(X[nns2[j,],])), nrow=k, ncol=n, byrow=TRUE)^2 )^2 ) )
          nns[i,j] <- nns[i,j]/(sqrt(Mi*Mj)*k) 
        }
      }
      nns <- t(apply(nns,1,rank))
      nns <- ( nns<=k+1 & nns>1 )
    }
    
    return (nns)
  }

calc_k <- 
  function(X,m,kmin=1,kmax=20,plotres=TRUE,parallel=FALSE,cpus=2,iLLE=FALSE){
    
    #set parameters
    N <- dim(X)[1]
    if( kmax>=N ) kmax <- N - 1 #more neighbourse than points doesnt make sense
    if( .Platform$OS.type=="windows" ) dev <- "nul" else dev <- "/dev/null"	
    
    #set up parallel computation
    if( parallel==TRUE) sfInit( parallel=TRUE, cpus=cpus ) else sfInit( parallel=FALSE ) 
    #require lle for every node
    options("warn"=-1)
    sfLibrary( lle )
    options("warn"=0)
    
    
    perform_calc <- function( k, X, m, iLLE=FALSE ){
      
      N <- dim(X)[1]
      
      #perform LLE
      sink( dev ) #surpress output
      Y <- lle(X,m,k,2,0,iLLE=iLLE)$Y
      sink()
      
      #distance matrix of original data
      Dx <- as.matrix(dist(X))
      
      #distance matrix of embedded data
      Dy <- as.matrix(dist(Y))
      
      #calculate correlation between original and embedded data for every data point
      rho <- c()
      for( i in 1:N ) rho <- c( rho, cor(Dx[i,],Dy[i,]) ) 
      
      #rho cant be cumulated, rho^2 can
      return( mean(1-rho^2) )
    }
    
    rho <- invisible( sfLapply( kmin:kmax, perform_calc, X, m, iLLE ) )
    rho <- unclass( unlist( rho ) )
    sfStop()	
    
    res <- data.frame( k=c(kmin:kmax), rho=rho )
    
    if( plotres ){
      par( mar=c(5,5,4,2)+0.1 )
      plot( res$k, res$rho, type="b", xlab="k", ylab=expression(1-rho^2), main="" )
      abline(h=min(res$rho,na.rm=TRUE),col="red")
      grid()
    } else cat( "best k:",head(res$k[order(res$rho)],3), "\n\n" )	
    return ( res )
  }

find_coords <-
  function(wgts,nns,N,n,m)
  {
    W <- wgts
    M <- t(diag(1,N)-W)%*%(diag(1,N)-W)
    
    #calculate the eigenvalues and -vectors of M
    e <- eigen(M)
    
    #choose the eigenvectors belonging to the 
    #m smallest not-null eigenvalues
    #and re-scale data
    Y <- e$vectors[,c((N-m):(N-1))]*sqrt(N)
    
    return(Y)
  }

find_nn_eps <-
  function(X,eps){
    
    #calculate distance matrix
    nns <- as.matrix(dist(X))
    
    #choose neighbours using eps environment
    nns <- nns<eps
    diag(nns) <- FALSE
    
    return(nns)
  }

find_weights <-
  function(nns,X,m,reg=2,ss=FALSE,p=0.5,id=FALSE,v=0.99)
  {
    #get dimensions
    N <- dim(X)[1]
    n <- dim(X)[2]
    
    #matrix of weights
    wgts <- 0*matrix(0,N,N)
    #subset selection vector
    s <- c()
    #intrinsic dim
    intr_dim <- c()
    
    for (i in (1:N)){
      #number of neighbours
      k <- sum(nns[i,])
      
      #no neighbours (find_nn_k(k=0) or eps-neighbourhood)
      if( k==0 ) next
      
      # calculate the differences  between xi and its neighbours
      Z <- matrix(c(t(X)) - c(t(X[i,])), nrow=nrow(X), byrow = TRUE)
      Z <- matrix( Z[nns[i,],], ncol=n, nrow=k )		
      
      #gram-matrix
      G <- Z%*%t(Z)
      
      #regularisation
      delta <- 0.1
      #calculate eigenvalues of G
      e <- eigen(G, symmetric=TRUE, only.values=TRUE)$values
      #skip if all EV are null
      if( all( e==0 ) ) next
      
      #choose regularisation method
      #see documentation
      if( reg==1 ){ 
        r <- delta*sum(head( e, n-m ))/(n-m)
      } else if( reg==2 ){
        r <- delta^2/k*sum(diag(G))
      } else r <- 3*10^-3
      
      #calculate intrinsic dimension
      if( id ){	
        tmp <- 1
        while( sum(e[1:tmp])/sum(e) <= v ) tmp <- tmp + 1
        intr_dim <- c( intr_dim, tmp )
      }
      
      #basic value for subset selection
      s <- c( s, sum(head( e, n-m ))/(n-m) )
      
      #use regularisation if more neighbourse than dimensions!
      if( k>n ) alpha <- r else alpha <- 0
      
      #regularisation
      G <- G + alpha*diag(1,k)
      
      #calculate weights
      #using pseudoinverse ginv(A): works better for bad conditioned systems
      if( k >= 2) wgts[i,nns[i,]] <- t(ginv(G)%*%rep(1,k)) else wgts[i] <- G
      wgts[i,] <- wgts[i,]/sum(wgts[i,])
    }
    
    #subset selection
    #see documentation
    if( ss ){
      s <- s/sum(s)
      Fs <- ecdf(s)
      choise <- sample( 1:N, round(p*N), replace=FALSE, prob=Fs( seq(0,max(s),length=N) ) )
      X <- X[choise,]
    } else choise <- 0
    
    #print intrinsic dimension
    if( id ) cat("intrinsic dim: mean=",mean(intr_dim),", mode=",names(sort(-table(intr_dim)))[1],"\n",sep="")
    #print(intr_dim)
    
    return(list(X=X, wgts=wgts, choise=choise,id=intr_dim))
  }

lle_rectangular <-
  function(N=40,k=5,v=0.9){
    
    #dimension
    n <- 500 
    t <- seq( 0, 1, length=n ) 
    
    #function to generate rectangular signal
    rects <- function(t, width, distance ){
      data <- t*0
      t0 <- t[1]
      t1 <- head(t,1)
      data[t>t0+distance & t<t0+distance+width] <- rnorm(1,1,0.01)
      return( data )
    }
    
    #set parameters
    width_range <- seq( 0.02, 0.35, length=N) #max(dist)+max(width) < 1 (see documentation)
    dist_range <- seq( 0.05, 0.6, length=N )
    
    #generate data
    X <- matrix( 0, ncol=n, nrow=N*N )
    for( i in 1:N ){
      for( j in 1:N ){
        w <- width_range[i]
        d <- dist_range[j]
        X[(i-1)*N+j,] <- rects( t, w, d )
      }
    }
    
    
    #perform lle
    m <- 2
    ss <- FALSE
    p <- 0.8
    reg <- 2
    id <- TRUE
    iLLE <- FALSE
    
    res <- lle(X=X,m=m,k=k,reg=reg,ss=ss,p=p,id=id,iLLE=iLLE,v=v)
    
    Y <- res$Y
    X <- res$X
    choise <- res$choise
    
    #plot
    plot_lle(Y,X,0,col=3,inter=TRUE)
  }

lle_scurve <-
  function(N=800,k=12,ss=FALSE,p=0.5,reg=2,iLLE=FALSE,v=0.8){
    
    #set dimensions
    n <- 3
    m <- 2
    
    #generate S-curve with noise
    X <- 0*matrix(rnorm(n*N),nrow=N)
    angle <- pi*(1.5*runif(N/2)-1)
    height <- 5*runif(N)
    sd <- 0.03
    X[,1] <- c(cos(angle),-cos(angle)) + rnorm(1,0,sd)
    X[,2] <- height + rnorm(1,0,sd)
    X[,3] <- c( sin(angle),2-sin(angle) ) + rnorm(1,0,sd)
    
    #lle
    res <- lle(X,m,k,reg=reg,ss=ss,p=p,iLLE=iLLE)
    Y <- res$Y
    X <- res$X
    choise <- res$choise
    
    #plot
    col <- c(angle,angle)
    if( ss==1 ) col <- col[choise]
    col <- ((col-mean(col))/min(col)+1)*80
    plot_lle(Y,X,0,col,"",10)
  }

lle_sound <-
  function(t=500,dt=20,k=25,reg=2,ss=FALSE,p=0.5,id=TRUE){
    
    #generate data from wavefile
    #result: dataset with 'dt' dimension and approx '(f/dt-t/dt)' samples
    i <- 1 
    invisible( data( lle_wave ) )
    f <- length(lle_wave)
    X <- t( matrix( lle_wave[1:t] ) )
    while( i+t+dt < f ){
      i <- i + dt
      X <- rbind( X, lle_wave[i:(i+t-1)] )
    }
    
    #get/set parameters
    n <- dim(X)[2]
    N <- dim(X)[1]
    m <- 8
    
    #perform lle
    res <- lle(X=X,m=m,k=k,reg=reg,ss=ss,p=p,id=id)
    
    Y <- res$Y
    X <- res$X
    choise <- res$choise
    
    #plot
    plot( res$id, type="l", main="found intrinsic dimension")
    lines( smooth.spline( res$id, df=5), col="red", lwd=2)
    
  }

lle_spiral <- 
  function(){
    
    #set dimensions/parameters
    n <- 3
    m <- 1
    N <- 600
    k <- 5
    reg <- 2
    ss <- FALSE
    p <- 0.5
    iLLE <- FALSE
    v <- 0.8
    
    #generate elementary vectors
    v1 <- c(1,rep(0, each=n-1)) 
    v2 <- c(0,1,rep(0, each=n-2))
    v3 <- c(0,0,1,rep(0,each=n-3))
    dt <- 0.005	
    
    #prepare rotation matrices
    phix <- 0/180*pi #rotation x
    phiy <- 0/180*pi #rotation y
    phiz <- 0/180*pi #rotation z
    Rx <- matrix( c( 1, 0, 0, 0, cos(phix), -sin(phix), 0, sin(phix), cos(phix) ), 3 )
    Ry <- matrix( c( cos(phiy), 0, sin(phiy), 0, 1, 0, -sin(phiy), 0, cos(phiy) ), 3 )
    Rz <- matrix( c( cos(phiz), -sin(phiz), 0, sin(phiz), cos(phiz), 0, 0, 0, 1 ), 3 )
    
    #prepare splitscreen-plot
    plotdim <- c(3,3)
    split.screen(plotdim)
    scr <- 1
    
    for( T in seq(16,1.5,length=prod(plotdim)) ){
      
      #set parameters
      om <- 16*pi/T
      X <- 0*matrix(rnorm(n*N),nrow=N)
      
      #generate spiral data
      for( i in (1:N)){
        t=dt*i
        r=t
        X[i,]= r*v1*cos(om*t)+ r*v2*sin(om*t) + v3*t^3/30
      }
      
      #rotation
      X <- t(Rx%*%t(X))
      X <- t(Ry%*%t(X))
      X <- t(Rz%*%t(X))
      
      #perform lle
      cat( "------------\nT=",T,"\n")
      res <- lle(X,m,k,reg,ss=ss,p=p,id=TRUE,iLLE=iLLE)
      
      #plot
      if( (scr+2)%%3==0 ) zt=NULL else zt="" 
      if( scr>6 ) xt=NULL else xt=""
      if( scr%%3==0 ) yt=NULL else yt=""
      screen(scr)
      scatterplot3d(X,cex.symbols=0.5, mar=rep(.8,4), angle=10, xlab="", ylab="", zlab="", 
                    type="p", lwd=2, highlight.3d=TRUE, cex.axis=0.8,
                    main="", xlim=c(-4,4), ylim=c(-4,4),
                    x.ticklabs=xt, y.ticklabs=yt, z.ticklabs=zt )
      if( scr < 8) text(-1.2,2,"m=1") else text(-1.2,2,"m=2") 
      
      scr <- scr + 1
    }
  }

lle_swissrole <-
  function(N=1500,k=10,ss=FALSE,p=0.5,reg=2,iLLE=FALSE,v=0.8){
    
    #set dimensions
    n <- 3
    m <- 2
    
    #generate swiss roll data
    X <- 0*matrix(0,N,n)
    tt <- (3*pi/2)*(1+2*runif(N))
    height <- 21*runif(N)
    X[,1] <- tt*cos(tt) 
    X[,2] <- height 
    X[,3] <- tt*sin(tt)
    
    #perform lle
    res <- lle(X,m,k,reg=reg,ss=ss,p=p,iLLE=iLLE,v=v)
    
    Y <- res$Y
    X <- res$X
    choise <- res$choise
    
    #plot
    if( ss==0 ) col <- (tt-min(tt)) else col <- (tt[choise]-min(tt[choise])) 
    col <- col/max(col)*200
    plot_lle(Y,X,0,col)
    
  }

plot_lle <-
  function(Y,X,print=FALSE,col=3,name=as.numeric(Sys.time()),angle=60,inter=FALSE){
    
    #interactive plot
    if( inter ){
      require(rgl)
      plot3d( Y, col=col)
      return( NULL )
    }
    
    #get dimensions
    require( scatterplot3d )
    N <- dim(X)[1]
    n <- dim(X)[2]
    m <- dim(Y)[2]
    if( is.null(m) ) m <- 1
    
    #set plot parameters
    pch <- 19
    cex <- 1
    
    #set colour palette of 600 rainbow colours
    palette( rainbow(600) )
    
    #set colour vector
    if( length(col) == 1 & is.numeric( col ) ) col <- seq(0,1,length=N)*100*col
    
    
    #print direction of output file (if used)
    if( print ){
      cat( "graphics saved to",getwd(),"\n" )
      pdf( file=paste(name, ".pdf", sep="" ), width=800, height=480 )
    } else dev.new()
    
    #generate splitscreen
    #when outputing to file, splitscreen generates warning, which is harmless
    options("warn"=-1)
    if( n<4 ) split.screen(c(1,2))
    options("warn"=0)	
    
    #plot original data (depending on its dimension)
    if( n<4 ) screen(1)
    if( n==2 ){ 
      plot( X[,1], X[,2], col=col, main=paste("raw data n=",n,", N=",N,sep=""), pch=pch, cex=cex )
    } else if( n==3 ){
      scatterplot3d( X[,1], X[,2], X[,3], color=col, main=paste("raw data n=",n,", N=",N,sep=""), pch=pch, cex.symbols=cex, angle=angle )
    } else {
      plot( c(1,3),c(1,3), type="n" )
      text(2,2,"data with dimension n>3\n cannot be displayed", cex=1)
    }
    
    #plot embedded data (depending on its dimension)
    if( n<4 ) screen(2)
    if( m==1 ){
      plot( Y, col=col, main=paste("embedded data m=",m,", N=",N,sep=""), pch=pch, cex=cex )
    } else if( m==2 ){
      plot( Y[,2], Y[,1], col=col, main=paste("embedded data m=",m,", N=",N,sep=""), pch=pch, cex=cex )
    } else if( m==3 ){
      scatterplot3d( Y[,3], Y[,2], Y[,1], color=col, main=paste("embedded data m=",m,", N=",N,sep=""), pch=pch, cex.symbols=cex, angle=angle )
    } else {
      plot( c(1,3),c(1,3), type="n" )
      text(2,2,"data with dimension n>3\n cannot be displayed", cex=1)
    }
    
    if( print==1 ) dev.off()	
  }

