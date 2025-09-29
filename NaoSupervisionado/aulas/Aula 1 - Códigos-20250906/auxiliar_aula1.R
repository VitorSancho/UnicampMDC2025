#### DMwR functions



SMOTE <- function(form,data,
                  perc.over=200,k=5,
                  perc.under=200,
                  learner=NULL,...
)
  
{
  
  # the column where the target variable is
  tgt <- which(names(data) == as.character(form[[2]]))
  minCl <- levels(data[,tgt])[which.min(table(data[,tgt]))]
  
  # get the cases of the minority class
  minExs <- which(data[,tgt] == minCl)
  
  # generate synthetic cases from these minExs
  if (tgt < ncol(data)) {
    cols <- 1:ncol(data)
    cols[c(tgt,ncol(data))] <- cols[c(ncol(data),tgt)]
    data <-  data[,cols]
  }
  newExs <- smote.exs(data[minExs,],ncol(data),perc.over,k)
  if (tgt < ncol(data)) {
    newExs <- newExs[,cols]
    data <- data[,cols]
  }
  
  # get the undersample of the "majority class" examples
  selMaj <- sample((1:NROW(data))[-minExs],
                   as.integer((perc.under/100)*nrow(newExs)),
                   replace=T)
  
  # the final data set (the undersample+the rare cases+the smoted exs)
  newdataset <- rbind(data[selMaj,],data[minExs,],newExs)
  
  # learn a model if required
  if (is.null(learner)) return(newdataset)
  else do.call(learner,list(form,newdataset,...))
}

smote.exs <- function(data,tgt,N,k)
  # INPUTS:
  # data are the rare cases (the minority "class" cases)
  # tgt is the name of the target variable
  # N is the percentage of over-sampling to carry out;
  # and k is the number of nearest neighours to use for the generation
  # OUTPUTS:
  # The result of the function is a (N/100)*T set of generated
  # examples with rare values on the target
{
  nomatr <- c()
  T <- matrix(nrow=dim(data)[1],ncol=dim(data)[2]-1)
  for(col in seq.int(dim(T)[2]))
    if (class(data[,col]) %in% c('factor','character')) {
      T[,col] <- as.integer(data[,col])
      nomatr <- c(nomatr,col)
    } else T[,col] <- data[,col]
  
  if (N < 100) { # only a percentage of the T cases will be SMOTEd
    nT <- NROW(T)
    idx <- sample(1:nT,as.integer((N/100)*nT))
    T <- T[idx,]
    N <- 100
  }
  
  p <- dim(T)[2]
  nT <- dim(T)[1]
  
  ranges <- apply(T,2,max)-apply(T,2,min)
  
  nexs <-  as.integer(N/100) # this is the number of artificial exs generated
  # for each member of T
  new <- matrix(nrow=nexs*nT,ncol=p)    # the new cases
  
  for(i in 1:nT) {
    
    # the k NNs of case T[i,]
    xd <- scale(T,T[i,],ranges)
    for(a in nomatr) xd[,a] <- xd[,a]==0
    dd <- drop(xd^2 %*% rep(1, ncol(xd)))
    kNNs <- order(dd)[2:(k+1)]
    
    for(n in 1:nexs) {
      # select randomly one of the k NNs
      neig <- sample(1:k,1)
      
      ex <- vector(length=ncol(T))
      
      # the attribute values of the generated case
      difs <- T[kNNs[neig],]-T[i,]
      new[(i-1)*nexs+n,] <- T[i,]+runif(1)*difs
      for(a in nomatr)
        new[(i-1)*nexs+n,a] <- c(T[kNNs[neig],a],T[i,a])[1+round(runif(1),0)]
      
    }
  }
  newCases <- data.frame(new)
  for(a in nomatr)
    newCases[,a] <- factor(newCases[,a],levels=1:nlevels(data[,a]),labels=levels(data[,a]))
  
  newCases[,tgt] <- factor(rep(data[1,tgt],nrow(newCases)),levels=levels(data[,tgt]))
  colnames(newCases) <- colnames(data)
  newCases
}


lofactor <- function(data,k) {
  
  data <- as.matrix(data)
  
  # obtain the k nearest neighbors and their distance from each observation
  distdata <- dist.to.knn(data,k)
  p <- dim(distdata)[2L]
  
  # calculate the local reachability density for each observation in data
  lrddata <- reachability(distdata,k)
  
  lof <- rep(0,p)
  
  # computer the local outlier factor of each observation in data
  for (i in 1:p) {
    nneigh <- distdata[2,i]-distdata[1,i]+1
    j <- seq(0,(nneigh-1))
    local.factor <- sum(lrddata[distdata[3+j,i]]/lrddata[i])/nneigh
    lof[i] <- local.factor
  }
  
  # return lof, a vector with the local outlier factor of each observation
  lof
}




# =====================================================
# Function that returns an object in which each column
# contains the indices of the first k neighours followed
# by the distances to each of these neighbours
# =====================================================
dist.to.knn <- function(dataset,neighbors) {
  
  numrow <- dim(dataset)[1L]
  mxNN <- neighbors*2+2
  knndist <- matrix(0,nrow=mxNN,ncol=numrow)
  
  
  for (i in 1:numrow) {
    # find obervations that make up the k-distance neighborhood for observation dataset[i,]
    neighdist <- knneigh.vect(dataset[i,],dataset,neighbors)
    
    x <- length(neighdist)
    if (x > mxNN) {
      knndist <- rbind(knndist,matrix(rep(0,(x-mxNN)*numrow),ncol=numrow))
      mxNN <- x
    }
    knndist[1:x,i] <- neighdist
  }
  
  return(knndist[1:mxNN,])
}



# =====================================================
# Function that returns the distance from a vector "x" to   
# its k-nearest-neighbors in the matrix "data"
# =====================================================
knneigh.vect <- function(x,data,k) {
  
  temp <- as.matrix(data)
  numrow <- dim(data)[1L]
  dimnames(temp) <- NULL
  
  # subtract rowvector x from each row of data
  difference<- scale(temp, x, FALSE)
  
  # square and add all differences and then take the square root
  dtemp <- drop(difference^2 %*% rep(1, ncol(data)))
  dtemp <- sqrt(dtemp)
  
  # order the distances
  order.dist <- order(dtemp)
  nndist <- dtemp[order.dist]
  
  # find distance to k-nearest neighbor
  # uses k+1 since first distance in vector is a 0
  knndist <- nndist[k+1]
  
  # find neighborhood
  # eliminate first row of zeros from neighborhood 
  neighborhood <- drop(nndist[nndist<=knndist])
  neighborhood <- neighborhood[-1]
  numneigh <- length(neighborhood)
  
  # find indexes of each neighbor in the neighborhood
  index.neigh <- order.dist[1:numneigh+1]
  
  # this will become the index of the distance to first neighbor
  num1 <- numneigh+3
  
  # this will become the index of the distance to last neighbor
  num2 <- numneigh+numneigh+2
  
  
  return(c(num1,num2,index.neigh,neighborhood))
}


# =====================================================
# Function that calculates the local reachability density
# of Breuing(2000) for each observation in a matrix, using
# a matrix (distdata) of k nearest neighbors computed by the function dist.to.knn2
# =====================================================
# Example run:
# data(iris)
# 
#
reachability <- function(distdata,k) {
  
  p <- dim(distdata)[2]
  lrd <- rep(0,p)
  
  for (i in 1:p) {
    j <- seq(3,3+(distdata[2,i]-distdata[1,i]))
    # compare the k-distance from each observation to its kth neighbor
    # to the actual distance between each observation and its neighbors
    numneigh <- distdata[2,i]-distdata[1,i]+1
    temp <- rbind(diag(distdata[distdata[2,distdata[j,i]],distdata[j,i]]),distdata[j+numneigh,i])
    
    #calculate reachability
    reach <- 1/(sum(apply(temp,2,max))/numneigh)
    lrd[i] <- reach
  }
  lrd
}