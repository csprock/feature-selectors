###########################################
### Efficient Feature Selection for LDA ###
###########################################

# source: "Efficient Feature Selection for Linear Discriminant Analysis and its Application to Face Recognition" (Lei, Liao and Li)


#EFS for LDA main function
# sequentially selects features which yield the largest Fisher separation improvements (see source paper for details)
# input: the output of the required.data() function, the correlation weight parameter L
# output: an ordered list of feature indices in the order they were selected
EFS.LDA<-function(req.data, L = 1)
{
  Jmat<-req.data[[1]]
  f.scores<-req.data[[2]]
  cor.matrix<-req.data[[3]]

  m<-length(f.scores)    #size of feature set
  
  #initialize indices of starting feature set
  initial.set<-which.max(f.scores)          #initial feature is one with largest Fisher score
  candidate.set<-setdiff(1:m, initial.set)  #initial candidate set are all features except the initial feature
  
  while (length(initial.set) < m)
  {
    temp<-unlist(lapply(X = candidate.set, FUN = minmax, initial = initial.set, L = L, Jmat = Jmat, f.scores = f.scores, cor.matrix = cor.matrix))
    toAdd<-candidate.set[which.max(temp)]
    initial.set<-append(initial.set, toAdd)
    candidate.set<-setdiff(candidate.set, toAdd)
  }
  names(initial.set)<-NULL
  return(initial.set)
}

#this function computes the required data used by the main feature selection function
#returns matrix of Fisher separation values, the Fisher scores (which are used to select an initial feature set
# and for computing the Fisher separation improvement scores), and the correlation matrix of the data set
required.data<-function(class.labels, data)
{
  
  #compute a matrix where the (i,j) entry is the Fisher separation value between the ith and jth features
  J_values<-matrix(0, dim(data)[2], dim(data)[2])
  for (i in 1:(dim(data)[2] -1))
  {
    for (j in (i+1):dim(data)[2])
    {
      J_values[i,j]<-J(class.labels , data[,i], data[,j])
    }
  }
  J_values<-J_values+t(J_values)
  f.scores<-apply(X = data, MARGIN = 2, FUN = fisher.score, class.labels = class.labels)
  cor.matrix<-cor(data)
  return(list(Jmat,f.scores,cor.matrix))
}


# computes the value of the min-max function
# inputs are the index of a candidate feature, the initial indices, the weight given to 
# the correlation coefficient, the matrix of Fisher separation values, individual Fisher scores, 
# and the correlation matrix of the data set
# returns the min-max
minmax<-function(candidate.idx, initial, L = 1, Jmat, f.scores, cor.matrix)
{
  Jscore<-vector()
  cor.temp<-vector()
  for (i in 1:length(initial))
  {
    Jscore[i]<-Jmat[candidate.idx, initial[i]] - f.scores[initial[i]]
    cor.temp[i]<-cor.matrix[candidate.idx, initial[i]]
  }
  return(min(Jscore) - L*max(cor.temp))
}

#computes the Fisher separation value for two features
J<-function(class.labels, feature, feature_prime)
{
  matrices<-scatter.matrices(class.labels, cbind(feature,feature_prime))
  Sb<-matrices[[1]]
  Sw<-matrices[[2]]
  A<-solve(Sw)%*%Sb
  E<-eigen(A)
  W<-E$vectors
  J<-det(t(W)%*%Sb%*%W)/det(t(W)%*%Sw%*%W)
  return(J)
}

#creates between-class and within-class scatter matrices Sb and Sw
scatter.matrices<-function(class.labels, data)
{
  classes<-unique(class.labels)
  Nc<-table(class.labels)
  N<-dim(data)[1]
  grand.mean<-apply(X = data, MARGIN = 2, FUN = mean)
  class.means<-matrix(0,nrow = length(classes), ncol = dim(data)[2])
  
  #class means in rows
  for (i in 1:length(classes))
  {
    class.means[i,]<-apply(X = data[which(class.labels == classes[i]), ], MARGIN = 2, FUN = mean)
  }
  
  #create Sb matrix
  Sb<-matrix(0,nrow = dim(data)[2], ncol = dim(data)[2])
  Sw<-matrix(0,nrow = dim(data)[2], ncol = dim(data)[2])
  for (i in 1:length(classes))
  {
    Sb<-Sb + Nc[i]*(class.means[i,] - grand.mean)%*%t(class.means[i,] - grand.mean)
    
    Sw.temp<-matrix(0,nrow = dim(data)[2], ncol = dim(data)[2])
    data.temp<-as.matrix(data[which(class.labels == classes[i]), ])
    for (j in 1:Nc[i])
    {
      Sw.temp<-Sw.temp + (data.temp[j,] - class.means[i,])%*%t(data.temp[j,] - class.means[i,])
    }
    Sw<-Sw + Sw.temp
    
  }
  Sb<-(1/N)*Sb
  Sw<-(1/N)*Sw
  return(list(Sb,Sw))
}
