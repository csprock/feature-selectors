#############################################
############## Laplacian Score ############## 
#############################################

### Source:  "Laplacian Score for Feature Selection"; He, Cai and Niyogi

#input:
# data: is a data frame or matrix with features in columns
# k: k for k-nearest neighbors
# o: shape parameter for Gaussian function

#output: sorted list of Lapacian scores and their indices
#features with small Lapacian scores are considered important

laplacian.rank<-function(data, k = 5, o = 1)
{
  #returns index of k nearest neighbors in data.row
  knn.index<-function(k, data.row)
  {
    temp<-sort(data.row, decreasing = FALSE, index.return = TRUE)
    return(temp$ix[2:(k+1)])
  }
  
  #gaussian kernel function
  gaussian.kernel<-function(o, distance)
  {
    return(exp((-distance^2)/o))
  }
  
  #computes laplacian score for input feature
  #inputs: feature vector, weight matrix S, diagonal matrix D
  #output: feature score of vector
  laplacian.score<-function(feature, S,D)
  {
    #create modified feature
    ones<-matrix(1,dim(S)[1],1)
    f.tilde<-feature - (as.numeric(t(feature)%*%D%*%ones)/as.numeric(t(ones)%*%D%*%ones))*ones
    
    #create lapacian matrix
    L<-D-S
    
    #compute numerator and denominator of Lapacian score for the feature
    num<-t(f.tilde)%*%L%*%f.tilde
    denom<-t(f.tilde)%*%D%*%f.tilde
    return(num/denom)
  }
  
  #compute pairwise distances between observations (in rows of data)
  distances<-as.matrix(dist(data, upper=TRUE, diag=TRUE))
  
  #initialize adjacency matrix
  adj_matrix<-matrix(0, dim(distances), dim(distances))
  #fill nearest neighbor matrix
  for (i in 1:dim(adj_matrix)[1])
  {
    idx<-knn.index(k,distances[,i])
    adj_matrix[idx,idx]<-1
  }

  #create weight matrix S
  S<-apply(X = distances, MARGIN = c(1,2), FUN = gaussian.kernel, o = o)
  S<-adj_matrix*S

  #create diagonal matrix D
  D<-matrix(0,dim(S)[1], dim(S)[1])
  diag(D)<-S%*%matrix(1,dim(S)[1],1)

  #apply l.score() to columns of data
  scores<-apply(X = data, MARGIN = 2, FUN = laplacian.score, S = S, D = D)
  #sort Lapacian scores
  scores<-sort(scores, decreasing = FALSE, index.return = TRUE)
  return(scores)
}









