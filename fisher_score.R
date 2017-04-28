#############################################
############## Fisher Score ################# 
#############################################
# this function calculates the Fisher Score of the features of a dataset and ranks them

###Source:  "Generalized Fisher Score for Feature Selection"; Gu, Lo and Han

#input:
# class.labels: vector of class labels
# data: data frame of features, with features in columns

#output: sorted list of Fisher scores along with their indices
#features with large Fisher scores are considered important
fisher.rank<-function(class.labels, data)
{
  #applies fisher.score() to the columns of data
  scores<-apply(data, MARGIN = 2, FUN = fisher.score, class.labels = class.labels)
  scores2<-sort(scores, decreasing = TRUE, index.return = TRUE)
  return(scores2)
}


#calculates Fisher score for given feature
#this function is called on by fisher.rank() and applied to the columns of a data frame
#inputs the class labels and the data vector/feature of the same length
fisher.score<-function(class.labels, feature)
{
  classes<-sort(unique(class.labels))
  class.sizes<-table(class.labels)
  grand.mean<-mean(feature)
  
  denom<-vector()
  num<-vector()
  
  for (i in 1:length(classes))
  {
    num[i]<-class.sizes[i]*(mean(feature[which(class.labels == classes[i])]) - grand.mean)^2
    denom[i]<-class.sizes[i]*var(feature[which(class.labels == classes[i])])
  }
  return(sum(num)/sum(denom))
}
