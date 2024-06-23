#X<-matrix(c(1,2,3,4),nrow=2)

normalizeMatrix <- function(X){
  normMat = norm(X,type='2')
  return(X/normMat)
}

### Example
#X %>% normalizeMatrix()
#