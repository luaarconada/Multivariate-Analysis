#     
#   Squared distance based on Gower's general similarity coefficient for two
#   (row) vectors x,y containing a mixture of continuous, categorical and 
#   binary variables.
#
#   p is a 3-vector, containing (in this order):
#       p[1] = Number of continuous variables
#       p[2] = Number of categorical variables
#       p[3] = Number of binary variables
#       sum(p) is thus the total number of variables. It should be equal to
#       both length(x)and length(y).
#
#   It is assumed that entries in x and y are coded as numerical 
#   (floating point).
#
#   Xr is a p[1]-vector, containing the ranges (max-min) of the continuous
#   variables. It has to be previously computed with 'Xrange'.
#
#   The returned object, d, is the SQUARED Gower distance between x and y.
#

gower<-function(x,y,p,Xr){
   nmatch<-0
   nposmatch<-0
   sc<-0
   nnegmatch<-0
   pc<-p[1]
   pq<-p[2]
   pb<-p[3]
# ------------------------------------------------------------------
#  Quantitative variables
# ------------------------------------------------------------------
   if (pc>0){
      xc<-x[1:pc]
      yc<-y[1:pc]
      sc<-sum(1-abs(xc-yc)/Xr)
      }
# ------------------------------------------------------------------
#  Categorical variables
# ------------------------------------------------------------------
   if (pq>0){
      xq<-x[(pc+1):(pc+pq)]
      yq<-y[(pc+1):(pc+pq)]
      nmatch<-sum(xq==yq)
      }
# ------------------------------------------------------------------
#  Binary variables
# ------------------------------------------------------------------
   if (pb>0){
       xb<-x[(pc+pq+1):sum(p)]
       yb<-y[(pc+pq+1):sum(p)]
       nposmatch<-sum(xb*yb)
       nnegmatch<-(sum((xb-yb)==0)-sum(xb*yb))
       }
   else{
      nposmatch<-0
      nnegmatch<-0
      }
# ------------------------------------------------------------------
   if (nnegmatch==sum(p))
        s<-0
   else
        s<-(sc+nmatch+nposmatch)/(sum(p)-nnegmatch)
# ------------------------------------------------------------------
   d<-(1-s)
 }