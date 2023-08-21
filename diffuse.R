dat1<-read.csv("/Users/apull/Downloads/cleaned_table.csv")
dat1<-dat1[1:8257,]

library(caret)
library(KRLS)
library(RSpectra)
library(Matrix)
install.packages("scatterplot3d")
library("scatterplot3d")
process <- preProcess(as.data.frame(dat1), method=c("range"))
norm_scale <- predict(process, as.data.frame(dat1)) #normalizing and scaling data so all points are between 0 and 1

eps=.5 #epsilon parameter relates to time scale
K<-gausskernel(norm_scale,eps) #kernal function
R<-colSums(K)
Rp<-diag(1/R)
Rp<-as(Rp,"dgCMatrix")
P<-Rp%*%K
P<-as(P,"dgCMatrix")
D_right<-diag((R)^0.5)
D_right<-as(D_right,"dgCMatrix")
D_left<-diag((1/R^0.5))
D_left<-as(D_left,"dgCMatrix")
P_T<-D_right%*%(P%*%D_left)
### Normalized transition matrix reprenting the probability of moving between two data points in one time step

eig=eigs_sym(P_T,k=22, which='LA')
eig$values
##eigenvectors and values of transition matrix
matplot(seq(1:length(eig$values)),eig$values, type='l',xlab='n',ylab='eigenvalue of P')
##plot of 1st 22 eigenvalues
I<-diag(1,nrow=nrow(P_T), ncol=(ncol(P_T)))
eigL=eigs_sym((4/(sqrt(2)*(eps)))*(P_T-I),k=22, which='LA')
eigL$values
#eigenvectors and eigenvalues of L matrix approximating differential operator
matplot(seq(1:length(eigL$values)),eigL$values, type='l',xlab='n', ylab='eignvalue of L')
#plot of eigenvalues of L


###Check all combinations of 21 variables for optimal lower dimensional representation
###Should look into search algorithms for finding least spectral loss
eig_tot=matrix(data=NA, nrow=22,ncol=22)
eigL_tot=matrix(data=NA, nrow=22,ncol=22)
for(i in 1:22){
  norm_shuffle_21=norm_scale[,-i]
  K<-gausskernel(norm_shuffle_21,eps)
  
  R<-colSums(K)
  Rp<-diag(1/R)
  Rp<-as(Rp,"dgCMatrix")
  P<-Rp%*%K
  P<-as(P,"dgCMatrix")
  D_right<-diag((R)^0.5)
  D_right<-as(D_right,"dgCMatrix")
  D_left<-diag((1/R^0.5))
  D_left<-as(D_left,"dgCMatrix")
  P_T<-D_right%*%(P%*%D_left)
  
  
  eig_21=eigs_sym(P_T,k=22, which='LA')
  eig_tot[i,]=eig_21$values
  plot(seq(1:length(eig_21$values)),eig$values, type='l')
  I<-diag(rowSums(P_T))
  eigL_21=eigs_sym((4/(eps))*(P_T-I),k=22, which='LA')
  eigL_tot[i,]=eigL_21$values
  plot(seq(1:length(eigL_21$values)),eigL$values, type='l')
}
matplot(as.data.frame(t(eig_tot)),type="l",xlab=('n'), ylab=("Eignevalue"))
###plot of 1st 21 eigenvalues for all combinations of 21 variables for P matrix

matplot(as.data.frame(t(eigL_tot)),type="l")
###plot of 1st 21 eigenvalues for all combinations of 21 variables for L matrix

spect_loss=numeric(21)
for(i in 1:21){
  spect_loss[i]=sqrt(sum(eig_tot[i,]-eig$values)^2)
}
spect_loss
min(spect_loss)

### Minimum spectral loss corresponds to variable that can be removed while losing the least amount of information

eig_graph<-rbind(eig_tot[which.min(spect_loss),],eig$values)

###returns eigenvalues that minimize spectral loss and combines with eigenvalues of full data into one matrix

matplot(as.data.frame(t(eig_tot)), type='l', title('Eigenvalues for 21 Dim Representaitons of P')) #graph of largest 22 eigenvalues of P
matplot(as.data.frame(t(eig_graph)), type=c('p','l'),lwd=1,col=c("blue","red"), xlab='n', ylab='eigenvalue')

### graph of eigenvalues for full data (red) and 21 dim data with 19th eigenvalue removed
###spectral loss determined that removing 19th eigenvalue has least effect on data

###Repeat to find 20 dim representation
eig_tot=matrix(data=NA, nrow=22,ncol=22)
eigL_tot=matrix(data=NA, nrow=22,ncol=22)
for(i in 1:21){
  norm_shuffle_20=norm_shuffle_21[,-i]
  K<-gausskernel(norm_shuffle_20,eps)
  
  R<-colSums(K)
  Rp<-diag(1/R)
  Rp<-as(Rp,"dgCMatrix")
  P<-Rp%*%K
  P<-as(P,"dgCMatrix")
  D_right<-diag((R)^0.5)
  D_right<-as(D_right,"dgCMatrix")
  D_left<-diag((1/R^0.5))
  D_left<-as(D_left,"dgCMatrix")
  P_T<-D_right%*%(P%*%D_left)
  
  
  eig_20=eigs_sym(P_T,k=22, which='LA')
  eig_tot[i,]=eig_20$values
  #plot(seq(1:length(eig_20$values)),eig$values, type='l')
  I<-diag(rowSums(P_T))
  eigL_20=eigs_sym((4/(eps))*(P_T-I),k=22, which='LA')
  eigL_tot[i,]=eigL_20$values
  #plot(seq(1:length(eigL_20$values)),eigL$values, type='l')
}
matplot(as.data.frame(t(eig_tot)),type="l",xlab=('n'), ylab=("Eignevalue"))
###plot of 1st 20 eigenvalues for all combinations of 20 variables for P matrix

matplot(as.data.frame(t(eigL_tot)),type="l")
###plot of 1st 20 eigenvalues for all combinations of 20 variables for L matrix

spect_loss=numeric(20)
for(i in 1:20){
  spect_loss[i]=sqrt(sum(eig_tot[i,]-eig$values)^2)
}
spect_loss
min(spect_loss)

eig_graph<-rbind(eig_tot[which.min(spect_loss),],eig$values)
#removes 18th coordinate
###returns eigenvalues that minimize spectral loss and combines with eigenvalues of full data into one matrix

matplot(as.data.frame(t(eig_tot)), type='l') #graph of largest 22 eigenvalues of P
matplot(as.data.frame(t(eig_graph)), type=c('p','l'),lwd=1,col=c("blue","red"), xlab='n', ylab='eigenvalue')

###Repeat to find 19 dim representation
eig_tot=matrix(data=NA, nrow=22,ncol=22)
eigL_tot=matrix(data=NA, nrow=22,ncol=22)
for(i in 1:20){
  norm_shuffle_19=norm_shuffle_20[,-i]
  K<-gausskernel(norm_shuffle_19,eps)
  
  R<-colSums(K)
  Rp<-diag(1/R)
  Rp<-as(Rp,"dgCMatrix")
  P<-Rp%*%K
  P<-as(P,"dgCMatrix")
  D_right<-diag((R)^0.5)
  D_right<-as(D_right,"dgCMatrix")
  D_left<-diag((1/R^0.5))
  D_left<-as(D_left,"dgCMatrix")
  P_T<-D_right%*%(P%*%D_left)
  
  
  eig_19=eigs_sym(P_T,k=22, which='LA')
  eig_tot[i,]=eig_19$values
  #plot(seq(1:length(eig_20$values)),eig$values, type='l')
  I<-diag(rowSums(P_T))
  eigL_19=eigs_sym((4/(eps))*(P_T-I),k=22, which='LA')
  eigL_tot[i,]=eigL_19$values
  #plot(seq(1:length(eigL_19$values)),eigL$values, type='l')
}
matplot(as.data.frame(t(eig_tot)),type="l",xlab=('n'), ylab=("Eignevalue"))
###plot of 1st 20 eigenvalues for all combinations of 20 variables for P matrix

matplot(as.data.frame(t(eigL_tot)),type="l")
###plot of 1st 20 eigenvalues for all combinations of 20 variables for L matrix

spect_loss=numeric(19)
for(i in 1:19){
  spect_loss[i]=sqrt(sum(eig_tot[i,]-eig_20$values)^2)
}
spect_loss
min(spect_loss)

eig_graph<-rbind(eig_tot[which.min(spect_loss),],eig$values)

###returns eigenvalues that minimize spectral loss and combines with eigenvalues of full data into one matrix

matplot(as.data.frame(t(eig_tot)), type='l', title('Eigenvalues for 21 Dim Representaitons of P')) #graph of largest 22 eigenvalues of P
matplot(as.data.frame(t(eig_graph)), type=c('p','l'),lwd=1,col=c("blue","red"), xlab='n', ylab='eigenvalue')
###Exploration of different eigenvector spaces

###Produce eigenvector plot of vectors 1 and 4 and color codes points based on graduation rate
grad_analysis<-cbind(eigL$vectors[,1],eigL$vectors[,4],norm_scale[,11])
for (i in 1:length(grad_analysis[,3])){
  if (norm_scale[i,11]>.90){
    grad_analysis[i,3]='green'
  }
  if (norm_scale[i,11]<=.90){
    grad_analysis[i,3]='yellow'
  }  
  if (norm_scale[i,11]<.8){
    grad_analysis[i,3]='red'
  }

  
}
plot(eigL$vectors[,4],eigL$vectors[,1], col=grad_analysis[,3],title("Graduation Rates eps=.1"))


grad_analysis<-cbind(eigL$vectors[,1],eigL$vectors[,4],norm_scale[,8])
for (i in 1:length(grad_analysis[,3])){
  
  if (norm_scale[i,8]>.19){
    grad_analysis[i,3]='red'
  }
  if (norm_scale[i,8]<=.19){
    grad_analysis[i,3]='orange'
  }  
  if (norm_scale[i,8]<.17){
    grad_analysis[i,3]='yellow'
  }
  if (norm_scale[i,8]<.15){
    grad_analysis[i,3]='green'
  }
  
  
}
plot(eigL$vectors[,1],eigL$vectors[,2], col=grad_analysis[,3],title("Graduation Rates eps=.1"))
scatterplot3d(eigL$vectors[,4],eigL$vectors[,2],eigL$vectors[,3], color=grad_analysis[,3])
