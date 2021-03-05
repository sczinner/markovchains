library(plot.matrix)
###### generate (discrete finite) markov chains
num_states<-4


chain_length = 100
num_chains = 10


chains<-matrix(data=NA,nrow=chain_length,ncol=num_chains)

#initial state
rrvector<-c(0,runif(num_states-1))
orderedvector<-sort(rrvector)#sort
pi<-c(orderedvector[2:num_states],1)-orderedvector
pi_cdf<-cumsum(pi)
rr<-runif(num_chains)
pi_realized<-matrix(data=pi_cdf,nrow = num_chains,ncol=num_states,byrow=TRUE)<=rr
chains[1,]=apply(pi_realized,1,function(x){(which(x==F))[1]})#the new state is the first false value in the row


randomProbMatrix<-function(n){
  rrmatrix<-matrix(data=c(rep(0,n),runif((n-1)*n)),ncol=n,nrow=n)
  orderedmatrix<-t(apply(rrmatrix,1,sort))# sorts rows
  probmatrix<-cbind(orderedmatrix[,2:n],rep(1,n))-orderedmatrix#rows of this matrix will contain all positive, and sum to 1
  return (probmatrix)
}

PP<-randomProbMatrix(num_states)


for(nn in 2:chain_length){
  trans_probs<-PP[chains[nn-1,,drop=F],,drop=F]
  trans_cdf<-t(apply(trans_probs,1,cumsum))
  rr<-runif(num_chains)
  trans_realized<-trans_cdf<=rr
  chains[nn,]=apply(trans_realized,1,function(x){(which(x==F))[1]})#the new state is the first false value in the row
}

matplot(chains, type='l', lty=1, col=1:num_chains, ylim=c(0,num_states), ylab='state', xlab='time')

###### fitting
pi_hat<-vector(length=num_states)#estimated initial state vector
PP_hat<-matrix(data=NA,nrow=num_states,ncol=num_states)#estimated transition matrix

pi_hat<-apply(matrix(chains[1,], nrow=num_states, ncol=num_chains, byrow=TRUE)==(1:num_states),1,sum)/num_chains

#P(Xn in j| Xn-1 in i)=P(Xn in j, Xn-1 in i)/P(Xn-1 in i) so
Nij<-matrix(data=0,nrow=num_states,ncol=num_states)#number of times Xn is in j where Xn-1 is in i
Ni<-vector(mode="numeric",length=num_states)#number of time Xn-1 is in i

for(nn in 2:chain_length){
  j<-chains[nn,]
  i<-chains[nn-1,]
  #non-vector friendly way
  for(chain in 1:length(i)){
      Nij[i[chain],j[chain]]=Nij[i[chain],j[chain]]+1
      Ni[i[chain]]=Ni[i[chain]]+1
  }
  
}
PP_hat<-Nij/matrix(data=Ni,ncol=num_states,nrow=num_states,byrow=FALSE)


plot(PP_hat)
plot(PP)
plot(as.matrix(pi))
plot(as.matrix(pi_hat))

