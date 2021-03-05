PP<-t(matrix(data = c(0,1,0,
                      0,0,1,
                      0.5,0,0.5),ncol = 3))

chain_length = 10
num_chains = 10


chains<-matrix(data=NA,nrow=chain_length,ncol=num_chains)

chains[1,]=1#starting of all of them

for(nn in 2:chain_length){
  trans_probs<-PP[chains[nn-1,,drop=F],,drop=F]
  trans_cdf<-t(apply(trans_probs,1,cumsum))
  rr<-runif(num_chains)
  trans_realized<-trans_cdf<=rr
  chains[nn,]=apply(trans_realized,1,function(x){(which(x==F))[1]})#the new state is the last false value in the row
}

matplot(chains, type='l', lty=1, col=1:num_chains, ylim=c(0,4), ylab='state', xlab='time')