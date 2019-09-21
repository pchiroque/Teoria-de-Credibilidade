#...........................................................................
#   BS Model 
#...........................................................................
#
# Author : Pamela Ch. Solano
#          24/03/2018
#          Universidade Federal do Rio de Janeiro 
#          Departamento de M?todos estat?sticos
#          Rio de Janeiro, Ilha do Governador, Fund?o 
#          Email : pchiroque@gmail.com
#                : pamela@dme.ufrj.br
#
#...........................................................................

rm(list=ls(all=TRUE))
root <- getwd()
setwd(root)


for(u in 1){
#################### Dados ###################

  v=c(729,786,872,951,1019,1631,1802,2090,2300,2368,
      796,827,874,917,944,3152,3454,3715,3859,4198,400,420,422,424,440)
  s=c(583,1100,262,837,1630,99,1298,326,463,895,
      1433,496,699,1742,1038,1765,4145,3121,4129,3358,40,0,169,1018,44)
x <- s/v
    
DataMBS <- list(
  v=c(729,786,872,951,1019,1631,1802,2090,2300,2368,
          796,827,874,917,944,3152,3454,3715,3859,4198,400,420,422,424,440),
  i=rep(1:5,each=5),
  x = x
)

#################### Modelo Poisson-Gamma ###################

    
BS <- function (){
  mu0 ~dnorm(0, 0.001)
  tau ~dunif(0, 10)
  sigma ~dunif(0,100)
  invtau2 <- 1/ (tau*tau)
  sigma2<-sigma*sigma
  tau2<-tau*tau
  invsigma2 <- 1/(sigma*sigma)
  kappa<-sigma2/tau2
  for(k in 1:25){
    invsigmax[k] <-v[k]*invsigma2
  }
  
  for(k in 1:5) {
    theta[k] ~dnorm(mu0, invtau2)
    vp[k]<-sum(v[ ( ((k-1)*5) +1):( 5*k)])
    omega[k] <-vp[k] /(vp[k] + (sigma2/tau2))
  }
  for( j  in  1:25) {
    x[j]  ~ dnorm(theta[i[j]], invsigmax[j])
  }
}

require(R2WinBUGS)
write.model(BS, file.path(getwd(), "BS.txt"))

parMBS <- c("mu0","sigma2","tau2", "theta","omega")


require(R2jags)


n.iter=60000;
n.burnin=40000;
n.thin=2;

#################### M_P-G ###################

ptmMBS <- proc.time()    
MBS <- jags(data = DataMBS, #inits = inits,
           model.file = "BS.txt",
           parameters.to.save = parMBS,
           n.chains = 1,
           n.iter=n.iter,
           n.burnin=n.burnin,
           n.thin=n.thin,
           DIC=TRUE)
proc.time() 

manter= c("parMBS","MBS","ptmMBS") #quero

rm(list=setdiff(ls(), manter)) #remove 

save(list = ls(), file=paste("modelBS_.rda",sep="") )#salvo .rda ? bem leve
}


MBS$BUGSoutput$summary
test <- as.mcmc(MBS)

plot(density(MBS$BUGSoutput$sims.list$mu0,bw = 0.5),main="",xlab=expression(mu[0]))
points(quantile(MBS$BUGSoutput$sims.list$mu0,0.025),0,pch=20,col=2)
points(quantile(MBS$BUGSoutput$sims.list$mu0,0.975),0,pch=20,col=2)
points(mean(MBS$BUGSoutput$sims.list$mu0),0,pch=20,col=4)

plot(density(MBS$BUGSoutput$sims.list$sigma2,bw = 30),main="",xlab=expression(sigma^2))
points(quantile(MBS$BUGSoutput$sims.list$sigma2,0.025),0,pch=20,col=2)
points(quantile(MBS$BUGSoutput$sims.list$sigma2,0.975),0,pch=20,col=2)
points(mean(MBS$BUGSoutput$sims.list$sigma2),0,pch=20,col=4)

plot(density(MBS$BUGSoutput$sims.list$tau2,bw = .6),main="",xlab=expression(tau^2),xlim=c(0,6))
points(quantile(MBS$BUGSoutput$sims.list$tau2,0.025),0,pch=20,col=2)
points(quantile(MBS$BUGSoutput$sims.list$tau2,0.975),0,pch=20,col=2)
points(mean(MBS$BUGSoutput$sims.list$tau2),0,pch=20,col=4)

omega.cad <- expression(omega[1],omega[2],omega[3],omega[4],omega[5])

for(i in 1:5){
  plot(density(MBS$BUGSoutput$sims.list$omega[,i],bw = 0.2),main="",xlab=bquote(~.(omega.cad[[i]])))
  points(quantile(MBS$BUGSoutput$sims.list$omega[,i],0.025),0,pch=20,col=2)
  points(quantile(MBS$BUGSoutput$sims.list$omega[,i],0.975),0,pch=20,col=2)
  points(mean(MBS$BUGSoutput$sims.list$omega[,i]),0,pch=20,col=4)
}

theta.cad <- expression(mu(theta[1]),mu(theta[2]),mu(theta[3]),mu(theta[4]),mu(theta[5]))

for(i in 1:5){
  plot(density(MBS$BUGSoutput$sims.list$theta[,i],bw = 0.2),main="",xlab=bquote(~.(theta.cad[[i]])))
  points(quantile(MBS$BUGSoutput$sims.list$theta[,i],0.025),0,pch=20,col=2)
  points(quantile(MBS$BUGSoutput$sims.list$theta[,i],0.975),0,pch=20,col=2)
  points(mean(MBS$BUGSoutput$sims.list$theta[,i]),0,pch=20,col=4)
}

MBS$BUGSoutput$DIC





