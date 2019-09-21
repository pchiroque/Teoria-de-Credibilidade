#...........................................................................
#
# Projeto: Onix 
# link : http://www2.susep.gov.br/menuestatistica/Autoseg/menu1.aspx
# Variáveis usadas: Paseio nacional/GM chevrolet onix/ RJ Met RJ/sexo/idade/2016
# Programas de estimação do prêmio
#
#...........................................................................
#
# Author : Pamela Ch. Solano
#          10/09/2019
#          Universidade Federal do Rio de Janeiro 
#          Departamento de Métodos estatísticos
#          Rio de Janeiro, Ilha do Governador, Fund?o 
#          Email : pchiroque@gmail.com
#                : pamela@dme.ufrj.br
#
#...........................................................................

rm(list = ls())

root <- paste(getwd(),"/",sep = "")
Onix <- read.table(paste(root,"Dados.txt",sep = ""),header = TRUE, sep = "\t", dec = ".")
attach(Onix)
summary(Onix)

plot(cbind(Onix$ISMedia,Onix$PremioMedio)~Onix$SexoCondutor)
plot(Onix[,c(8,9,10,12,14,16)])

boxplot(Onix$ISMedia~Onix$Ano)
hist(Onix$Expostos)
table(Onix$FreqIncencioeRoubo)

colnames(Onix)

#### Onix-Kappa-conhecido ####

# Para 
# Fixo : Paseio nacional / / GM chevrolet onix // Modelo ano (2015) //RJ Met RJ // sexo //idade //
# Por anos : 2011 - 2016

# v = expostos
# s = indeniz Incêncido e Roubo

### Exer: Indeniz Colisao - Sexo e Faixa etárea

k=1300
v=rbind(c(2,53,272,290),
        c(4,168,985,1244),
        c(7,174,977,1300),
        c(5,151,926,1205),
        c(7,190,1154,1642)
)
s=rbind(c(13512,129379,272187,260276),
        c(0,80107,830011,1198200),
        c(0,215812,817763,1291340),
        c(0,0,441551,880669),
        c(0,128515,812771,1303998)
)

Premio.k=function(k,v,s){
  n_g=dim(v)[1];n_a=dim(v)[2];
  X=s/v;v.=s.=matrix(NA,n_g,1);
  for(i in 1:n_g){v.[i]=sum(v[i,]);s.[i]=sum(s[i,])}
  Xi=s./v.;v..=sum(v.);wi=(v.)/(v.+k);w=sum(wi);
  mu.0=sum((wi/w)*Xi)
  Est.hom=wi*Xi+(1-wi)*mu.0;Premioi=((Est.hom)*v.);
  Premio=sum((Est.hom)*v.)
  return(list(wi=wi,Est.hom=Est.hom,mu.0=mu.0,Premioi=Premioi,Premio=Premio))
}
Premio.k(k,v,s)#wi,Est.hom,mu.0,Premioi,Premio


# .....................................................
# .....................................................

# Repare que não se considerou a tabela FIP 
# HW : Fazer comparação do Prémio considerando desvalorização e não.

# .....................................................
# .....................................................


#### Onix-Kappa-Desconhecido ####

v=rbind(c(2,53,272,290),
        c(4,168,985,1244),
        c(7,174,977,1300),
        c(5,151,926,1205),
        c(7,190,1154,1642)
)

s=rbind(c(13512,129379,272187,260276),
        c(0,80107,830011,1198200),
        c(0,215812,817763,1291340),
        c(0,0,441551,880669),
        c(0,128515,812771,1303998)
)

k.est=function(v,s){
  n_g=dim(v)[1];
  n_a=dim(v)[2]
  X=s/v;  v.=s.=matrix(NA,n_g,1);s.est=matrix(NA,n_g,n_a);
  si=matrix(NA,n_g,1)
  for(i in 1:n_g){v.[i]=sum(v[i,]);s.[i]=sum(s[i,])}
  Xi=s./v.;v..=sum(v.);  X.bar=sum((v./v..)*Xi);
  c=((n_g-1)/n_g)*((sum(((v./v..))*(1-(v./v..))))^(-1))
  d=sum((n_g/(n_g-1))*((Xi-X.bar)^2)*(v./v..))
  for(j in 1:n_a){
    s.est[,j]=(v[,j]*((X[,j]-Xi)^(2)))
  }
  for(i in 1:n_g){
    si[i]=sum(s.est[i,])/(n_a-1)
  }
  sigma.2=mean(si); tau.2=c*(d-(n_g*(sigma.2/v..)))
  k=sigma.2/tau.2
  return(rbind(k,sigma.2,tau.2)  )
}

k.est(v,s)
Premio.k(k.est(v,s)[1],v,s)

#...........................................................................
#   BS Model: Full Bayesian - Onix 2015  ####
#...........................................................................

rm(list=ls(all=TRUE))

root <- paste(getwd(),"/",sep = "")
setwd(root)


for(u in 1){
  #################### Dados ###################
  
  v=c(2,53,272,290,
      4,168,985,1244,
      7,174,977,1300,
      5,151,926,1205,
      7,190,1154,1642)
  
  s= c(13512,129379,272187,260276,
       0,80107,830011,1198200,
       0,215812,817763,1291340,
       0,0,441551,880669,
       0,128515,812771,1303998)
  
  x <- s/v
  
  DataMBS <- list(
    v=c(2,53,272,290,
        4,168,985,1244,
        7,174,977,1300,
        5,151,926,1205,
        7,190,1154,1642),
    i=rep(1:5,each=4),
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
    for(k in 1:20){
      invsigmax[k] <-v[k]*invsigma2
    }
    
    #Estima??o dos pr?mios theta's e do fatores de credibilidades alpha's#
    for(k in 1:5) {
      theta[k] ~dnorm(mu0, invtau2)
      vp[k]<-sum(v[ ( ((k-1)*4) +1):( 4*k)])
      omega[k] <-vp[k] /(vp[k] + (sigma2/tau2))
    }
    for( j  in  1:20) {
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

  root <- paste(getwd(),"/",sep = "")
  
  setwd(root)
  
  save(list = ls(), file=paste(root,"modelOnixBS.rda",sep="") )#salvo .rda ? bem leve
}


MBS$BUGSoutput$summary
test <- as.mcmc(MBS)
plot(test)

plot(density(MBS$BUGSoutput$sims.list$mu0,bw = 0.5),main="",xlab=expression(mu[0]))
points(quantile(MBS$BUGSoutput$sims.list$mu0,0.025),0,pch=20,col=2)
points(quantile(MBS$BUGSoutput$sims.list$mu0,0.975),0,pch=20,col=2)
points(mean(MBS$BUGSoutput$sims.list$mu0),0,pch=20,col=4)

plot(density(MBS$BUGSoutput$sims.list$sigma2,bw = 30),main="",xlab=expression(sigma^2))
points(quantile(MBS$BUGSoutput$sims.list$sigma2,0.025),0,pch=20,col=2)
points(quantile(MBS$BUGSoutput$sims.list$sigma2,0.975),0,pch=20,col=2)
points(mean(MBS$BUGSoutput$sims.list$sigma2),0,pch=20,col=4)

plot(density(MBS$BUGSoutput$sims.list$tau2,bw = .6),main="",xlab=expression(tau^2))#,xlim=c(0,6))
points(quantile(MBS$BUGSoutput$sims.list$tau2,0.025),0,pch=20,col=2)
points(quantile(MBS$BUGSoutput$sims.list$tau2,0.975),0,pch=20,col=2)
points(mean(MBS$BUGSoutput$sims.list$tau2),0,pch=20,col=4)

omega.cad <- expression(omega[1],omega[2],omega[3],omega[4],omega[5])

for(i in 1:5){
  plot(density(MBS$BUGSoutput$sims.list$omega[,i],bw = 0.0002),main="",xlab=bquote(~.(omega.cad[[i]])))
  points(quantile(MBS$BUGSoutput$sims.list$omega[,i],0.025),0,pch=20,col=2)
  points(quantile(MBS$BUGSoutput$sims.list$omega[,i],0.975),0,pch=20,col=2)
  points(mean(MBS$BUGSoutput$sims.list$omega[,i]),0,pch=20,col=4)
}

theta.cad <- expression(mu(theta[1]),mu(theta[2]),mu(theta[3]),mu(theta[4]),mu(theta[5]))

for(i in 1:5){
  plot(density(MBS$BUGSoutput$sims.list$theta[,i],bw = 0.6),main="",xlab=bquote(~.(theta.cad[[i]])))
  points(quantile(MBS$BUGSoutput$sims.list$theta[,i],0.025),0,pch=20,col=2)
  points(quantile(MBS$BUGSoutput$sims.list$theta[,i],0.975),0,pch=20,col=2)
  points(mean(MBS$BUGSoutput$sims.list$theta[,i]),0,pch=20,col=4)
}

MBS$BUGSoutput$DIC
