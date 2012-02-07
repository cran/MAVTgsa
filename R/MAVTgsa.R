######### Hotelling's T^2 statistic using Shrinkage covariance matrix estimates ###############
######### x:data matrix; Row: sample; Column: variables(genes)                  ###############
######### y: vector defining two-group of the samples                           ###############


library(corpcor)

Hott2 <- function(x, y, var.equal=TRUE){
    if(!is.null(y)){
         cl <- as.factor(y)
         lab <- levels(cl)
         ind1 <- which(cl==lab[1])
         ind2 <- which(cl==lab[2])
     }
     if(is.null(ind1) | is.null(ind2))stop("Error: classes 1 and 2 are undefined.")

    data1 <- x[ind1,]
    data2 <- x[ind2,]
    nn <- dim(data1)[1]
    nd <- dim(data2)[1]
    SN <- cov.shrink(data1,verbose=FALSE, lambda.var=0)
    SD <- cov.shrink(data2,verbose=FALSE, lambda.var=0)
    xbar1 <- colMeans(data1)
    xbar2 <- colMeans(data2)
    xdiff <- xbar1 - xbar2

    if (var.equal){
       S <- ((nd - 1) * SD + (nn - 1) * SN)/(nd + nn - 2)
       t2 <- ((nd * nn)/(nd + nn)) * (xdiff %*% solve(S) %*% xdiff)
    }
    else{
       S <- SN/nn + SD/nd
       t2 <- xdiff %*% solve(S) %*% xdiff
    }

    return(as.numeric(t2))
}

######### Ordinary least square(OLS) statistic                                  ###############
######### x:data matrix; Row: sample; Column: variables(genes)                  ###############
######### y: vector defining two-group of the samples  
Tols <- function (x, y ) {
if(!is.null(y)){
         cl <- as.factor(y)
         lab <- levels(cl)
         ind1 <- which(cl==lab[1])
         ind2 <- which(cl==lab[2])
     }
     if(is.null(ind1) | is.null(ind2))stop("Error: classes 1 and 2 are undefined.")
    data1 <- x[ind1,]
    data2 <- x[ind2,]
    nn <- dim(data1)[1]
    nd <- dim(data2)[1]
    #-------------O'Brien's OLS----------------
#    mk <- colMeans(x)
#    s1 <- cov(data1)
#    s2 <- cov(data2)
#    sk2 <- ((nn-1)*s1+(nd-1)*s2)/(nn+nd-2)  #pooled sample covariance matrix
#    sk <- sqrt(diag(sk2))                   #pooled sd  
#    data1s <- colMeans(t((t(data1)-mk)/sk)) #mean of standardize variable 
#    data2s <- colMeans(t((t(data2)-mk)/sk)) 
#    zdif <- data1s-data2s                    #difference between two groups
#    s1s <- cov(t((t(data1)-mk)/sk))          
#    s2s <- cov(t((t(data2)-mk)/sk))
#    sks <- ((nn-1)*s1s+(nd-1)*s2s)/(nn+nd-2)  # pooled sample covariance matrix of zdif
#    J <- rep(1,length(zdif)) 
#    tols <- sqrt((nd * nn)/(nd + nn))*((t(J)%*%zdif)/(t(J)%*%sks%*%J)^0.5)
#--------------------------------------
  olsf <- function(x){
    t.test(x[ind1],x[ind2],var.equal=TRUE)$statistic
  }
 
  t.emat <- apply(x,2,olsf)
  cor.emat <- cor(x)
  tols <- sum(t.emat)/sqrt(sum(cor.emat))
   

    return(abs(as.numeric(tols)))
}




######### Wilks Lambda for n-group comparisons                           ###############
######### Y:data matrix; Row: sample; Column: variables(genes)               ###############
######### Class: vector defining the clinical outcome of the samples         ###############
######### type: type of contrast;                                            ###############
######### base: which group is considered the baseline group for Dunnett contrasts ###############
library(MASS)
library(multcomp)
"design.matrix" <-
function(factors)
{
### Keep the initial order
  n<-length(factors)
  fac<-factor(factors)
  lev<-levels(fac)
  l<-length(lev)
  if(l<2)
    stop("Should have at least two groups")
  X<-matrix(0,n,l)
  for(i in 1:l)
    X[factors==lev[i],i]<-1
  X
}



library(MASS)
ma.estimate <- function (Y, X) {
    ginv(t(X) %*% X) %*% t(X) %*% Y
}



Wilksn <- function(Y, class,type = c("Tukey", "Dunnett", "Sequence"), base=1){
   X <- design.matrix(class)
   hatB <- ma.estimate(Y, X)
   type <- match.arg(type)
   lab <- unique(class)
   #indnames<-paste("ind",1:length(lab),sep="")
   E=0
   for(i in 1:length(lab)){
   ind0 <- which(class==lab[i])
   data0<- Y[ind0,]
   n0   <- length(ind0)
   S0   <- cov.shrink(data0,verbose=FALSE)
   E    <- E+(n0-1)*S0
   }

   k<- length(lab)
    names(Y) <- paste("T", 1:k, sep="")
    CM <- c()
    rnames <- c()
    if (!is.null(names(Y)))
        varnames <- names(Y)
    else varnames <- 1:length(Y)
    kindx <- 1:k
    switch(type, Dunnett = {
        for (i in kindx[-base]) CM <- rbind(CM, as.numeric(kindx == i) - as.numeric(kindx == base))
        rnames <- paste(varnames[kindx[-base]], "-", varnames[base])
		},
		          Tukey = {
        for (i in 1:(k - 1)) {
            for (j in (i + 1):k) {
                CM <- rbind(CM, as.numeric(kindx == j) - as.numeric(kindx == i))
                rnames <- c(rnames, paste(varnames[j], "-", varnames[i]))
            }}}, 
	              Sequence = {
        for (i in 2:k) {
            CM <- rbind(CM, as.numeric(kindx == i) - as.numeric(kindx == i - 1))
            rnames <- c(rnames, paste(varnames[i], "-", varnames[i - 1]))
        }})


   L <- CM
   H <- t(hatB) %*% t(L) %*% ginv(L %*% ginv(t(X) %*% X) %*% t(L)) %*% L %*% hatB
   temp <- solve(E) %*%H
   stat <- 1/det(temp+diag(1,dim(temp)[1],dim(temp)[1]))


   for (i in 1:nrow(L)){
   D  <- matrix(0,nrow(L),ncol(L))
   D[i,]=1
   L0  <- D*L
   H0 <- t(hatB) %*% t(L0) %*% ginv(L0 %*% ginv(t(X) %*% X) %*% t(L0)) %*% L0 %*% hatB
   temp0 <- solve(E) %*%H0
   stat0 <- 1/det(temp0+diag(1,dim(temp0)[1],dim(temp0)[1]))
   stat=c(stat, stat0)

   }
   names(stat)=c("ANOVA", rnames)
   return(stat)
}


library(MASS)
library(multcomp)
library(corpcor)

MAVTn <- function(DATA, GS, MCP=1 , alpha=0.01, nbPerm=5000){


 # DATA : expression data with rows=genes, columns=samples
 #        Note that the first row is the group of sample
 #
 # GS : gene sets
 #      -> a data matrix with rows=genes,
 #                        columns= gene sets,
 #                        GS[i,j]=1 if gene i in gene set j
 #                        GS[i,j]=0 otherwise
 #
 #
 # alpha: signinificant level
 #
 #
 # MCP: type of multiple comparison methods;
 #      Dunnett = 1, Tuckey = 2, Sequential Group = 3 
 # 
     cl <- as.numeric(DATA[1,])
	 DATA <- DATA[-1,]
	 genes <- rownames(DATA)
	 k <- length(levels(factor(cl)))
     n.Samples  <- ncol(DATA)
	 n.GeneSets <- ncol(GS)
     GeneSets.sizes <- apply(GS,2, function(z) sum(z==1))
     
	 if (k>=3) {
	 base <- 1
	 if (MCP==1){type <- "Dunnett"}
	 if (MCP==2){type <- "Tukey"}
	 if (MCP==3){type <- "Sequence"}
     if ((1*(MCP==1)+1*(MCP==2)+1*(MCP==3))==0) stop("Error: MCP must be 1, 2 or 3")
	 } 

     # observed statitic for each gene set
     if (k==2) {
      stat.ols.obs <- apply(GS, 2, function(z) Tols(t(DATA[which(z==1),]), cl))
      stat.hott.obs <- apply(GS, 2, function(z) Hott2(t(DATA[which(z==1),]), cl))
      }
     if (k>=3) stat.obs <- apply(GS, 2, function(z) Wilksn(t(DATA[which(z==1),]), cl, type , base))

     # stats obtained on 'permuted' data

     stat.ols.permut <- matrix(NA,nbPerm,n.GeneSets)
	 stat.hott.permut <- matrix(NA,nbPerm,n.GeneSets)

     if (k>=3){
     stat.temp.permut <- matrix(NA,nbPerm*nrow(stat.obs),n.GeneSets)
     }
      for(i in 1:nbPerm) {
         ind <- sample(n.Samples)
         p.data <- DATA[,ind]
         if (k==2){
           stat.ols.permut[i,] <- apply(GS, 2, function(z) Tols(t(p.data[which(z==1),]), cl))
           stat.hott.permut[i,] <- apply(GS, 2, function(z) Hott2(t(p.data[which(z==1),]), cl))
           }

         if (k>=3){
         k.p <- apply(GS, 2, function(z) Wilksn(t(p.data[which(z==1),]), cl, type , base))
           stat.temp.permut[c(seq(1,nrow(stat.obs)))+nrow(stat.obs)*(i-1),] <- k.p
           }
         }



     if (k==2) {
        GeneSets.ols.pval <- apply(t(stat.ols.permut) >= stat.ols.obs, 1, sum)/nbPerm
		GeneSets.hott.pval <- apply(t(stat.hott.permut) >= stat.hott.obs, 1, sum)/nbPerm
     s.pvalue <- apply(DATA,1, function(z) unlist(summary(aov(z~cl)))["Pr(>F)1"])
     nb.ols.sign <- which(GeneSets.ols.pval<=alpha)
     sg.ols.pvalue <- apply(as.matrix(GS[,nb.ols.sign]),2,function(z) c(list("GS size"=sum(z==1),"p-value"=round(s.pvalue[z==1],4))))
     nb.hott.sign <- which(GeneSets.hott.pval<=alpha)
     sg.hott.pvalue <- apply(as.matrix(GS[,nb.hott.sign]),2,function(z) c(list("GS size"=sum(z==1),"p-value"=round(s.pvalue[z==1],4))))
     ols.star <- c()
     hott.star <- c()
     for (i in 1:n.GeneSets){
     if (GeneSets.ols.pval[i]<=alpha) ols.star[i] <- "*"
     else ols.star[i] <- " "

     if (GeneSets.hott.pval[i]<=alpha) hott.star[i] <-"*"
     else hott.star[i] <- " "
     }
     ols.sign <- data.frame(ols.star)
     colnames(ols.sign) <- " "
     hott.sign <- data.frame(hott.star)
     colnames(hott.sign) <- " "     
        res <- as.data.frame(cbind("GS size"              = GeneSets.sizes,
		                           "OLS p-value"          = GeneSets.ols.pval,
                                       ols.sign,
                                   "T_square p-value"     = GeneSets.hott.pval,
                                        hott.sign))
     result <- list("p value"=res,"Singinificant gene set (one-sided)"=sg.ols.pvalue,"Singinificant gene set (two-sided)"=sg.hott.pvalue)
     }
     if (k>=3){
         stat.temp.pval <- matrix(NA,nrow(stat.obs),n.GeneSets)
         for (i in 1:nrow(stat.obs)){
         t.pval=apply(t(stat.temp.permut[seq(i,nrow(stat.temp.permut),nrow(stat.obs)),])<= stat.obs[i,],1,sum)/nbPerm
           stat.temp.pval[i,] <- t.pval
           }
     r1<-(1:(k*(k-1)/2))
     r2<-NULL
     r3<-NULL
     for(i in 1:k){
     r2<-c(r2,rep(i,k-i))
     if (i+1<=k)
     r3<-c(r3,((i+1):k))
     }
     names=paste(rownames(k.p)," p-value",sep="")
     dimnames(stat.temp.pval) <- list(names)
     s.pvalue <- apply(DATA,1, function(z) unlist(summary(aov(z~cl)))["Pr(>F)1"])
     nb.sign <- which(stat.temp.pval[1,]<=alpha)
     sg.pvalue <- apply(as.matrix(GS[,nb.sign]),2,function(z) c(list("GS size"=sum(z==1),"p-value"=round(s.pvalue[z==1],4))))
     
     star <- c()
     for (i in 1:n.GeneSets){
     if (stat.temp.pval[1,i]<=alpha) star[i] <- "*"
     else star[i] <- " "
     }
     sign <- data.frame(star)
     colnames(sign) <- " "
     res <- as.data.frame(cbind("GS size"              = GeneSets.sizes,
                                "MANOVA p-value"   =t(stat.temp.pval)[,1],
                                sign,
                                t(stat.temp.pval)[,-1]     ))
   
     rownames(res)=colnames(GS)
     result <- list("p value"=res, "singinificant gene set"=sg.pvalue)
     }
	 #Single gene t-test

	 
	 
	 
   
	return(result)
}


