library(lpSolve)
noloci=10
noinds=100
#SNP genotype matrix
geno=matrix(sample(c(0,1),noloci*2*noinds,replace=TRUE),dim=c(1000,2,1000))

x=seq(-0.3,0.3,0.001)

y=rnorm(x,0,0.1)

#
effects=rnorm(1000,0,0.01)

#Marker effect genotypic matrix
mem=geno*effects

#1 Data a
#Selecting mate-pair constraint
#Output is dataijk 
datacons<-function(noloci,noinds){
	one<-seq(1,noinds*noloci*2,noloci*2)
	two<-seq(noloci*2,noinds*noloci*2,noloci*2)
	three<-cbind(one,two)
	cons=array(0,c(noinds,noloci*2*noinds))
		for(i in 1:nrow(cons)){
			cons[i,three[i,1]:three[i,2]]=1
		}
		sels=-diag(noinds)*(noloci*2)
		ret=list(cons,sels)
		return(ret)
	}

#2 X matrix
xmat=function(noloci,noinds){
		#Data variable matrix
	it=array(0,dim=c(noloci*2,noloci*2*noinds))
		for(i in 1:(noloci*2)){
			it[i,seq(i,noloci*2*noinds,20)]=1
		}
		#X Matrix variable
			x=diag(20)
			itx=list(x,it)
			return(itx)
	}

#3 Transition matrix
tmat=function(noloci){
	iter=(nrow(tm)/4)-4
	tm<-array(0,dim=c(noloci*16,noloci*16))
	ind1=1:4
	ind2=17:20
	for(i in 1:iter){
		a=ind1+4*(i-1)
		b=ind2+4*(i-1)
		tm[i,a]=1
		tm[i,b]=-1
	}
	return(tm)
}


