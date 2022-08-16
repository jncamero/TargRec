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
#dataijk ≤ noinds*noloci*Sk
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

dc<-datacons(noloci,noinds)

#2 X matrix
#Each value in X is the sum along 3rd dimension of data matrix
xmat=function(noloci,noinds){
		#Data variable matrix
		#LHS
	it=array(0,dim=c(noloci*2,noloci*2*noinds))
		for(i in 1:(noloci*2)){
			it[i,seq(i,noloci*2*noinds,20)]=1
		}
		#RHS
			x=diag(20)
			itx=list(it,x)
			return(itx)
	}

xm<-xmat(noloci,noinds)

#3 Transition matrix
#sumj(Tijk)=1
#TM:noloci*16/4 - 4
#TM2:16*noloci,noloci
#T11k,T21k,T31k,T41k,T22k,T12k,T32k,T42k,T33k,T13k,T23k,T43k,T44k,T14k,T24k,T34k

#Need to modify new linkage group constraints in tm2
#1) T111 + T221 + T331 + T441 = 1
tmat=function(noloci,lgend){
	iter=(noloci*16/4)-4
	tm<-array(0,dim=c(iter,noloci*16))
	ind1=1:4
	ind2=17:20
	for(i in 1:iter){
		a=ind1+4*(i-1)
		b=ind2+4*(i-1)
		tm[i,a]=1
		tm[i,b]=-1
	}

	tm2<-array(0,dim=c(ncol(tm)/16,noloci*16))
		a=seq(1,ncol(tm),16)
		b=seq(16,ncol(tm),16)
		
		for(i in 1:(noloci)){
			if(i%in%lgend){
				tm2[i,c(1,5,9,13)]=1
				}else{
			tm2[i,a[i]:b[i]]=1
		}
	}
	return(list(tm,tm2))
}

tm=tmat(noloci,lgend)
#Markers at end of linkage group
lgend<-1



