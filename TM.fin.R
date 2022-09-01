#Transition matrix optimization only
library(Rglpk)
noloci=1000
noinds=2

#linkage group 'beginning' loci
#lgend=c(1,50)
lgend=c(1,300,600,900)
#Number of feasible recombinations
norec=8

###################################################
#Building Constraint Matrices
###################################################

#Simulated marker effect matrix
effects=rnorm(noloci,2,0.01)
geno=array(sample(c(0,1),noloci*2*noinds,replace=TRUE),dim=c(noloci,2,noinds))
efmat=array(rep(effects,noinds),dim=c(noloci,2,noinds))

#Marker effect genotypic matrix
data=geno*effects
data=cbind(data[,,1],data[,,2])

#Transition matrix
#3 Transition matrix
#sumj(Tijk)=1
#TM:noloci*16/4 - 4
#TM2:16*noloci,noloci
#T11k,T21k,T31k,T41k,T22k,T12k,T32k,T42k,T33k,T13k,T23k,T43k,T44k,T14k,T24k,T34k

#Constraints of row1 to row2; inflow must equal outflow

one=array(0,32)
one[1:4]=1
one[16+c(1,6,10,14)]=-1

two=array(0,32)
two[5:8]=1
two[16+c(2,5,11,15)]=-1

three=array(0,32)
three[9:12]=1
three[16+c(3,7,9,16)]=-1

four=array(0,32)
four[13:16]=1
four[16+c(4,8,12,13)]=-1

tmat<-function(one,two,three,four){
#Constraint: Inflow must equal outflow
	output=array(0,dim=c(4,noloci*16))
	output[1,1:32]=one
	output[2,1:32]=two
	output[3,1:32]=three
	output[4,1:32]=four
	write.table(file='output.csv',sep=",",col.names=FALSE,row.names=FALSE,output)

	op2=list()
	for(i in 2:(noloci-1)){
		print(i)
		start.index=16*(i-1)+1
		hold=array(0,dim=c(4,noloci*16))
		hold[1,(start.index):(start.index+32-1)]=one
		hold[2,(start.index):(start.index+32-1)]=two
		hold[3,(start.index):(start.index+32-1)]=three
		hold[4,(start.index):(start.index+32-1)]=four
		#output<-rbind(output,hold)
		write.table(file='output.csv',sep=",",append=TRUE,col.names=FALSE,row.names=FALSE,hold)
	}
	output<-read.csv('output.csv',header=FALSE)
	oo=lgend*16-16+1
	tt=oo+15

#Constraint: that only 1 'decision' can occur at a locus (i.e. recombine vs. not recombine)
		tm2<-array(0,dim=c(ncol(output)/16,noloci*16))
		tmzeros<-array(0,dim=c(ncol(output)/16,noloci*16))
		a2=seq(1,ncol(output),16)
		b2=seq(16,ncol(output),16)
		
		#Don't start on transition state at "1st locus"
		tm2[1,c(1,5,9,13)]=1

		for(i in 2:(noloci)){
			#If no loci upstream on linkage group, f
			#if(i%in%lgend){
				#tm2[i,c(1,5,9,13)+16*(i-1)]=1
				#}else{
			tm2[i,(16*(i-1)+1):(16*i)]=1
		}
	#}

#Targeted recombination constraint: no more than norec recombinations
#Accounting for linkage group ends
	trc=rep(c(0,1,1,1),4*noloci)
			for(i in 1:length(oo)){
				trc[oo[i]:tt[i]]=0
			}

#LHS
lhs=rbind(output,tm2,trc)
sign=c(rep("==",nrow(output)),rep("==",nrow(tm2)),"<=")
rhs=c(rep(0,nrow(output)),rep(1,nrow(tm2)),norec)
	return(list(lhs,sign,rhs))
}

tm=tmat(one,two,three,four)

OF=c()
for(i in 1:nrow(data)){
	for(j in 1:ncol(data)){
			OF=c(OF,rep(data[i,j],4))
		}
	}

const=tm[[1]]
const.dir=tm[[2]]
const.rhs=tm[[3]]
xa=rep("B",ncol(const))
system.time({result=Rglpk_solve_LP(OF, const, const.dir, const.rhs, types = xa, control = list("verbose" =
TRUE), max = TRUE)})

########################################################################################################
#Viewing Optimal Solution
########################################################################################################

x=c()
res=result$solution
for(i in 1:noloci){
x<-rbind(x,res[1:16])
res=res[-c(1:16)]
}
a<-rowSums(x[,1:4])
b<-rowSums(x[,5:8])
c<-rowSums(x[,9:12])
d<-rowSums(x[,13:16])
fin=cbind(a,b,c,d)

opthap=data*fin
sum(data*fin)
#Total number of recombinations 

	oo=lgend*16-16+1
	tt=oo+15

	del=c()
	for(i in length(oo)){
	del=c(del,oo[i]:tt[i])
	}

res=result$solution
res[del]=0
sum(res)

x=c()
for(i in 1:noloci){
x<-rbind(x,res[1:16])
res=res[-c(1:16)]
}
a<-rowSums(x[,1:4])
b<-rowSums(x[,5:8])
c<-rowSums(x[,9:12])
d<-rowSums(x[,13:16])
fin=cbind(a,b,c,d)

#return the targeted rec
cnt=0
recpts=c()
for(i in c(2,3,4,6,7,8,10,11,12,14,15)){
counter=counter+sum(which(x[,i]==1)%in%lgend)
print(which(x[,i]==1)%in%lgend)
recpts=c(recpts,which(x[,i]==1))
}
targrec=recpts[which(!recpts%in%lgend)]

pos<-c()
chr<-c()
for(i in 1:length(targrec)){
qw=max(which(lgend<targrec[i]))
chr=c(chr,qw)
pp=sort((qw-1)*300)
pos=c(pos,targrec[i]-pp)
}

recpts=data.frame(Chr=chr,TC.locus.BP=pos)
###################################################
#Brute force optimization of mate-pair
#Integer programming to determine optimal recombination points
###################################################
#Simulated marker effect matrix
effects=rnorm(noloci,2,0.01)
geno=array(sample(c(0,1),noloci*2*noinds,replace=TRUE),dim=c(noloci,2,noinds))
efmat=array(rep(effects,noinds),dim=c(noloci,2,noinds))

#Marker effect genotypic matrix
datahold=geno*effects

combs<-t(combn(noinds,2))
ofval<-c()
for(i in 1:nrow(combs)){
print(i)
data=cbind(datahold[,,combs[i,1]],datahold[,,combs[i,2]])

OF=c()
for(j in 1:nrow(data)){
	for(k in 1:ncol(data)){
			OF=c(OF,rep(data[j,k],4))
		}
	}

result=Rglpk_solve_LP(OF, const, const.dir, const.rhs, types = xa, max = TRUE)
ofval<-c(ofval,result$optimum)
if(result$optimum>=max(ofval)){
	y<-result
	print(y$optimum)
}
}

