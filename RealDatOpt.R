load("/Users/gnyku/Desktop/JohnCameron/NewMarkerData/C1.11.RData")
marker.effects=read.csv("/Users/gnyku/Desktop/JohnCameron/Randomized Trials Effects_GCA/C1.11.csv")
data1=parents_geno.new[,1]
data2=parents_geno.new[,2]

#Number of feasible recombinations
norec=7

#This function takes as input genotypic data for inbred maize line, where
#1 is a homozygote for allele A, -1 homozygote for allele B, and 0 is the heterozygote
geno<-function(data){
parent=c()
marker=c()
for(i in 1:length(data)){
	if(data[i]==1){
		parent=rbind(parent,c(1,1))
		marker=c(marker,row.names(parents_geno.new[i,]))
			}
	if(data[i]==-1){
		parent=rbind(parent,c(-1,-1))
		marker=c(marker,row.names(parents_geno.new[i,]))	
			}	
	if(data[i]==0){
		parent=rbind(parent,sample(c(-1,1),2,replace=FALSE))
		marker=c(marker,row.names(parents_geno.new[i,]))	
		}
	}	
return(data.frame(marker=marker,parent=parent))
}


a=geno(data1)
b=geno(data2)

marks=intersect(a$marker,b$marker)
parA=a[which(a$marker%in%marks),]
parB=b[which(b$marker%in%marks),]
pars=cbind(parA,parB[,2:3])

marker.effects=marker.effects[which(marker.effects[,1]%in%pars$marker),]
marker.effects=marker.effects[match(pars[,1],marker.effects[,1]),]

#####################################################################################
#####################################################################################
#sudo apt-get install glpk-utils libglpk-dev glpk-doc
#Transition matrix optimization only
library(Rglpk)
#noloci=3000
#noinds=2

#linkage group 'beginning' loci
#lgend=c(1,300,600,900,1200,1500,1800,2100,2400,2700)


###################################################
#Building Constraint Matrices
###################################################

#Simulated marker effect matrix
#effects=rnorm(noloci,2,0.01)
#geno=array(sample(c(0,1),noloci*2*noinds,replace=TRUE),dim=c(noloci,2,noinds))
#efmat=array(rep(effects,noinds),dim=c(noloci,2,noinds))

effects=marker.effects[,2]
geno=pars[,2:5]
noloci=nrow(pars)
noinds=2
#Marker effect genotypic matrix
data=geno*effects
#data=cbind(data[,,1],data[,,2])
map.tog=map.tog[which(map.tog$Marker%in%pars[,1]),]
map.tog=map.tog[match(map.tog[,1],pars[,1]),]

map.tog$Chrom=gsub(" ","",map.tog$Chrom)
e1<-max(which(as.character(map.tog$Chrom)==1))
e2<-max(which(as.character(map.tog$Chrom)==2))
e3<-max(which(as.character(map.tog$Chrom)==3))
e4<-max(which(as.character(map.tog$Chrom)==4))
e5<-max(which(as.character(map.tog$Chrom)==5))
e6<-max(which(as.character(map.tog$Chrom)==6))
e7<-max(which(as.character(map.tog$Chrom)==7))
e8<-max(which(as.character(map.tog$Chrom)==8))
e9<-max(which(as.character(map.tog$Chrom)==9))
e10<-max(which(as.character(map.tog$Chrom)==10))
lgend<-c(1,e1,e2,e3,e4,e5,e6,e7,e8,e9)

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



#Object function coefficients take into account the state of allele (i.e. 1 or -1)
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
#Optimal solution of decison variables
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

#Optimal haplotype
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
#res[del]=0
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
cnt=cnt+sum(which(x[,i]==1)%in%lgend)
print(which(x[,i]==1)%in%lgend)
recpts=c(recpts,which(x[,i]==1))
}
targrec=recpts[which(!recpts%in%lgend)]

recpts=map.tog[targrec,]
recpts=recpts[order(recpts$Chrom),]

hap.inds=sort(c(targrec,lgend,nrow(data)))

lbc=matrix(hap.inds[1:2],ncol=2)
for(i in 2:(length(hap.inds)-1)){
a=hap.inds[i]+1
b=hap.inds[i+1]
lbc=rbind(lbc,c(a,b))
}

haps<-c()
chr<-c()
pos<-c()
#Getting individual haplotype values
for(i in 1:nrow(lbc)){
a<-colSums(data[lbc[i,1]:lbc[i,2],])
chr=c(chr,map.tog$Chrom[lbc[i,1]])
pos=c(pos,map.tog$Position[lbc[i,2]])
haps=rbind(haps,a)
}
haps=data.frame(haps)

hap2<-data.frame(chr,pos,cbind(haps,apply(haps,1,which.max)))
names(hap2)=c("Chr","Pos","AH1","AH2","BH1","BH4","Max")

############################################################################
#Sampling GEBV DH lines (F2 derived)
############################################################################
library(Rcpp)
library(RcppArmadillo)
sourceCpp("/Users/gnyku/Desktop/Introgression-tools/Meiosis.cpp")

#Setting initial recombination frequencies
rcr=c()
#Setting any recombination distance greater than 50cM to 50cM (unlinked)

#Setting chromosome ends to being unlinked with next chromosome
ens=c()
for(i in 1:10){
ens<-c(ens,max(which(map.tog$Chrom==i)))
pos<-map.tog$Position[which(map.tog$Chrom==i)]


beg=c(50)
	for(j in 2:length(pos)){
		beg=c(beg,pos[j]-pos[j-1])
	}
rcr=c(rcr,beg)
}
ens=ens[1:(length(ens)-1)]
rr=rcr/100
rr[ens]=0.5
ens=ens+1

rr[ens]=0.5

data=as.matrix(data)
oo=colSums(meiosis(10000,rr,data[,c(1,3)],data[,c(1,3)]))
hist(oo)

#################################################################
#################################################################
recpts.list=list()
sols<-c()
for(norec in 1:50){
const.rhs[2852]=norec
xa=rep("B",ncol(const))
system.time({result=Rglpk_solve_LP(OF, const, const.dir, const.rhs, types = xa, control = list("verbose" =
TRUE), max = TRUE)})

#####################################
x=c()
#Optimal solution of decison variables
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

#Optimal haplotype
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
#res[del]=0
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
cnt=cnt+sum(which(x[,i]==1)%in%lgend)
print(which(x[,i]==1)%in%lgend)
recpts=c(recpts,which(x[,i]==1))
}
targrec=recpts[which(!recpts%in%lgend)]

recpts=map.tog[targrec,]
recpts=recpts[order(recpts$Chrom),]

hap.inds=sort(c(targrec,lgend,nrow(data)))

lbc=matrix(hap.inds[1:2],ncol=2)
for(i in 2:(length(hap.inds)-1)){
a=hap.inds[i]+1
b=hap.inds[i+1]
lbc=rbind(lbc,c(a,b))
}

haps<-c()
chr<-c()
pos<-c()
#Getting individual haplotype values
for(i in 1:nrow(lbc)){
a<-colSums(data[lbc[i,1]:lbc[i,2],])
chr=c(chr,map.tog$Chrom[lbc[i,1]])
pos=c(pos,map.tog$Position[lbc[i,2]])
haps=rbind(haps,a)
}
haps=data.frame(haps)

hap2<-data.frame(chr,pos,cbind(haps,apply(haps,1,which.max)))
names(hap2)=c("Chr","Pos","AH1","AH2","BH1","BH4","Max")


#####################################
sols=c(sols,result$optimum)
recpts.list=c(recpts.list,list(recpts))
}


plot(1:length(sols),sols,xlab="Number of optimal recombinations",ylab="Sum of marker effects",pch=20)



d