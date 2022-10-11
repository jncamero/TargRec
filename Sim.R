#Reqs for Rgplk solver
#sudo apt-get install glpk-utils libglpk-dev glpk-doc

###########################################################
#####SAMPLING QTL positions, effects, and marker positions
###########################################################
library(rrBLUP)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("/Users/gnyku/Desktop/Introgression-tools/Meiosis.cpp")
sourceCpp("/mnt/Introgression-tools/Meiosis.cpp")
nochrom=10
#####
#heritability
h2<-0.25
# RR model for GS sim
nl<-100000   #no. hypothetical loci
no.qtl<-2000 #no. QTL
no.markers<-2000 #no. markers
#unique founder haplotypes
nohaps=10
#################################################################
#2 Draw marker, and QTL positions and sort
QTL<-sample(1:nl,no.qtl,replace=FALSE)
markers<-sample(1:nl,no.markers,replace=FALSE)
loci<-sort(unique(c(markers,QTL)))

###############################################################
#Indices for accessing 1) QTL and 2) marker loci 
#QTL and marker sets can overlap (i.e. intersection is not the null set)
qtl_index<-which(loci%in%QTL)
marker_index<-which(loci%in%markers)
#Sampline QTL effects
QTL_effects<-rnorm(no.qtl,0,0.1) 

################################################################
#Calculating genetic map-distances
################################################################
#1 Calculate subset all qtl and marker indicies into their respective linkage groups
#length of each chromosome
cl<-nl/10
lg<-seq(0,nl,cl) #linkage group intevals
#DISTANCE B/W hypothetical loci
key<-(100/cl)*0.01	#100 cM/ 200000 loci = 1 loci every 5*10^-4 cM, so c=5*10-4*0.01, because 1cM corresponds to c=0.01

#########################################################
#Calculating recombination rate between markers, and QTL
chr<-list()
for(i in 2:length(lg)){
m2<-sort(markers[markers>lg[i-1]&markers<=lg[i]])
qtl2<-sort(QTL[QTL>lg[i-1]&QTL<=lg[i]])
chr[[i-1]]=sort(unique(c(m2,qtl2)))
}

unqloc=length(unique(loci))
#1 Draw QTL positions (unique)
founder=array(rbinom(unqloc*nohaps,prob=0.5,size=1),dim=c(unqloc,nohaps))

#Recombination rate vector
rr_l<-list()
for(i in 1:nochrom){ #For i in 1:10 chromosome
y<-chr[[i]]
z<-y
z[1]=0
for(j in 2:length(y)){
z[j]<-(y[j]-y[j-1])*key	
}
#new linkage group
z[1]=0.5
rr_l[[i]]=z
}
##############################################
# Assembling recombination rate vector
rr<-c()
for(i in 1:nochrom){
rr<-c(rr,rr_l[[i]])
}

#########################################
g_val<-(t(founder[qtl_index,])%*%QTL_effects)

#Error sd to meet heritability requirements
error_sd<-sqrt((var(g_val)-(var(g_val)*h2))/h2)
phenotypes<-g_val + rnorm(ncol(founder),mean=0,sd=error_sd)
#3b from schematic: Phenotype 15,000 progeny,select best 200

#################################################
#Random mating to create linkage
#################################################
nogen=20
noprog=3000

F1=founder[,sample(1:10,100,replace=TRUE)]
cl<-cbind(sample(1:(ncol(founder)),noprog,replace=TRUE),sample(1:(ncol(founder)),noprog,replace=TRUE))
input=F1
for(i in 1:nogen){
	cl2=cl-1
	output=c()
	for(j in 1:nrow(cl)){
		a=input[,c(cl2[j,1]*2-1,cl2[j,1])]
		b=input[,c(cl2[j,2]*2-1,cl2[j,2])]
		output=cbind(output,meiosis(1,rr,a,b))
		print(dim(output))
	}
	cl<-cbind(sample(1:(ncol(output)/2),noprog,replace=TRUE),sample(1:(ncol(output)/2),noprog,replace=TRUE))
	input=output
	print(paste('generation',i))
	print(dim(input))
}
input1=input

###################################################
#Selection to deviate random allele frequencies
###################################################
nogen=10
#Number of individuals to select/ number of crosses
vb=noprog*.01
#Number of progen per cross
pny=ceiling(noprog/vb)

for(i in 1:nogen){
gval=t(getgeno(input[qtl_index,]))%*%QTL_effects
error_sd<-sqrt((var(g_val)/h2)-(var(g_val)))
phens=gval + rnorm(length(gval),0,error_sd)
ax=cbind(1:length(phens),phens)
ax2=ax[order(ax[,2],decreasing=TRUE),]
cl=array(cbind(sample(ax2[1:(vb*2),1],replace=FALSE)),dim=c(vb,2))
	output=c()
	for(j in 1:nrow(cl)){
		inds=sort(c(c(cl[1,1],cl[1,2])*2,c(cl[1,1],cl[1,2])*2-1))
		ab=input[,sort(inds)]
		output=cbind(output,meiosis(pny,rr,ab[,1:2],ab[,3:4]))
		print(c(i,j))
	}
	print(mean(gval))
	input=output
}
train=output
###################################################
#Creating Testing Set
###################################################
vb=10
pny=10

gval=t(getgeno(input[qtl_index,]))%*%QTL_effects
error_sd<-sqrt((var(g_val)/h2)-(var(g_val)))

phens=gval + rnorm(length(gval),mean=0,sd=error_sd)
ax=cbind(1:length(phens),phens)
ax2=ax[order(ax[,2],decreasing=TRUE),]
cl=array(cbind(sample(ax2[1:(vb*2),1],replace=FALSE)),dim=c(vb,2))

	output=c()
	for(j in 1:nrow(cl)){
		inds=sort(c(c(cl[1,1],cl[1,2])*2,c(cl[1,1],cl[1,2])*2-1))
		ab=input[,sort(inds)]
		output=cbind(output,meiosis(pny,rr,ab[,1:2],ab[,3:4]))
	}
test=output

##################################################
#Getting genotypes and phenotypes for training and testing in prediction model
##################################################
train.geno=c()
inds=cbind(seq(1,ncol(train),2),seq(2,ncol(train),2))
for(i in 1:nrow(inds)){
train.geno=cbind(train.geno,rowSums(train[,inds[i,]]))	
}

test.geno=c()

inds=cbind(seq(1,ncol(test),2),seq(2,ncol(test),2))
for(i in 1:nrow(inds)){
test.geno=cbind(test.geno,rowSums(test[,inds[i,]]))	
}

#########################################################
#Phenotypes
#########################################
g_val.train<-(t(train.geno[qtl_index,])%*%QTL_effects)
g_val.test<-(t(test.geno[qtl_index,])%*%QTL_effects)
g_tot=c(g_val.train,g_val.test)
#Error sd to meet heritability requirements
#error_sd<-sqrt((var(g_val)-(var(g_val)*h2))/h2)
error_sd<-sqrt((var(g_tot)/h2)-(var(g_tot)))
train.pheno<-g_val.train + rnorm(ncol(train.geno),mean=0,sd=error_sd)
test.pheno<-g_val.test + rnorm(ncol(test.geno),mean=0,sd=error_sd)

####################################################
#Training model
#####################################################
me=mixed.solve(train.pheno, Z=t(train.geno[marker_index,]))
#Estimated marker effects
meff=me$u

#Estimated 
ff=meff%*%test.geno[marker_index,]
cor(t(ff),test.pheno)

###################################################################################
###################################################################################
#Optimization
###################################################################################
###################################################################################
#Start for loop here outside of all possible combinations

#Parameters
norec=20

#Parental genotypes
a=test[,1:2]
b=test[,3:4]

cind=c(which(rr==0.5),nrow(test))

chr=c()
for(i in 1:(length(cind)-1)){
	if(i<(length(cind)-1)){
	chr=c(chr,rep(i,length(cind[i]:(cind[i+1]-1))))

	}else{
		chr=c(chr,rep(i,length(cind[i]:(cind[i+1]))))

	}
}

map.tog=data.frame(marker=c(1:length(chr)),Chrom=chr,Position=rr)
map.tog.markers=map.tog[marker_index,]

lgend1=c()
lgend2=c()
for(i in unique(map.tog$Chrom)){
lgend1=c(lgend1,min(which(map.tog.markers$Chrom==i)))
lgend2=c(lgend2,max(which(map.tog.markers$Chrom==i)))
}

pars=cbind(a,b)
#####################################################################################
###################################################
#Building Constraint Matrices
###################################################
effects=QTL_effects
geno=pars[marker_index,]
#noloci=nrow(pars)
noinds=2
#Estimated marker effects genotypic matrix
data=geno

for(i in 1:ncol(geno)){
data[,i]=geno[,i]*meff
}

data2=geno

for(i in 1:ncol(geno)){
data2[,i]=geno[,i]*QTL_effects
}

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
#Optimal solution of decison variables
#x is the optimal transition matrix
x=c()
res=result$solution
for(i in 1:nrow(data)){
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
	oo=lgend1*16-16+1
	tt=oo+15

	del=c()
	for(i in length(oo)){
	del=c(del,oo[i]:tt[i])
	}

res=result$solution
#res[del]=0
sum(res)

x=c()
for(i in 1:nrow(geno)){
x<-rbind(x,res[1:16])
res=res[-c(1:16)]
}
a<-rowSums(x[,1:4])
b<-rowSums(x[,5:8])
c<-rowSums(x[,9:12])
d<-rowSums(x[,13:16])
fin=cbind(a,b,c,d)

#return the Optimal recombination points
cnt=0
recpts=c()
for(i in c(2,3,4,6,7,8,10,11,12,14,15,16)){
cnt=cnt+sum(which(x[,i]==1)%in%lgend1)
print(which(x[,i]==1)%in%lgend1)
recpts=c(recpts,which(x[,i]==1))
}
targrec=recpts[which(!recpts%in%lgend1)]

recpts=map.tog.markers[targrec,]
recpts=recpts[order(recpts$marker),]

hap.inds=unique(sort(c(targrec,lgend1,lgend2)))

#BUILDING
#Chromosome and index position in map file of recombination points
zz=cbind(hap.inds,as.numeric(map.tog.markers$Chrom[hap.inds]))
output=hapint(zz,map.tog.markers,lgend1,lgend2,data)

#############################################################
#Loop
#############################################################
vals=c()
truval=c()
set=t(combn(100,2))
fina=list()
for(dp in 1:nrow(set)){

a=test[,(2*set[dp,1]-1):(2*set[dp,1])]
b=test[,(2*set[dp,2]-1):(2*set[dp,2])]
pars=cbind(a,b)
geno=pars[marker_index,]
####################################################
#Estimated marker effects x genotypic matrix
data=geno

for(i in 1:ncol(geno)){
data[,i]=geno[,i]*meff
}

OF=c()
for(i in 1:nrow(data)){
	for(j in 1:ncol(data)){
			OF=c(OF,rep(data[i,j],4))
		}
	}

xa=rep("B",ncol(const))
system.time({result=Rglpk_solve_LP(OF, const, const.dir, const.rhs, types = xa, control = list("verbose" =
TRUE), max = TRUE)})

x=c()
res=result$solution
for(i in 1:nrow(data)){
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
	oo=lgend1*16-16+1
	tt=oo+15

	del=c()
	for(i in length(oo)){
	del=c(del,oo[i]:tt[i])
	}

res=result$solution
#res[del]=0
sum(res)

x=c()
for(i in 1:nrow(geno)){
x<-rbind(x,res[1:16])
res=res[-c(1:16)]
}
a<-rowSums(x[,1:4])
b<-rowSums(x[,5:8])
c<-rowSums(x[,9:12])
d<-rowSums(x[,13:16])
fin=cbind(a,b,c,d)

#return the Optimal recombination points
cnt=0
recpts=c()
for(i in c(2,3,4,6,7,8,10,11,12,14,15,16)){
cnt=cnt+sum(which(x[,i]==1)%in%lgend1)
print(which(x[,i]==1)%in%lgend1)
recpts=c(recpts,which(x[,i]==1))
}
targrec=recpts[which(!recpts%in%lgend1)]

recpts=map.tog.markers[targrec,]
recpts=recpts[order(recpts$Chrom,recpts$marker),]

hap.inds=sort(c(targrec,lgend1,nrow(data)))
#BUILDING
#Chromosome and index position in map file of recombination points
zz=cbind(hap.inds,as.numeric(map.tog.markers$Chrom[hap.inds]))

output=hapint(zz,map.tog.markers,lgend1,lgend2,data)

#GEBV
b=sum(data*fin)
vals=b
write.csv(file='vals.out.csv',row.names=FALSE,c(a,b))
#tru<-truth(pars,map.tog,recpts)
#TRUE GENETIC VALUE
tru<-sum(truth())
fina=c(fina,list(fin))
truval=c(truval,tru)
}

#################################################
#Calculating GEBV value
lgenda=lgend1
lgendb=lgend2
out=c()
chr=c()
for(i in 1:nochrom){
doot=pars[which(map.tog.markers$Chrom==i),]
haps=subset(recpts,Chrom==i)$marker
ind=sort(unique(c(which(map.tog.markers$marker%in%haps),lgenda[i],lgendb[i])))
ind2=sort(c(ind,ind[2:(length(ind)-1)]+1))
print(length(ind2)/2)
nh=c()			
first=seq(1,length(ind2),2)
last=seq(1,length(ind2),2)+1
	for(j in 1:length(seq(1,length(ind2),2))){

			if(!is.null(colSums(data[ind2[first[j]]:ind2[last[j]],]))){
			print("YES")
			ans=colSums(data[ind2[first[j]]:ind2[last[j]],])
			}else{ans=data[ind2[first[j]]:ind2[last[j]],]
			}
			nh=rbind(nh,ans)
			#chr=c(chr,rep(i,length())))
		}

out=rbind(out,nh)
}

#########################################################################
#########################################################################
#########################################################################
truth=function(){
lgenda=c()
lgendb=c()
for(i in unique(map.tog$Chrom)){
lgenda=c(lgenda,min(which(map.tog$Chrom==i)))
lgendb=c(lgendb,max(which(map.tog$Chrom==i)))
}

#Calculating True Genetic value
out=c()
chr=c()
qtlcnt=1
for(i in 1:nochrom){

doot=pars[which(map.tog$Chrom==i),]
haps=subset(recpts,Chrom==i)$marker
#ind=sort(unique(c(which(map.tog.markers$marker%in%haps),lgenda[i],lgendb[i])))
#Get loci interval of QTL in haplotype
ind=sort(c(map.tog.markers[c(which(map.tog.markers$marker%in%haps)),]$marker,lgenda[i],lgendb[i]))

#If there are more than 2 QTL in the haplotype
if(length(ind)>2){
ind2=sort(c(ind,ind[2:(length(ind)-1)]-1))
}else{ind2=ind}
print(length(ind2)/2)
nh=c()			
first=seq(1,length(ind2),2)
last=seq(1,length(ind2),2)+1

	for(j in 1:length(seq(1,length(ind2),2))){
		rws=qtl_index[which(qtl_index>=ind2[first[j]]&qtl_index<=ind2[last[j]])]
		
		if(length(rws)==1){
		hapvals=t(t(pars[rws,])*QTL_effects[which(qtl_index>=min(rws)&qtl_index<=max(rws))])
		nh=rbind(nh,t(hapvals))
		}else if(length(rws)>=2){
		hapvals=t(t(pars[rws,])%*% QTL_effects[which(qtl_index>=min(rws)&qtl_index<=max(rws))])
		nh=rbind(nh,hapvals)
			}else if(length(rws)==0){
		nh=rbind(nh,rep(0,4))
		}
		}
	out=rbind(out,nh)
	chr=c(chr,rep(i,length(ind2)/2))			
	}

ll=output[[1]]

gh=c()
for(i in 1:nrow(ll)){
gh=c(gh,out[i,ll$Max[i]])
}
return(gh)
}
gh=truth()

###########################################################################
###########################################################################
#############PARALLEL SIMULATION CODE######################################
###########################################################################
###########################################################################
library(parallel)
library(snow)
library(tictoc)
detectCores()

truval=c()
set=t(combn(100,2))
fina=list()

cty<- makeCluster(20, type = "SOCK")
clusterEvalQ(cty, library(Rglpk))
clusterExport(cty, list("nochrom","truth","hapint","getgeno","tmat","map.tog",
"map.tog.markers","test","qtl_index","QTL_effects",
"marker_index","meff","effects","set",
"const","const.rhs","const.dir",
"lgend1","lgend2","opt"))

opt=function(dp){
a<-test[,(2*set[dp,1]-1):(2*set[dp,1])]
b<-test[,(2*set[dp,2]-1):(2*set[dp,2])]
pars<-cbind(a,b)
geno<-pars[marker_index,]
####################################################
#Estimated marker effects x genotypic matrix
data<-geno

for(i in 1:ncol(geno)){
data[,i]=geno[,i]*meff
}

OF=c()
for(i in 1:nrow(data)){
	for(j in 1:ncol(data)){
			OF=c(OF,rep(data[i,j],4))
		}
	}

xa=rep("B",ncol(const))
system.time({result=Rglpk_solve_LP(OF, const, const.dir, const.rhs, types = xa, control = list("verbose" =
TRUE), max = TRUE)})

x=c()
res=result$solution
for(i in 1:nrow(data)){
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
	oo<-lgend1*16-16+1
	tt<-oo+15

	del=c()
	for(i in length(oo)){
	del=c(del,oo[i]:tt[i])
	}

res=result$solution
#res[del]=0
sum(res)

x=c()
for(i in 1:nrow(geno)){
x<-rbind(x,res[1:16])
res=res[-c(1:16)]
}
a<-rowSums(x[,1:4])
b<-rowSums(x[,5:8])
c<-rowSums(x[,9:12])
d<-rowSums(x[,13:16])
fin=cbind(a,b,c,d)

#return the Optimal recombination points
cnt=0
recpts=c()
for(i in c(2,3,4,6,7,8,10,11,12,14,15,16)){
cnt=cnt+sum(which(x[,i]==1)%in%lgend1)
print(which(x[,i]==1)%in%lgend1)
recpts=c(recpts,which(x[,i]==1))
}
targrec=recpts[which(!recpts%in%lgend1)]

recpts=map.tog.markers[targrec,]
recpts=recpts[order(recpts$Chrom,recpts$marker),]

hap.inds=sort(c(targrec,lgend1,nrow(data)))
#BUILDING
#Chromosome and index position in map file of recombination points
zz=cbind(hap.inds,as.numeric(map.tog.markers$Chrom[hap.inds]))

output=hapint(zz,map.tog.markers,lgend1,lgend2,data)

#GEBV
val=sum(data*fin)
#write.csv(file='vals.out.csv',row.names=FALSE,c(a,b))
#tru<-truth(pars,map.tog,recpts)
#TRUTH
tru<-sum(truth(output,map.tog,map.tog.markers,recpts,pars,qtl_index,QTL_effects,nochrom))
return(list(val,tru))
}

tic();
output=parSapply(cty,1:nrow(set),opt);
toc();

stopCluster(cty)
#####################################################################
#####################################################################
#####################################################################
#FUNCTIONS
#####################################################################
#####################################################################
#####################################################################

truth=function(output,map.tog,map.tog.markers,recpts,pars,qtl_index,QTL_effects,nochrom){
lgenda=c()
lgendb=c()
for(i in unique(map.tog$Chrom)){
lgenda=c(lgenda,min(which(map.tog$Chrom==i)))
lgendb=c(lgendb,max(which(map.tog$Chrom==i)))
}

#Calculating True Genetic value
out=c()
chr=c()
for(i in 1:nochrom){

doot=pars[which(map.tog$Chrom==i),]
haps=subset(recpts,Chrom==i)$marker
#ind=sort(unique(c(which(map.tog.markers$marker%in%haps),lgenda[i],lgendb[i])))
#Get loci interval of QTL in haplotype
ind=sort(c(map.tog.markers[c(which(map.tog.markers$marker%in%haps)),]$marker,lgenda[i],lgendb[i]))

#If there are more than 2 QTL in the haplotype
if(length(ind)>2){
ind2=sort(c(ind,ind[2:(length(ind)-1)]-1))
}else{ind2=ind}
print(length(ind2)/2)
nh=c()			
first=seq(1,length(ind2),2)
last=seq(1,length(ind2),2)+1

	for(j in 1:length(seq(1,length(ind2),2))){
		rws=qtl_index[which(qtl_index>=ind2[first[j]]&qtl_index<=ind2[last[j]])]
		
		if(length(rws)==1){
		hapvals=t(t(pars[rws,])*QTL_effects[which(qtl_index>=min(rws)&qtl_index<=max(rws))])
		nh=rbind(nh,t(hapvals))
		}else if(length(rws)>=2){
		hapvals=t(t(pars[rws,])%*% QTL_effects[which(qtl_index>=min(rws)&qtl_index<=max(rws))])
		nh=rbind(nh,hapvals)
			}else if(length(rws)==0){
		nh=rbind(nh,rep(0,4))
		}
		}
	out=rbind(out,nh)
	chr=c(chr,rep(i,length(ind2)/2))			
	}

ll=output[[1]]

gh=c()
for(i in 1:nrow(ll)){
gh=c(gh,out[i,ll$Max[i]])
}
return(gh)
}
gh=truth()

hapint=function(zz,map.tog,lgend1,lgend2,data){
chroms=c(1:10)
lbc1=c()
lbc2=c()
for(j in chroms){

	waw=subset(zz,zz[,2]==j)[,1]
	if(length(waw)>1){
	inds=unique(sort(c(waw,lgend2[j])))
	ss=sort(c(inds,inds[-c(1,length(inds))]+1))
	bld1=c()
	bld2=c()
	for(i in 1:(length(ss)/2)){
		bld1<-rbind(bld1,map.tog[ss[1:2],]$Position)
		bld2<-rbind(bld2,ss[1:2])
		ss=ss[-c(1:2)]
	}
	}else{
		a<-subset(map.tog,Chrom==j)
		a$Position=cumsum(a$Position)
		mi=min(a$Position)
		ma=max(a$Position)
		fi=min(which(map.tog$Chrom==j))
		la=max(which(map.tog$Chrom==j))
		bld1<-c(mi,ma)
		bld2<-c(fi,la)
	}
print(paste("CHROMOSOME",j))
print(bld2)
print(lbc2)
lbc1=rbind(lbc1,bld1)
lbc2=rbind(lbc2,bld2)

haps<-c()
chr<-c()
pos<-c()
#Getting individual haplotype values
for(i in 1:nrow(lbc2)){
	if(length(lbc2[i,1]:lbc2[i,2])>1){
a<-colSums(data[lbc2[i,1]:lbc2[i,2],])
}else{
	a<-data[lbc2[i,1]:lbc2[i,2],]
	}
chr=c(chr,map.tog$Chrom[lbc2[i,1]])
#pos=c(pos,map.tog$Position[lbc2[i,2]])
haps=rbind(haps,a)
print(i)
}
haps=data.frame(haps)

hap2<-data.frame(chr,cbind(haps,apply(haps,1,which.max)))
names(hap2)=c("Chr","AH1","AH2","BH1","BH4","Max")

stop=map.tog.markers$Position[lbc2[,2]]
start=map.tog.markers$Position[lbc2[,1]]
start[1]=1
hap2$start=start
hap2$stop=stop
}
return(list(hap2))
}

getgeno=function(test){
test.geno=c()

inds=cbind(seq(1,ncol(test),2),seq(2,ncol(test),2))
for(i in 1:nrow(inds)){
test.geno=cbind(test.geno,rowSums(test[,inds[i,]]))	
}
return(test.geno)
}

tmat<-function(one,two,three,four){
#Constraint: Inflow must equal outflow
	output=array(0,dim=c(4,nrow(geno)*16))
	output[1,1:32]=one
	output[2,1:32]=two
	output[3,1:32]=three
	output[4,1:32]=four
	write.table(file='output.csv',sep=",",col.names=FALSE,row.names=FALSE,output)

	op2=list()
	for(i in 2:(nrow(geno)-1)){
		print(i)
		start.index=16*(i-1)+1
		hold=array(0,dim=c(4,nrow(geno)*16))
		hold[1,(start.index):(start.index+32-1)]=one
		hold[2,(start.index):(start.index+32-1)]=two
		hold[3,(start.index):(start.index+32-1)]=three
		hold[4,(start.index):(start.index+32-1)]=four
		#output<-rbind(output,hold)
		write.table(file='output.csv',sep=",",append=TRUE,col.names=FALSE,row.names=FALSE,hold)
	}
	output<-read.csv('output.csv',header=FALSE)
	oo=lgend1*16-16+1
	tt=oo+15

#Constraint: that only 1 'decision' can occur at a locus (i.e. recombine vs. not recombine)
		tm2<-array(0,dim=c(ncol(output)/16,nrow(geno)*16))
		tmzeros<-array(0,dim=c(ncol(output)/16,nrow(geno)*16))
		a2=seq(1,ncol(output),16)
		b2=seq(16,ncol(output),16)
		
		#Don't start on transition state at "1st locus"
		tm2[1,c(1,5,9,13)]=1

		for(i in 2:(nrow(geno))){
			#If no loci upstream on linkage group
			tm2[i,(16*(i-1)+1):(16*i)]=1
		}

#Targeted recombination constraint: no more than norec recombinations
#Accounting for linkage group ends
	trc=rep(c(0,1,1,1),4*nrow(geno))
			for(i in 1:length(oo)){
				trc[oo[i]:tt[i]]=0
			}
#LHS
lhs=rbind(output,tm2,trc)
sign=c(rep("==",nrow(output)),rep("==",nrow(tm2)),"<=")
rhs=c(rep(0,nrow(output)),rep(1,nrow(tm2)),norec)
	return(list(lhs,sign,rhs))
}
################################################
#END FUNCTIONS
################################################
nn=array(0,length(rr))
nn[qtl_index]=QTL_effects
cc=cbind(map.tog$Chrom,nn)

qtl=map.tog[qtl_index,]
hh=subset(qtl,Chrom==1)
p2=t(pars[which(map.tog$Chrom==10),])%*%nn[which(map.tog$Chrom==10)]
p3=t(pars[which(map.tog$Chrom==1),])%*%nn[which(map.tog$Chrom==10)]

t(getgen(meiosis(100,rr,a,b)))%*%QTL_effects
