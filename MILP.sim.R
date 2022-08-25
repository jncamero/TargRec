library(lpSolve)
noloci=3
noinds=5
#linkage group 'endpoint' loci
lgend=c(1)

#Number of feasible recombinations
norec=2

###################################################
#Simulated marker effect matrix
effects=rnorm(noloci,2,0.01)
geno=array(sample(c(0,1),noloci*2*noinds,replace=TRUE),dim=c(noloci,2,noinds))
efmat=array(rep(effects,noinds),dim=c(noloci,2,noinds))

#Marker effect genotypic matrix
data=geno*effects

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

		cons<-rbind(cons,array(0,ncol(cons)))
		conszeros<-array(0,dim=c(nrow(cons),ncol(cons)))
		selzeros=array(0,dim=c(noinds+1,noinds))

		sels=-diag(noinds)*(noloci*2)

		#sels=rbind(sels,1)
		#LHS
		sels2=rbind(sels,1)
		##########################################################
		#Constraint variables for selecting only 1 individual for parent 1
		#Constraint variables for selecting only 1 individual for parent 1
		#Just integer variables, no marker effects
		#[data variable for Parent 1][data variable for Parent 2][S variable][Zeros]
		#[data variable for Parent 1][data variable for Parent 2][S variable][Zeros]
		a1=cbind(cons,conszeros,sels2,selzeros)
		a2=cbind(conszeros,cons,selzeros,sels2)
		lhs=rbind(a1,a2)
		###################################################
		#RHS
		#
		rhs.half=c(rep(0,nrow(cons)-1),1)
		rhs=rep(rhs.half,2)
		signs1=c(rep("=",(nrow(cons)-1)),"=")
		signs2=c(rep("=",(nrow(cons)-1)),"=")
		signs=c(signs1,signs2)

		#nrow is 2*noinds
		return(list(lhs,signs,rhs))
	}
#Output is:
#1 LHS (data for selection of 2 parents)
#2 Signs of constraints
#3 RHS 

dc<-datacons(noloci,noinds)
###########################################################
#Marker effects as coefficients on data decision variable matrix
#Inputs: number of loci, number of individuals to select from, marker effects
#X<=sumk(dataijk)
xmat=function(noloci,noinds,data){
		#X matrix marker effects for P1 and P2
		#LHS
		#[Marker effect variables for P1][]
	it=array(0,dim=c(noloci*2,noloci*2*noinds))
		for(i in 1:(noloci*2)){
			if(i<=noloci){
				it[i,seq(i,noloci*2*noinds,noloci*2)]=data[i,1,]
			}
					if(i>noloci){
						print(i)
						ah=i-noloci
							it[i,seq(i,noloci*2*noinds,noloci*2)]=data[ah,2,]

				}
		}
		#LHS variable structure
		#[Parent1][Parent2][S1][S2][X1][X2]
		#[Parent1][Parent2][S1][S2][X1][X2]
		#RHS
			dat0<-array(0,dim=c(nrow(it),ncol(it)))
			s=array(0,dim=c(noinds*2,noinds*2))
			x=-diag(noloci*2)
			x0=array(0,dim=c(noloci*2,noloci*2))			
			

			z=rbind(cbind(x,x0),cbind(x0,x))
			it2=rbind(cbind(it,dat0),cbind(dat0,it))
			#S variable from previous function, picking 2 parents
			s=array(0,dim=c(nrow(it2),noinds*2))
			#LHS variables: Data, S, 
			lhs=cbind(it2,s,z)
			sign=rep("=",nrow(lhs))
			rhs=rep(0,nrow(lhs))
			return(list(lhs,sign,rhs))
	}

xm<-xmat(noloci,noinds,data)
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

	op2=list()
	for(i in 2:(noloci-1)){
		start.index=16*(i-1)+1
		hold=array(0,dim=c(4,noloci*16))
		hold[1,(start.index):(start.index+32-1)]=one
		hold[2,(start.index):(start.index+32-1)]=two
		hold[3,(start.index):(start.index+32-1)]=three
		hold[4,(start.index):(start.index+32-1)]=four
		output<-rbind(output,hold)
	
		}

	oo=lgend*16-16+1
	tt=oo+15

#Constraint: that only 1 'decision' can occur at a locus (i.e. recombine vs. not recombine)
		tm2<-array(0,dim=c(ncol(output)/16,noloci*16))
		tmzeros<-array(0,dim=c(ncol(output)/16,noloci*16))
		a2=seq(1,ncol(output),16)
		b2=seq(16,ncol(output),16)
		
		for(i in 1:(noloci)){
			#If no loci upstream on linkage group, f
			if(i%in%lgend){
				tm2[i,c(1,5,9,13)]=1
				}else{
			tm2[i,a2[i]:b2[i]]=1
		}
	}

#Targeted recombination constraint: no more than norec recombinations
#Accounting for linkage group ends
	trc=rep(c(0,1,1,1),4*noloci)
			for(i in 1:length(oo)){
				trc[oo[i]:tt[i]]=0
			}

#LHS
lhs=rbind(output,tm2,trc)
sign=c(rep("=",nrow(output)),rep("=",nrow(tm2)),"<=")
rhs=c(rep(0,nrow(output)),rep(1,nrow(tm2)),norec)
	return(list(lhs,sign,rhs))
}


tm=tmat(one,two,three,four)

############################################
#Constraints on value of X by T
xt<-function(tm,xm){
xx=diag(noloci*4)
tt=array(0,dim=c(4,16*noloci))
#First locus of each linkage group
tt[1,1]=-1
tt[2,5]=-1
tt[3,9]=-1
tt[4,13]=-1

zz=array(0,dim=c(nrow(xx),16*noloci))

for(i in 5:nrow(zz)){
zz[i,((i-1)*4+1):(i*4)]=-1
}

zz[1:4,]=tt[1:4,]
#ph: [Parent 1][Parent 2][S1][S2]
ph=array(0,dim=c(nrow(xx),2*(2*noinds*noloci+noinds)))
lhs=cbind(ph,xx,zz)
rhs=rep(0,nrow(lhs))
sign=rep("<=",nrow(lhs))

return(list(lhs,sign,rhs))
}

gg=xt(tm,xm)

########################################
#Compiling variables: LHS, SIGN, RHS
#dc sets the value of xij as sum of coefficients (marker-effect) at

compile<-function(dc,xm,tm,gg){
a<-ncol(gg[[1]])-ncol(dc[[1]])
b<-ncol(gg[[1]])-ncol(xm[[1]])
c<-ncol(gg[[1]])-ncol(tm[[1]])

#Variable block matrix structure
#[Meff Data for parent 1][S variable][Zeros]
#[Meff Data for Parent 2][Zeros][S variable]
datamat=cbind(dc[[1]],array(0,dim=c(nrow(dc[[1]]),a)))
xm2=cbind(xm[[1]],array(0,dim=c(nrow(xm[[1]]),b)))
tm2=cbind(array(0,dim=c(nrow(tm[[1]]),c)),tm[[1]])

lhs.big=rbind(datamat,xm2,tm2,gg[[1]])
sign.big=c(dc[[2]],xm[[2]],tm[[2]],gg[[2]])
rhs.big=c(dc[[3]],xm[[3]],tm[[3]],gg[[3]])
return(list(lhs.big,sign.big,rhs.big))
}

fin=compile(dc,xm,tm,gg)
const.rhs=fin[[3]]
const.dir=fin[[2]]
const=fin[[1]]

OF=array(0,ncol(fin[[1]]))
a=(ncol(fin[[1]])+1)-(noloci*16+noloci*2*2)
b=a:(a+noloci*2*2-1)
OF[b]=1

bin=which(OF!=1)

result=lp(direction = "max", OF, const, const.dir, const.rhs,
    transpose.constraints = TRUE, binary.vec=bin)