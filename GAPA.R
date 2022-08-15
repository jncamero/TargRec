
library(GA)
library(Rcpp)

sourceCpp("/Users/gnyku/Desktop/Introgression-tools/Meiosis.cpp")
#1Sample
chrlen=201
nochr=2
nop=2

#Parents
parents=matrix(sample(c(1,0),chrlen*nochr*4,replace=TRUE,prob=c(0.5,0.5)),c(chrlen*nochr,nochr*nop))

#Creating true haplotype markers
#3 SNPs per haplotype
haplen=chrlen/3
haps=c(1,cumsum(rep(3,haplen)))

hapinds<-cbind(seq(1,max(haps),3),seq(3,max(haps),3))

ind=3

haplist=list()
for(i in 1:nrow(hapinds)){
curhap=parents[hapinds[i,1]:hapinds[i,2],]
haplist=c(haplist,list(curhap[,which.max(colSums(parents[hapinds[i,1]:hapinds[i,2],]))]))
}

epieffects=rnorm(nrow(hapinds),mean=3,sd=1)

hapsinpop<-function(pop,hapinds){		
	haplist=list()

	for(i in 1:nrow(hapinds)){
		hapsin=c()
		curhap=pop[hapinds[i,1]:hapinds[i,2],]
			for(j in 1:ncol(curhap)){
				x<-as.character(curhap[,j])
				y<-paste(x,collapse="")
				hapsin=c(hapsin,y)	
				}
				print(length(unique(hapsin)))
			haplist=c(haplist,list(unique(hapsin)))

	}
	return(haplist)
}

#Calculate genetic value of each plant
phencalc<-function(plant){
score=sum(plant)
epi=c()
for(i in 1:nrow(hapinds)){
curhap=plant[hapinds[i,1]:hapinds[i,2],]
ll=sum(rowSums(curhap==haplist[[i]])>=1)==nrow(curhap)
print(ll)
epi=c(epi,ll)
}
output=list(as.numeric(epi)%*%epieffects,score)
return(output)
}

rr=rep(0.005,nrow(parents))
rr[1]=0.5
rr[201]=0.5
pro=meiosis(1000,rr,parents[,1:2],parents[,3:4])