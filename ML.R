library(tools)
library(gtools)
library(entropy)
filepath = "/home/karthik/Documents/CSE674/Project/Reqset/Cursive/Cursive_all"
#filepath = "/home/karthik/Documents/CSE674/Project/Reqset/Handprint/Handprint_all"
setwd(filepath)
ldf <- list() # creates a list
listcsv <- dir(pattern = "*.txt") # creates the list of all the csv files in the directory
numfiles = length(listcsv)
testdata = c(0,0,0,0,0,0,0,0,0,0,0,0)
cnames = colnames(testdata)
for (k in 1:numfiles){
  filename = listcsv[k]
  temp=read.csv(file=filename,header=F)
  count = length(readLines(filename))
  temp1 = matrix(NA,(nrow(temp)-1),ncol(temp))
  temp1 = temp[2:nrow(temp),(ncol(temp)-12):(ncol(temp)-1)]
  testdata=rbind(testdata,temp1)
}
grade5hand = na.omit(testdata[2:nrow(testdata),])
notcleanedmat = as.matrix(grade5hand)
for (u in 12:1){
  colnames(notcleanedmat)[u] <- paste('D',u)
}
#save(grade5hand,file="grade4cursive.Rda")   #### Save Data Set

for (j in 1:ncol(grade5hand)) {
  grade5hand[,j][grade5hand[,j]==max(grade5hand[,j][grade5hand[,j]==max(grade5hand[,j])] )] <- ((max(grade5hand[,j][grade5hand[,j]!=max(grade5hand[,j])])))
  grade5hand[,j][grade5hand[,j]==-1] <- ((max(grade5hand[,j][grade5hand[,j]!=max(grade5hand[,j])])))
} 
#save(grade5hand,file="grade4handcleaned.Rda")   ### Data Cleaning
matrixis = as.matrix(grade5hand)
for (u in 12:1){
  colnames(matrixis)[u] <- paste('D',u)
}

entrop = matrix(NA,1,12)
for (g in 1:ncol(matrixis)){
  entrop[1,g] = entropy(matrixis[,g],method="ML")
}
for (u in 12:1){
  colnames(entrop)[u] <- paste('D',u)
}

meanmat = matrix(NA,1,12)
for (g in 1:ncol(matrixis)){
  meanmat[1,g] = mean(matrixis[,g])
}
for (u in 12:1){
  colnames(meanmat)[u] <- paste('D',u)
}

save(entrop,file="grade3hand_entropy.Rda")
save(meanmat,file="grade3hand_mean.Rda")

Chimat = matrix(0,12,12)
for (i in 1:12){
  for (j in 1:12){
    Xsq <- chisq.test(matrixis[,i],matrixis[,j])
    #Chimat[i,j] <- Xsq$p.value
    Chimat[i,j] <- Xsq$statistic
  }
}  ### Pearson's Chi-square Test
Chimat = round(Chimat,10)
for (u in 12:1){
  colnames(Chimat)[u] <- paste('D',u)
}

Chimaxmat = matrix(NA,12,12)
for (i in 1:11){
  for (j in (i+1):12){
    Chimaxmat[i,j] <- Chimat[i,j]
  }
}  
for (i in 1:12){
  for (j in 1:i){
    Chimaxmat[i,j] <- 0
  }
} ### Generate Chimat with only above diagonal elements
for (u in 12:1){
  colnames(Chimaxmat)[u] <- paste('D',u)
}

maxpairs <- which(Chimaxmat >= sort(Chimaxmat, decreasing=T)[10], arr.ind=TRUE)  ### Get top 10 elements row and col indices
for (i in 1:nrow(maxpairs)){
A=data.frame(x1=c(matrixis[,maxpairs[i,1]]),x2=c(matrixis[,maxpairs[i,2]])) 
AA=data.frame(xx1=c(matrixis[,maxpairs[i,2]]),xx2=c(matrixis[,maxpairs[i,1]]))  ### get pairs from original data
Node1 = maxpairs[i,1]
Node2 = maxpairs[i,2]
Nodee1 = maxpairs[i,2]
Nodee2 = maxpairs[i,1]
B = assign(paste("paircount", i, sep = ""), table(A$x1, A$x2))
B1 = assign(paste("paircountt", i, sep = ""), table(AA$xx1, AA$xx2))
C = assign(paste("pairdf", i, sep = ""),as.data.frame(B))
C1 = assign(paste("pairdff", i, sep = ""),as.data.frame(B1))
D = assign(paste("pairmat", i, sep = ""),data.matrix(C))
D1 = assign(paste("pairmatt", i, sep = ""),data.matrix(C1))
S = sum(D[,3])
SS = sum(D1[,3])
E = matrix(NA,nrow(D),(ncol(D)+4))
E1 = matrix(NA,nrow(D1),(ncol(D1)+4))
E[,1:3] = D
E1[,1:3] = D1
for (j in 1:nrow(D)){ 
E[j,4] <- (D[j,3]/S)
} 
for (h in 1:nrow(D1)){ 
  E1[h,4] <- (D1[h,3]/SS)
} ## Inserts column into E with Freq/Total value
E = cbind(E,Node1,Node2)
E1 = cbind(E1,Nodee1,Nodee2)
#print ("here")
F = E 
F1 = E1  ### Get frequency of different combinations and (freq/total_count) with node pairs
for (l in 1:nrow(E)){ 
  G = data.frame(x3=c(matrixis[,Node2]))
  Cnt = table(G$x3,G$x3)
  Cnt1 = as.data.frame(Cnt)
  Cnt2 = data.matrix(Cnt1)
  S1 = sum(Cnt2[,3])
  G1 = table(G$x3==E[l,2])
  G2 = as.data.frame(G1)
  countofval = S1-(G2[1,2])
  cpd = countofval/S1
  #print (cpd)
  F[l,5]=cpd 
  F[l,6]=(F[l,4]/F[l,5])
  F[l,7]=log10(F[l,6])
}

for (j in 1:nrow(E1)){ 
  GG = data.frame(x3=c(matrixis[,Nodee2]))
  Cntt = table(GG$x3,GG$x3)
  Cntt1 = as.data.frame(Cntt)
  Cntt2 = data.matrix(Cntt1)
  SS1 = sum(Cntt2[,3])
  GG1 = table(GG$x3==E1[j,2])
  GG2 = as.data.frame(GG1)
  countofvall = SS1-(GG2[1,2])
  cpdd = countofvall/SS1
  F1[j,5]=cpdd
  F1[j,6]=F1[j,4]/F1[j,5]
  F1[j,7]=log10(F1[j,6])
}
F = na.omit(F)
F = F[!rowSums(!is.finite(F)),]
F1 = na.omit(F1)
F1 = F1[!rowSums(!is.finite(F1)),]
assign(paste("log", i, sep = ""),F)
assign(paste("logg", i, sep = ""),F1)
assign(paste("freqmat", i, sep = ""),cbind(F[1,8:9],sum(F[,7])))
assign(paste("freqmatt", i, sep = ""),cbind(F1[1,8:9],sum(F1[,7])))
}

H = rbind(freqmat1,freqmat2,freqmat3,freqmat4,freqmat5,freqmat6,freqmat7,freqmat8,freqmat9,freqmat10)
HH = rbind(freqmatt1,freqmatt2,freqmatt3,freqmatt4,freqmatt5,freqmatt6,freqmatt7,freqmatt8,freqmatt9,freqmatt10)

colnames(H)[2] <- 'logloss1'
colnames(H)[1] <- 'nodepairs'
colnames(HH)[2] <- 'logloss2'
colnames(HH)[1] <- 'nodepairs'

save(H,file="grade3hand_logloss1.Rda")
save(HH,file="grade3hand_logloss2.Rda")

####### Markov network ###################
invmat = solve(Chimat)
inversemat = round(invmat,10)
inversemat1 = inversemat
inversemat2 = matrix(0,12,12)
for (i in 1:12){
  for (j in i:12){
 if (inversemat1[i,j] > 0.0003 && inversemat1[i,j] < 0.0005) {
    inversemat2[i,j] <- inversemat1[i,j]
} 
}
}  




