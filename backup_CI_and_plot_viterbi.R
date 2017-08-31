#plotviterbi

install.packages("moveHMM")
require(moveHMM)

plotSat(subsample[1:300,],zoom = 9,states = testvit)


#tabelle konfidenzintervalle

fthis <- round(confints(k6,3,0.05),digits = 2)
for(i in seq(1,72,by=3)){
  print(c(fthis[i,],fthis[i+1,],fthis[i+2,]))
}

fthis <- fthis[31:72,]
fthis

for(i in seq(1,42,by=6)){
  print(c(fthis[i,],fthis[i+1,],fthis[i+2,],fthis[i+3,],fthis[i+4,],fthis[i+5,]))
}

