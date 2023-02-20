#hmmbuilder
HMMEmissionCounter<-function(corematrix)
{
#instantiate hmm dataframe
  state<-c('Core','Miscore','ComplexCore','Spacer','Noise')
  stateInstances<-c(0,0,0,0,0)
  CoreEmission<-c(0,0,0,0,0)
  MisCoreEmission<-c(0,0,0,0,0)
  ComplexCoreEmission<-c(0,0,0,0,0)
  SpacerEmission<-c(0,0,0,0,0)
  NoiseEmission<-c(0,0,0,0,0)
  
  
  hmm <- data.frame(state,stateInstances,CoreEmission,MisCoreEmission,ComplexCoreEmission,SpacerEmission,NoiseEmission)
  


#loop through and for each look at the succeding entry and fill in hmm dataframe

for (i in 1:nrow(corematrix)-1)
{

  #in the future take into account coresize
coresize<-(corematrix[i,][,'start']-corematrix[i,][,'end'])/4

corestate<-toString(corematrix[i,][,'state'])
emittedstate<-toString(corematrix[i+1,][,'state'])




if (corestate=='core')
{
hmm[1,][2]=hmm[1,][2]+1


if (emittedstate=='core')
hmm[1,][3]=hmm[1,][3]+1
else if (emittedstate=='miscore' || emittedstate=='miscoregroup')
hmm[1,][4]=hmm[1,][4]+1
else if (emittedstate=='complexcore')
hmm[1,][5]=hmm[1,][5]+1
else if (emittedstate=='spacer')
hmm[1,][6]=hmm[1,][6]+1
else if (emittedstate=='noise')
hmm[1,][7]=hmm[1,][7]+1
}

else if(corestate=='miscore' || corestate=='miscoregroup')
{
hmm[2,][2]=hmm[2,][2]+1

if (emittedstate=='core')
hmm[2,][3]=hmm[2,][3]+1
else if (emittedstate=='miscore' || emittedstate=='miscoregroup')
hmm[2,][4]=hmm[2,][4]+1
else if (emittedstate=='complexcore')
hmm[2,][5]=hmm[2,][5]+1
else if (emittedstate=='spacer')
hmm[2,][6]=hmm[2,][6]+1
else if (emittedstate=='noise')
hmm[2,][7]=hmm[2,][7]+1
}


else if (corestate=='complexcore')
{
hmm[3,][2]=hmm[3,][2]+1

if (emittedstate=='core')
hmm[3,][3]=hmm[3,][3]+1
else if (emittedstate=='miscore' || emittedstate=='miscoregroup')
hmm[3,][4]=hmm[3,][4]+1
else if (emittedstate=='complexcore')
hmm[3,][5]=hmm[3,][5]+1
else if (emittedstate=='spacer')
hmm[3,][6]=hmm[3,][6]+1
else if (emittedstate=='noise')
hmm[3,][7]=hmm[3,][7]+1
}

else if (corestate=='spacer')
{
hmm[4,][2]=hmm[4,][2]+1

if (emittedstate=='core')
hmm[4,][3]=hmm[4,][3]+1
else if (emittedstate=='miscore' || emittedstate=='miscoregroup')
hmm[4,][4]=hmm[4,][4]+1
else if (emittedstate=='complexcore')
hmm[4,][5]=hmm[4,][5]+1
else if (emittedstate=='spacer')
hmm[4,][6]=hmm[4,][6]+1
else if (emittedstate=='noise')
hmm[4,][7]=hmm[4,][7]+1
}



else if (corestate=='noise')
{
hmm[5,][2]=hmm[5,][2]+1

if (emittedstate=='core')
hmm[5,][3]=hmm[5,][3]+1
else if (emittedstate=='miscore' || emittedstate=='miscoregroup')
hmm[5,][4]=hmm[5,][4]+1
else if (emittedstate=='complexcore')
hmm[5,][5]=hmm[5,][5]+1
else if (emittedstate=='spacer')
hmm[5,][6]=hmm[5,][6]+1
else if (emittedstate=='noise')
hmm[5,][7]=hmm[5,][7]+1
}



}

hmm

}

