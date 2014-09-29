#! /bin/bash

root -q -l -b MutualInfoSingleVar.C+\(\"2*fjet1QGtagSub2+fjet1QGtagSub1\",\"fjet1QGTagComb\",120,-2.5,3.5\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1QGtagSub1\",\"fjet1QGtagSub1\",110,-1.1,1.1\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1QGtagSub2\",\"fjet1QGtagSub2\",110,-1.1,1.1\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1QGtag\",\"fjet1QGtag\",60,-0.1,1.1\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1PullAngle\",\"fjet1PullAngle\",140,-3.5,3.5\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1Pull\",\"fjet1Pull\",50,0,0.01\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1MassTrimmed\",\"fjet1MassTrimmed\",210,-20,400\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1MassPruned\",\"fjet1MassPruned\",210,-20,400\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1MassSDbm1\",\"fjet1MassSDbm1\",210,-20,400\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1MassSDb2\",\"fjet1MassSDb2\",210,-20,400\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1MassSDb0\",\"fjet1MassSDb0\",210,-20,400\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1QJetVol\",\"fjet1QJetVol\",60,-0.5,2.5\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1C2b2\",\"fjet1C2b2\",100,0,1\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1C2b1\",\"fjet1C2b1\",100,0,1\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1C2b0p5\",\"fjet1C2b0p5\",100,0,1\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1C2b0p2\",\"fjet1C2b0p2\",100,0,1\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1Tau2\",\"fjet1Tau2\",75,0,150\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1Tau1\",\"fjet1Tau1\",150,0,300\)
root -q -l -b MutualInfoSingleVar.C+\(\"fjet1Tau2/fjet1Tau1\",\"fjet1Tau2Tau1\",100,0,1\)
