#! /bin/bash

#  "fjet1MassTrimmed   210   -20   400"
#  "fjet1Tau2/fjet1Tau1   100   0   1"
#  "fjet1MassSDbm1   210   -20   400"


vars=(
  "fjet1QGtagSub2   110   -1.1   1.1"
  "fjet1MassSDb2   210   -20   400"
  "fjet1QJetVol   60   -0.5   2.5"
  "fjet1QGtagSub1   110   -1.1   1.1"
  "fjet1MassPruned   210   -20   400"
  "fjet1Tau2   75   0   150"
  "fjet1Tau1   150   0   300"
  "fjet1PullAngle   140   -3.5   3.5"
  "fjet1Pull   50   0   0.01"
  "fjet1C2b1   100   0   1"
  "fjet1C2b0p5   100   0   1"
  "fjet1C2b0p2   100   0   1"
  "2*fjet1QGtagSub2+fjet1QGtagSub1   120   -2.5   3.5"
  "fjet1QGtag   60   -0.1   1.1"
  "fjet1MassSDb0   210   -20   400"
  "fjet1C2b2   100   0   1"
)

for var in "${vars[@]}";do
    cp bestConfig.txt testConfig.txt
    echo $var >> testConfig.txt
    root -q -l -b MutualInfoNVars_reducedBias.C+\(\"testConfig.txt\"\)
done
