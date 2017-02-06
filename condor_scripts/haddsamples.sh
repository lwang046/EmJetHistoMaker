#!/bin/env bash
OUTPUTDIR=$1
echo $OUTPUTDIR
for d in QCD_HT500to700 QCD_HT700to1000 QCD_HT1000to1500 QCD_HT1500to2000 QCD_HT2000toInf
do
    echo "./haddnorm/haddnorm eventCountPreTrigger histo$OUTPUTDIR/$d.root"
done
