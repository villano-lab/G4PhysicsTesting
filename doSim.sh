#!/bin/sh

FILE=$1
NUM=$2

RFILE=`basename ${FILE}`

IDSTRING=`echo ${RFILE} |awk 'BEGIN{FS="."}{print $1;}'`

echo ${RFILE}
echo ${IDSTRING}

#run the appropriate simulation
NeutReflectometry -src cascade -d ${NUM} -cascadein ${FILE} -set 0 -ngen 8000000 -otype txt -usestringlabel ${IDSTRING} -nomac -umncascade

#make a text file for the escape info
ls /data/chocula/villaa/cascadeSimData/Geant4/NeutReflectDiag*_Sourcecascade_neutron_${IDSTRING}_000_*.txt |awk '{system("./hitEvents.sh "$0);}' > hitevents_${IDSTRING}.txt

#use the escape info to write the tree
root -l -b -q addEscapeInfo.C\(\"${FILE}\",\"hitevents_${IDSTRING}.txt\",true\); 

#remove the escape info file
#rm hitevents_${IDSTRING}.txt
