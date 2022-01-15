#!/bin/sh

FILE=$1

#cat ${FILE} |awk '{if(NR>2 && $7>0){array[$1]++}}END{for(i in array){ print i;}}' |sort -n 
cat ${FILE} |awk '{
  if(NR>2 && $7>0){
    array[$1]++;
  } 
  
  #if(NR>2 && $3==100001){
  if(NR>2){
    det[$1]=$2
   }
}END{
  for(i in det){ 
    printf("%d\t%d",i,det[i]); 
    if(i in array){
      printf("\t1\n");
    }
    else{
      printf("\t0\n");
    }
  }
}' |sort -n 

