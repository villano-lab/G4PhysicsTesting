#!/bin/sh

# use this one for standard beryllium
for i in $(seq 1 55) 
do
  for j in 1 2 3 4 5  
  do
    ./condorNeutReflect -n 1 -set 18 -src 88Y -d -6 -dr ${i} -dt ${j} -bepure -ngen 1000000 -exe
  done
done

# use this one for standard beryllium
#for i in $(seq 1 55) 
#do
#  for j in 1.17	2.34 3.51 4.69  
#  do
#    ./bsubNeutReflect -n 1 -set 13 -src 124Sb -d -3 -dr ${i} -dt ${j} -bepure -ngen 1000000 -exe
#  done
#done

# use this one for BeO 
#for i in $(seq 1 55) 
#do
#  for j in 2 4 6 8 
#  do
#    ./bsubNeutReflect -n 1 -set 13 -src 124Sb -d -3 -dr ${i} -dt ${j} -ngen 1000000 -exe
#  done
#done
