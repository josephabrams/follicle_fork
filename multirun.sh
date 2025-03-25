#!/bin/bash
for i in {1..3..1}
do 
  echo "running $i"
  ./liver_spheroid 1 2 3 $i
done
