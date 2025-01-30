#!/bin/bash
for i in {1..10..1}
do 
  echo "running $i"
  ./liver_spheroid $i
done
