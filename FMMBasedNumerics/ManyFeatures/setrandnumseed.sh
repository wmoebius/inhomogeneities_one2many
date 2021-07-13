#!/bin/bash

# create random number seeds, argument is number of random number seeds

for i in $( seq 1 $1 )
do
   # find 32bit unsigned integer random number
   # http://unix.stackexchange.com/questions/268952/using-dev-random-dev-urandom-to-generate-random-data
   # http://www.thegeekstuff.com/2012/08/od-command/
   randnum=$(od -vAn -N4 -tu4 < /dev/urandom)
   echo "${i} ${randnum}"
done
