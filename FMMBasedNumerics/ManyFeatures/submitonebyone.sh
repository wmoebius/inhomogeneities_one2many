#!/bin/bash

cd ManyFeatures_$1
echo "ManyFeatures_$1" >> submitted.log 2>&1
./submit.sh >> submitted.log 2>&1
cd ..
