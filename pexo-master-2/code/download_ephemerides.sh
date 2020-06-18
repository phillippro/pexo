#!/bin/bash
DE=de$1
newDir=../data/
wget -r --no-parent ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/$DE/
#mkdir $newDir
mv ssd.jpl.nasa.gov/pub/eph/planets/ascii/$DE ../data/
rm -r ssd.jpl.nasa.gov

