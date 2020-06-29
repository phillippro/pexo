#!/bin/bash
#for target in HD100623 HD115617 HD191408 HD217987 HD102365 HD128621 HD191849 HD26965 HD10700 HD156026 HD202560 GL433 GL551 GL699
for target in HD100623 
do 
     echo $target
     R CMD BATCH --no-save --no-restore "--args -m emulate -c TR -t ../input/${target}.tim -p ../input/${target}.par -v BJDtdb,RvRes -o ../results/${target}_PexoBary0.txt" pexo.R test.out &
    echo R CMD BATCH --no-save --no-restore "--args -m fit -c TR -t ../input/${target}.tim -p ../input/${target}.par -v 'BJDtdb RvRes' -o ../results/${target}_PexoBary0.txt" pexo.R test.out &
done 