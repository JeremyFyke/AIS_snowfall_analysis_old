#!/bin/bash

for c in low mean high; do
    p=output/$c
    for m in $(seq -f "%02g" 1 12); do
       ncra -O -F -d time,$m,,12 -v ICEFRAC $p/ssticetemp_fv0.9x1.25.nc $p/TM_"$m".nc
       ncra -A -F -d time,$m,,12 -v SST     $p/ssticetemp_fv0.9x1.25.nc $p/TM_"$m".nc
       ncra -A -F -d time,$m,,12 -v PRECIP  $p/ssticetemp_fv0.9x1.25.nc $p/TM_"$m".nc       
    done
    ncrcat -O $p/TM_*.nc $p/TM.nc
    rm $p/TM_*.nc
done
