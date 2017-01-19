#!/bin/bash

for n in {1..27}; do
   echo $n
   ncra -O output/high/Basin"$n"_ssticetemp_fv0.9x1.25.nc high.nc
   ncra -O output/low/Basin"$n"_ssticetemp_fv0.9x1.25.nc low.nc
   ncdiff -O high.nc low.nc output/diffs/$n.nc
   rm high.nc low.nc
done
