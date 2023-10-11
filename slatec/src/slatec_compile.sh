#! /bin/bash

for f in *.o; do
 rm "$f"
done

for f in *.mod; do
 rm "$f"
done

gfortran *.f -W -fbackslash -c

for f in *.o; do
 ar cr libslatec.a "$f"
done