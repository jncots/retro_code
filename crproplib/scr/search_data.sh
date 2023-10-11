#!/bin/bash
(set -o igncr) 2>/dev/null && set -o igncr; # this comment is needed for Windows system

dir_search=$1
cd "$dir_search"


calc_dir0=($(echo */))

if [[ "$calc_dir0" =~ "*/" ]]; then
   calc_dir=("$dir_search")
else
  declare -a calc_dir
 i=0
 for calc in "${calc_dir0[@]}"
 do
  calc_dir[$i]="$dir_search"$calc
  i=$(( i + 1))
 done
 calc_dir=("$dir_search" "${calc_dir[@]}") 
fi

mpref=$2
declare -a data_files # array of files with data

i=0
for calc in "${calc_dir[@]}"
do
# echo 'calc=' "$calc"
 cd "$calc"
 calc_file=($(ls)) 
 for cf in "${calc_file[@]}"
 do
 if [[ "$cf" =~ "$mpref" ]]; then
   i=$(( i + 1))
   path="$(pwd)/$cf"
   
   if [[ "$path" =~ "cygdrive" ]]
    then
      path="$(cygpath -w "$path")"
   fi
   data_files[$i]="$path"
  fi
 done
done

for arg in "${data_files[@]}"
do
 echo "$arg">>$3
done
