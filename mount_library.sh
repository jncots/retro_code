#!/bin/bash
(set -o igncr) 2>/dev/null && set -o igncr; # this comment is needed for Windows system
# set -x
# Reading parameters
declare -a arg_val # values of dir scr out param arguments
declare -a req_args=('source_dir=' 'source_files=' 'out_dir=' 'out_file=' 'param=')


for arg in "$@"
do
 i=0
 for ra in "${req_args[@]}"
 do
#  echo $arg
#  echo ${arg}
#  echo ${arg}
 if [[ "$arg" =~ "$ra" ]]; then
    # echo "arg = " ${arg}
    # echo "ra = " ${ra}
    arg0="${arg#*=}"
    if [[ "$arg0" =~ "cygdrive" ]]
    then
      arg0="$(cygpath -w "$arg0")"
    fi
    arg_val[$i]="$arg0"
    # echo $i
    # echo ${arg_val[$i]}
  fi
  i=$(( i + 1))
 done
done

# cd "$root_path"
# pwd
# for arg in "${arg_val[@]}"
# do
#  echo -----------------
#  echo "$arg"
# done

# Preparation of paths
path_of_script()
{
    local path=$(dirname "$(readlink -f "$0")") # Path of the file
    # In case if we are using cygwin we convert path:
    if [[ "$path" =~ "cygdrive" ]]
    then
      path="$(cygpath -w "$path")"
    fi
    echo "$path"
}

root_path=$(path_of_script) # path of this script
cd "$root_path"
source compilation_functions.sh
libs=("$root_path" "crproplib" "pulsarlib" "physlib" "mfufm" "slatec only_L_path")
lib_dependencies=$( create_paths_for_libraries "${libs[@]}" )
sep=$( path_separator "$root_path" )

# Source and output files
scr="${arg_val[1]}"
out="${arg_val[2]}$sep${arg_val[3]}"

echo "HEYYY" "${arg_val[@]}"

mkdir -p "${arg_val[2]}"

# Compiler options
compiler="gfortran"
compiler_options="-W -fbackslash"

# Compilation
cd "${arg_val[0]}"
compilation=$(echo "$compiler" "$compiler_options" "$scr" "$lib_dependencies" -o "$out")
compilation=${compilation//"\\"/"\\\\"} # changing one backslashe to two backslashes
#echo "$compilation"
eval "$compilation"

# Execution
rm ./*.mod 2>/dev/null
if [[ "$out" =~ ".exe" ]]
    then
      cmd /C "$out"
    else
      "$out"
    fi
rm $out 2>/dev/null
