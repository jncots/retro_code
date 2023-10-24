#!/bin/bash
(set -o igncr) 2>/dev/null && set -o igncr; # this comment is needed for Windows system

# This script contains functions for creating dependencies paths:
#
# create_library_paths_for_gcc() and create_paths_for_libraries()
#
# and compilations of libraries:
# library_compilation()
#


# path_separator() distinguishes between Windows backslashe \ 
# and Linux slash / separations in the path
# i.e. C:\something\something and /home/something/something
path_separator()
{
  if [[ "$1" == *[/]* ]]
  then
    local sep="/"
  else
    local sep="\\"
  fi
  echo "$sep"
}


create_library_paths_for_gcc()
{
  local root_path=$1
  local lib_name=$2
  local only_L_path=$3
  
  if [[ "$root_path" == *[/]* ]]
  then
    local sep="/"
  else
    local sep="\\"
  fi

  local lib_path="$root_path""$sep""$lib_name"
  local mod="-I""$lib_path""$sep""mod"
  local lib="-L""$lib_path"" -l""$lib_name"  
  
  if [ "$only_L_path" == "only_L_path" ]; then
      local result_path="$lib"
  else
      local result_path="$mod"" $lib" 
  fi
  echo "$result_path"
}

### Example of usage:
#root_path="C:\lib\something"
#lib_name="cool"
#result=$( create_library_paths_for_gcc "$root_path" "$lib_name" "only_L_path")
#echo $result

create_paths_for_libraries()
{
    local -a argv=("${@}")
    local -i totlen=${#argv[@]}
    local root_path="${argv[0]}"
    local -a lib_names=("${argv[@]:1:$totlen}")
    local all_paths=""
   
    for lib_name in "${lib_names[@]}"
    do
        if [[ "$lib_name" =~ " only_L_path" ]]; then # if there is substring "only_L_path"
          lib_name=${lib_name//" only_L_path"/}      # remove "only_L_path" substring
          lib_paths=$( create_library_paths_for_gcc "$root_path" "$lib_name" "only_L_path")
          all_paths="$all_paths"" $lib_paths"
        else
          lib_paths=$( create_library_paths_for_gcc "$root_path" "$lib_name")
        all_paths="$all_paths"" $lib_paths"
        fi
    done
    
    all_paths="${all_paths## }" # removing of the first space
    echo "$all_paths"
}


### Example of usage:
## source create_paths_for_libraries.sh # if it used in another script
#root_path="C:\work\library"
#declare -a libs=("$root_path" "mfufm" "slatec only_L_path" "physlib" "pulsarlib")
#lib_paths=$( create_paths_for_libraries "${libs[@]}" )
#echo $lib_paths

library_compilation()
{
    local -a argv=("${@}")
    local -i totlen=${#argv[@]}
    local root_path="${argv[0]}"
    local lib_name="${argv[1]}"
    local -a libs=("${argv[@]:2:$totlen}")
    if [[ $libs != "" ]]; then
      lib_dependencies=$( create_paths_for_libraries "$root_path" "${libs[@]}" )
    else
      lib_dependencies=""
    fi  


    local compiler="gfortran"
    local compiler_options="-c -W -fbackslash"
    # local compiler_options="-c -W -fbackslash -fallow-argument-mismatch" # for slatec


    local sep=$( path_separator "$root_path" )
    local lib_path="$root_path"$sep"$lib_name"
    local src="$lib_path"$sep"src"
    local obj="$lib_path"$sep"obj"
    local mod="$lib_path"$sep"mod"
    local modules_dir="-J"$mod

    mkdir -p "$src"
    mkdir -p "$obj"
    mkdir -p "$mod"

    local lib_archive="lib""$lib_name"".a"


    sep_line="___________________________________________________________\n"
    echo -e "$sep_line"
    echo root_path="$root_path"
    echo -e "$sep_line"
    echo lib_name="$lib_name"
    echo -e "$sep_line"
    echo lib_dependencies="$lib_dependencies"
    echo -e "$sep_line"


    # Cleaning
    cd "$obj"
    for file in *.o; do
        rm "$file"
    done

    # Compilation
    cd "$src"
    compilation=$(echo "$compiler" "$compiler_options" "$lib_dependencies" "$modules_dir" *.f90)
    # compilation=$(echo "$compiler" "$compiler_options" *.f) # for slatec
    compilation=${compilation//"\\"/"\\\\"} # changing one backslashe to two backslashes
    echo "$compilation"
    eval "$compilation"

    # # Moving new object files
    for file in *.o; do
        mv "$file" "$obj"
    done
     
    # Creation of the archive
    cd $obj
    ar cr "$lib_archive" *.o

    # Cleaning
    mv "$lib_archive" "$lib_path"
    cd "$root_path"
    
    echo -e "$sep_line"
    echo "$lib_name"" is compiled!"
    echo -e "$sep_line"
}


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

#lib_name="physlib"
#lib_dependencies=("mfufm" "slatec only_L_path" "physlib" "pulsarlib")
#library_compilation "$root_path" "$lib_name" "${lib_dependencies[@]}"