#!/bin/bash
(set -o igncr) 2>/dev/null && set -o igncr; # this comment is needed for Windows system

#
# This script uses library_compilation from compilation_functions.sh script
# to compile below descibed libraries using their dependencies.
#
source compilation_functions.sh

root_path=$(path_of_script) # path of this script

# lib_name="slatec"
# lib_dependencies=''
# library_compilation "$root_path" "$lib_name" "${lib_dependencies[@]}"

lib_name="mfufm"
lib_dependencies=''
library_compilation "$root_path" "$lib_name" "${lib_dependencies[@]}"

lib_name="physlib"
lib_dependencies=("mfufm")
library_compilation "$root_path" "$lib_name" "${lib_dependencies[@]}"

lib_name="pulsarlib"
lib_dependencies=("mfufm")
library_compilation "$root_path" "$lib_name" "${lib_dependencies[@]}"

lib_name="crproplib"
lib_dependencies=("mfufm" "pulsarlib")
library_compilation "$root_path" "$lib_name" "${lib_dependencies[@]}"

# $root_path"/backup_library.sh"