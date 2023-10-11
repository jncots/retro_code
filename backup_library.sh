#!/bin/bash
(set -o igncr) 2>/dev/null && set -o igncr; # this comment is needed for Windows system

backup_lib="backup_lib"
exclude="--exclude={*.o,*.mod}"
#exclude="--exclude={*.o,*.mod} --exclude='slatec'"

path=$(dirname "$(readlink -f "$0")")
if [[ "$path" =~ ":" ]]
then
    path="$(cygpath -u "$path")"
fi

lib_path="$path"
parent_dir="$(dirname "$lib_path")"
backup_dir="$parent_dir/""$backup_lib"
mkdir -p "$backup_dir"

timestamp=$(date "+%Y%m%d-%H%M%S")
tar_file="$backup_dir/library_$timestamp.tar.gz"
command=$(echo tar "$exclude" -czf $tar_file ./)
command=${command//"\\"/"\\\\"} # changing one backslashe to two backslashes
eval "$command"


