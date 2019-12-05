#!/usr/bin/env zsh
cases_arr=(timoshenkobeam sphericaljoint hingejoint fixedjoint pullingmuscle snake flagella muscularsnake walker elbow wing)
make clean
for my_case in ${cases_arr[@]}; do make -j8 "$my_case" CC="g++-8" && make clean; done
