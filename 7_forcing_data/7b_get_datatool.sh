#!/bin/bash
# This script is used to get the datatool from the repository
dt_url="https://github.com/kasra-keshavarz/datatool/"
dt_des="/globalhome/wmk934/HPC/datatool/"
rm -rf $dt_des
git clone $dt_url $dt_des