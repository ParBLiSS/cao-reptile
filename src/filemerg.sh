# 
#   Author: Ankit Shah <shah.29ankit@gmail.com>
#           Sriram P C <sriram.pc@gmail.com>

#   Copyright 2011-2012 Ankit Shah, Sriram P C, Srinivas Aluru

#   This file is part of Parallel Reptile (version 1.1)

#   PReptile is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   PReptile is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU Lesser General Public License for more details.

#   You should have received a copy of the GNU Lesser General Public License
#   along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.



#!/bin/bash
processors=$1
inputfileprefix=$2
outputfile=$3
#echo $processors
#echo $inputfileprefix
#echo $outputfile
rm $outfile
for i in $(seq 0  $(($processors-1)))
do
    cat $inputfileprefix$i >> $outputfile
done

