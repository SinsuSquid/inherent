#!/bin/bash

ROOT="/mnt/CALC_SSD/bgkang/inherent/tool"
CWD=`pwd`

for temp in "300K" "400K" "500K"; do
	for ensemble in "ensemble00" "ensemble01" "ensemble02" "ensemble03"; do
		cd ${CWD}/${temp}/${ensemble}
		# ${ROOT}/fskt ./true.lammpstrj /dev/null /dev/null ./fskt_true.out /dev/null
		# ${ROOT}/fskt ./inherent.lammpstrj /dev/null /dev/null ./fskt_inherent.out /dev/null
		python ${ROOT}/energy.py ./ 3001 ./energy_profile.out
	done
done
