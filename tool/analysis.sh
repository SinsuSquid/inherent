#!/bin/bash

ROOT="/mnt/CALC_SSD/bgkang/inherent/tool"
CWD=`pwd`

for ensemble in "ensemble00" "ensemble01" "ensemble02" "ensemble03"; do
	cd ${CWD}/${ensemble}
	# ${ROOT}/fskt ./true.lammpstrj /dev/null /dev/null ./fskt_true.out /dev/null
	# ${ROOT}/fskt ./inherent.lammpstrj /dev/null /dev/null ./fskt_inherent.out /dev/null
	# python ${ROOT}/energy.py ./ 3001 ./energy_profile.out
	${ROOT}/config_dist ./inherent.lammpstrj ./config_dist.out
done
