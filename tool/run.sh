#!/bin/bash

ROOT="/mnt/CALC_SSD/bgkang/inherent/test/inherent/"
WORK=`pwd`

for item in ${WORK}/parsing/*.lammps_data; do
	name=`basename ${item}`
	id=`echo ${name} | tr -d -c 0-9`
	rm -rdf ${WORK}/DATA/${id}
	mkdir ${WORK}/DATA/${id}
	cp ${ROOT}/BASE/* ${WORK}/DATA/${id}
	cp ${item} ${WORK}/DATA/${id}/true.lammps_data
	cd ${WORK}/DATA/${id}/
	sbatch slurm_inherent.sh &> /dev/null
done
