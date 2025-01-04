#!/bin/bash

WORK=`pwd`

rm -f ${WORK}/inherent.lammpstrj
touch ${WORK}/inherent.lammpstrj

for id in {1..3001}; do
	cd ${WORK}/DATA/${id}/
	cat inherent.lammpstrj >> ${WORK}/inherent.lammpstrj
done
