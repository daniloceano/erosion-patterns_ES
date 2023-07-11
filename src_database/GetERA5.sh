#!/bin/bash -l

# Adapted from https://dreambooker.site/2019/10/03/Initializing-the-WRF-model-with-ERA5-pressure-level/

CODEDIR=../met_data/scripts
DATADIR=../met_data/data

mkdir ../met_data/
mkdir ../met_data/ERA5/
mkdir $CODEDIR
mkdir $DATADIR

# Domain limts 
Nort=0
West=-40
Sout=-60
East=-10

while IFS=, read -r name start end f; do

	NAME=$name
	
	YY1=`echo $start | cut -c1-4`
	MM1=`echo $start | cut -c6-7`
	DD1=`echo $start | cut -c9-10`
	
	YY2=`echo $end | cut -c1-4`
	MM2=`echo $end | cut -c6-7`
	DD2=`echo $end | cut -c9-10`

	DATE1=$YY1$MM1$DD1
	DATE2=$YY2$MM2$DD2
	
	echo 'Downloading data for system:' $NAME' dates from '$DATE1' to '$DATE2
	echo 'North '$Nort
	echo 'West '$West
	echo 'South '$Sout
	echo 'East '$East 

	# generate and run python script
	sed -e "s/DATE1/${DATE1}/g;s/NAME/${NAME}/g;s/DATE2/${DATE2}/g;s/Nort/${Nort}/g;s/West/${West}/g;s/Sout/${Sout}/g;s/East/${East}/g;" GetERA5-pl.py > GetERA5-${NAME}-pl.py
	python GetERA5-${NAME}-pl.py

	# move donwloaded data
	mkdir -p ${DATADIR}/
	mv ERA5-${NAME}-pl.nc ${DATADIR}/

	# move the generated files
	mv GetERA5-${NAME}-pl.py ${CODEDIR}/

done < ../database/dates-limits.csv

exit 0
