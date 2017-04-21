#!/usr/bin/env bash
. ~/scripts/ban-4.0/bash_lib/jlib.sh
TOTAL_L=$(shapefasta $1 | grep ">" | grep "\.L" | wc -l)
TOTAL_S=$(shapefasta $1 | grep ">" | grep "\.S" | wc -l)
[ -f $3 ] || touch $3
echo "miRNA_name	Percent_on_L	Percent_on_S	Seed"
echo "miRNA_name	Percent_on_L	Percent_on_S	Seed" >> $3
while read -r LINE;
do
	arrLINE=(${LINE/"\t"/ / })
	NAME=${arrLINE[1]}
	SEED=${arrLINE[0]}
	TARGET_L=$(shapefasta $1 | grep -B 1 "$SEED" | grep ">" | grep "\.L" | wc -l)
	TARGET_S=$(shapefasta $1 | grep -B 1 "$SEED" | grep ">" | grep "\.S" | wc -l)
	PERCENT_TARGET_L=$(echo $TARGET_L*100/$TOTAL_L | bc -l)
	PERCENT_TARGET_S=$(echo $TARGET_S*100/$TOTAL_S | bc -l)
	echo "$NAME	$PERCENT_TARGET_L	$PERCENT_TARGET_S	$SEED"
	echo "$NAME	$PERCENT_TARGET_L	$PERCENT_TARGET_S	$SEED" >> $3
done < $2
