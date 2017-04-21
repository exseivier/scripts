#!/usr/bin/env bash

TRANS=$1
OUTPUT=$2
echo "Counting L genes"
M=$(grep "\.L|" $1 | wc -l)
echo "Counting S genes"
N=$(grep "\.S|" $1 | wc -l)

echo "$M	$N"

echo "ID	M	N	K	X	Y" > $2
for file in *;
do
	echo "Analising $file..."
	X=$(grep "\.L" $file | wc -l)
	Y=$(grep "\.S" $file | wc -l)
	K=$(echo "$X + $Y" | bc -l)
	echo "$file	$M	$N	$K	$X	$Y" >> $2
done
echo "Success!"
