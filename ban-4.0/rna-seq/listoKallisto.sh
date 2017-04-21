#!/usr/bin/env bash

cd $PRWD
[ "$KALL_OUT" != "" ] || KALL_OUT="kall_out"
[ -d "$KALL_OUT" ] || mkdir $KALL_OUT
echo -n "Transcriptome fasta file ?: "; read INFILE; [ "$INFILE" != "" ] || (echo "File name is required"; exit 1)
INDEX_OP_COMMANDS=""
tjob=$(echo "module load kallisto/0.43.0
cd $PRWD
kallisto index -i $PRWD/$INDEX/$INFILE \
$INDEX_OP_COMMANDS \
$PRWD/$TRANSCRIPTOMES/$INFILE" | qsub -N "Kall_Index" -l nodes=1:ppn=8,mem=20gb,vmem=20gb -V -q default)

while IFS='' read -r LINE || [[ -n "$LINE" ]]; do
	[[ "$LINE" =~ \#.* ]] &&  continue
	arrLINE=(${LINE//"|"/ })
	PRONAME=${arrLINE[0]}
	FILE=$(echo ${arrLINE[1]} | sed 's/;/ /g')
	INDEX_OP_COMMANDS=$(echo ${arrLINE[2]} | sed 's/;/ /g')
	COUNTS_OP_COMMANDS=$(echo ${arrLINE[3]} | sed 's/;/ /g')
	echo "
module load kallisto/0.43.0
cd $PRWD/$READS
kallisto quant -i $PRWD/$INDEX/$INFILE -o $PRWD/$KALL_OUT/$PRONAME \
$COUNTS_OP_COMMANDS \
$FILE
cd $PRWD" | qsub -N $PRONAME -l nodes=1:ppn=8,vmem=20gb,mem=20gb -W depend=afterok:$tjob -V -q default
done < $PRWD/$COMMANDS/$KALLISTO_CMD


