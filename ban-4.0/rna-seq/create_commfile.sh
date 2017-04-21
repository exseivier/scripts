#!/usr/bin/env bash

key=$1

function create_kallisto {
	OC_index=";"
	OC_quant=";"
	[[ -e "$PRWD/$COMMANDS/$KALLISTO_CMD" ]] && > $PRWD/$COMMANDS/$KALLISTO_CMD
	while IFS="" read -r LINE || [[ -n "$LINE" ]];
	do
		files=$(ls $PRWD/$READS/trim${LINE%.sra}*_P.fastq.gz)
		afiles=(${files// / })
		files=""
		for item in ${afiles[@]};
		do
			parts=(${item//"/"/ })
			item=${parts[${#parts[@]}-1]}
			files=${files}${item}";"
		done
		echo ${LINE%.sra}"|"${files}"|"${OC_index}"|"${OC_quant} >> $PRWD/$COMMANDS/$KALLISTO_CMD
	done < $PRWD/$READS/$SRAFILES
	}
function create_trimmomatic {
	COMMS="HEADCROP:12;LEADING:28;TRAILING:28;SLIDINGWINDOW:4:28;MINLEN:36"
	[[ -e "$PRWD/$COMMANDS/$COMMAND_FILE" ]] && > $PRWD/$COMMANDS/$COMMAND_FILE
	while IFS="" read -r LINE || [[ -n "$LINE" ]];
	do
		files=$(ls $PRWD/$READS/${LINE%.sra}*.fastq.gz)
		afiles=(${files// / })
		files=""
		for item in ${afiles[@]};
		do
			parts=(${item//"/"/ })
			item=${parts[${#parts[@]}-1]}
			files=${files}${item}"|"
		done
		echo ${LINE%.sra}.sam"|"${files}${COMMS} >> $PRWD/$COMMANDS/$COMMAND_FILE
	done < $PRWD/$READS/$SRAFILES
	}


case $key in
	--trim)
	create_trimmomatic
	;;
	--kallisto)
	create_kallisto
	;;
	*)
	echo "Unknown parameter $1"
	;;
esac
shift

#__EOF__//~~CAT
