#!/usr/bin/env bash

usage(){
	
	echo -e "\t\t*******************************"
	echo -e "\t\t****      Target-Scan      ****"
	echo -e "\t\t****         PIPE          ****"
	echo -e "\t\t****       TUTORIAL        ****"
	echo -e "\t\t*******************************"
	echo -e "\n\n"
	echo -e "\tUsage:"
	echo -e "\t\tts-tools -i < utr.fasta > -m < mirnaome.fasta > --type1 'L' --type2 'S' --window 100 --percent 25\n"
	echo -e "\tOptions:\n"
	echo -e "\t-i\t\tUTR sequences in fasta format."
	echo -e "\t-m\t\tmature miRNA sequences in fasta format."
	echo -e "\t--type1\t\tParalog gene 1 label."
	echo -e "\t--type2\t\tParalog gene 2 label."
	echo -e "\t--window\tWindow size in hypergeometric test of target-scan results."
	echo -e "\t--percent\tPercentage of top gene targets used for target-scan results analysis.\n"
	echo -e "\tWARNING:\n"
	echo -e "\t\tExport a environmental variable \"\$TSCS_HOME\" with the path where parameters' files folder is before executing ts-tools.\n"
}


while getopts "hi:m:-:" opt
do
	if [ -z $opt ]; then
		echo "arguments were not specified"
		clear
		usage
		exit 0
	fi
	case $opt in
		h)
			clear
			usage
			exit 0
		;;
		i)
			INPUT=$OPTARG
		;;
		m)
			MIRNA=$OPTARG
		;;
		-)
			case $OPTARG in
				type1)
					TYPE1=${!OPTIND}
					OPTIND=$(( $OPTIND + 1 ))
				;;
				type2)
					TYPE2=${!OPTIND}
					OPTIND=$(( $OPTIND + 1 ))
				;;
				window)
					WINDOW=${!OPTIND}
					OPTIND=$(( $OPTIND + 1 ))
				;;
				percent)
					PERCENT=${!OPTIND}
					OPTIND=$(( $OPTIND + 1 ))
				;;
				help)
					clear
					usage
					exit 0
				;;
			esac
		;;
		*)
			clear
			usage
			exit 1
		;;
	esac
done

INPUTFILE=$(basename $INPUT '/')
MIRNAFILE=$(basename $MIRNA '/')

echo "Preprossesing input files..."
#echo "
fa2ints $INPUTFILE
miRNAfamily-build $MIRNAFILE

echo "Finding miRNA targets..."
#echo "
targetScan ${MIRNAFILE%.*}.ints ${INPUTFILE%.*}.ints ${INPUTFILE%.*}.outs

echo "Emulating file..."
#echo "
target2context ${INPUTFILE%.*}.outs > ${INPUTFILE%.*}.outsc++

echo "Linking parameters..."
#echo "
ln -s ${TSCS_HOME}/TargetScan7_context_scores/Agarwal_2015_parameters.txt Agarwal_2015_parameters.txt
ln -s ${TSCS_HOME}/TargetScan7_context_scores/All_cell_lines.AIRs.txt All_cell_lines.AIRs.txt
ln -s ${TSCS_HOME}/TargetScan7_context_scores/TA_SPS_by_seed_region.txt TA_SPS_by_seed_region.txt

echo "Calculating context score..."
#echo "
targetScan-ctx++ ${MIRNAFILE%.*}.intsc++ ${INPUTFILE%.*}.ints ${INPUTFILE%.*}.outsc++ \
				 ${TSCS_HOME}/TargetScan7_context_scores/ORF_8mer_counts_sample.txt \
				 ${TSCS_HOME}/TargetScan7_context_scores/ORF_Sequences_sample.lengths.txt \
				 ${INPUTFILE%.*}.ts

echo "Separating results..."
#echo "
targetScan-final ${INPUTFILE%.*}.ts

echo "Preprocessing data..."
#echo "
mkdir ${INPUTFILE%.*}_out
mv *.data ${INPUTFILE%.*}_out
cd ${INPUTFILE%.*}_out
targetScan-pp "$TYPE1" "$TYPE2"
cd ..

echo "Sattistic analysis & ploting..."
#echo "
targetScan-stats ${INPUTFILE%.*}_out $WINDOW $PERCENT

echo "Unliking parameters..."
#echo "
unlink Agarwal_2015_parameters.txt
unlink All_cell_lines.AIRs.txt
unlink TA_SPS_by_seed_region.txt

echo "Everything is OK... ;)"

