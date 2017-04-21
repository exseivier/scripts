#!/usr/bin/env bash
<<MESSAGE
	AUTHOR: Javier Montalvo-Arredondo
	CONTACT: javier.montalvo [at] cinvestav.mx
	VERSION: 0

USAGE:
	Just write the program name and press ENTER, jejeje!
	This script has an interactive data-input interface

MESSAGE


clear
echo "WELCOME!"
echo "Running Sylamer!"
echo ""
DATE=`date +%Y-%m-%d:%H:%M:%S`

# This is the Current Working Directory where inputs are stored, and
# where results will be placed.
echo -n "Set project working directory: "; read PROJECT_WD
if [ "$PROJECT_WD" == "" ]; then
	PROJECT_WD=`pwd`
	cd $PROJECT_WD
	ls -hl
	echo "You are in $PROJECT_WD"
else
	cd $PROJECT_WD
	ls -hl .
	echo "You are in $PROJECT_WD"
fi

# This is the name of the results directory.
# Creates this directory
WORKING_DIR="sylResults"
[ -d "$WORKING_DIR" ] || mkdir $WORKING_DIR

#ls -hl $PRWD/$TRANSCRIPTOMES

# This name will be used as specific directory name.
# Every result from sylamer of every analysed contrast
# will be allocated to $PROJECT_WD/$WORKING_DIR/$PROJECT_NAME/
echo -n "Project NAME ?: "; read PROJECT_NAME
if [ "$PROJECT_NAME" == "" ]; then
	PROJECT_NAME="syl"
fi

# This suffix name identifies every resulting file from sylamer
# It should emphasis the contrast used in the statistic analysis.
echo ""
echo -n "Output file suffix ?: "; read SUFFIX

# This is the name of the query sequences (UTRs/mRNA/Promoters/etc...)
# file used in the sylamer analysis
echo -n "FASTA ?: "; read FASTA_FILE
if [ "$FASTA_FILE" == "" ]; then
	echo "Fasta file of UTRs is needed!"
	exit 1
fi

# This variable stores the name of the header file.
# This file contais the headers of the query sequences ordered
# ascendently or descendently, depending on the expression level 
echo ""
echo -n "HD ?: "; read HEADER_FILE
if [ "$HEADER_FILE" == "" ]; then
	echo "The ordered headers file of UTRs is needed!"
	exit 1
fi

# If you want to restrict the sylamer searching, you could specify
# a words file that contains the specific words you want to search for.
# The format is a tab delimited file: the words are placed in the first column,
# and comments are placed in the second column
echo ""
echo -n "WORDS ?: "; read WORDS_FILE
if [ "$WORDS_FILE" == "" ]; then
	WORDS_COMM=""
	echo "WARNING!"
	echo "No words file was selected!"
else
	WORDS_COMM="-words ${WORDS_FILE}"
fi

# This is the word length that sylamer uses in the analysis
echo ""
echo -n "Word length ?: "; read WORD_LEN
if [ "$WORD_LEN" == "" ]; then
	WORD_LEN=6
fi

# This is the amount of sequences sylamer will increase each step
echo ""
echo -n "Growing window ?: "; read GROW
if [ "$GROW" == "" ]; then
	GROW=100
fi

# In this variable you could store additional parameters.
# Those parameters must fit the format required by sylamer
echo ""
echo -n "Optional commands ?: "; read OPTIONAL

# Specify the main title of the graph
echo -n "Title ?: "; read TITLE
# Specify the x axis label
echo -n "X axis label ?: "; read XLAB

# Sending the task to cluster
clear
echo ""
echo "Task information"
echo "Project name: ${PROJECT_NAME}_${DATE}"
[ -d "$WORKING_DIR/$PROJECT_NAME" ] || mkdir $WORKING_DIR/$PROJECT_NAME
CMD="cd $PRWD/$TRANSCRIPTOMES; sylamer -fasta $FASTA_FILE -o $PROJECT_WD/$WORKING_DIR/$PROJECT_NAME/${FASTA_FILE%.*}_${SUFFIX}.syl -universe $PROJECT_WD/$HEADER_FILE -k $WORD_LEN $WORDS_COMM -grow $GROW $OPTIONAL
echo '#'$TITLE > tmp.txt
echo '#'$XLAB >> tmp.txt
cat $PROJECT_WD/$WORKING_DIR/$PROJECT_NAME/${FASTA_FILE%.*}_${SUFFIX}.syl >> tmp.txt
mv tmp.txt $PROJECT_WD/$WORKING_DIR/$PROJECT_NAME/${FASTA_FILE%.*}_${SUFFIX}.syl"
echo "$CMD"
taskID=$(echo "$CMD" | qsub -N ${PROJECT_NAME}_${DATE} -l nodes=1:ppn=8,mem=20gb,vmem=20gb -V -q default)
echo ""
echo "Task was sent to queue, task ID is: ${taskID}"
echo "runSyl has finished!"
