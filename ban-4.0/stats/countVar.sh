#!/usr/bin/env bash

echo ""
echo "Processing $1..."
echo ""
echo "Total genes: " | tr -d "\n"
grep "gene" $1 | wc -l

echo "Total L genes: " | tr -d "\n"
grep "gene" $1 | cut -f9 | cut -d ";" -f2 | grep "\.L" | wc -l

echo "Total S genes: " | tr -d "\n"
grep "gene" $1 | cut -f9 | cut -d ";" -f2 | grep "\.S" | wc -l

echo "Total orphan genes: " | tr -d "\n"
grep "gene" $1 | cut -f9 | cut -d ";" -f2 | grep -v "\.L" | grep -v "\.S" | wc -l

echo "Total genes with +1 transcript: " | tr -d "\n"
grep "mRNA" $1 | cut -f9 | cut -d ";" -f2 | sort | uniq -c | grep -v " 1" | wc -l

echo "L genes with +1 transcript: " | tr -d "\n"
grep "mRNA" $1 | cut -f9 | cut -d ";" -f2 | sort | uniq -c | grep -v " 1" | grep "\.L" | wc -l

echo "S genes with +1 transcript: " | tr -d "\n"
grep "mRNA" $1 | cut -f9 | cut -d ";" -f2 | sort | uniq -c | grep -v " 1" | grep "\.S" | wc -l

echo "Orphan genes with +1 transcript: " | tr -d "\n"
grep "mRNA" $1 | cut -f9 | cut -d ";" -f2 | sort | uniq -c | grep -v " 1" | grep -v "\.L" | grep -v "\.S" | wc -l

echo ""
echo "Processing $2..."
echo ""

echo "total genes with +1 transcript for alternative file: " | tr -d "\n"
grep "gene" $2 | wc -l

echo "L genes with +1 transcript for alternative file: " | tr -d "\n"
grep "gene" $2 | grep "\.L" | wc -l

echo "S genes with +1 transcript for alternative file: " | tr -d "\n"
grep "gene" $2 | grep "\.S" | wc -l

echo "Orphan genes with +1 transcript for alternative file: " | tr -d "\n"
grep "gene" $2 | grep -v "\.L" | grep -v "\.S" | wc -l

