

utr2fasta_2.py:

load_genome(): GENOME data structure
HASH[STR:STR]
genome[Sequence name] +------> [Sequence]

main(): gene_models data structure
HASH[STR:(INT|STR|ARRAY[HASH[STR:STR]])
gene_models[geneID]	+----------> [Gene name]-STR
					|
					+----------> [Start]-INT
					|
					+----------> [End]-INT
					|
					+----------> [Orientation]-STR
					|
					+----------> [Locus name]-STR
					|
					+----------> [Genetic element]-ARRAY	+ ----> [Start]-INT
															|
															+-----> [End]-INT
					|----------------------------------|	|
								First layer					+-----> [ID]-STR
															|
															+-----> [Parent]-STR

															|-------------------|
																Second layer



