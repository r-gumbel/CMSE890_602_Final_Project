# target output files for the whole workflow
rule all:
	input:
		"/mnt/home/gumbelri/Snakemake/ca48pb208_101V.dat"

rule dataMcReadyRead:
	input:
		R1 = 
