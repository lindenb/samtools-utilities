#!/bin/sh
dialog --stdout --title "SAM FLAGS" \
	--nocancel \
	--separate-output \
	--checklist "SELECT FLAGS" 20 60 11 \
	1	"Read Paired" off \
	2	"Read mapped in proper pair" off \
	4	"Read unmapped" off \
	8	"Mate unmapped" off \
	16	"Read reverse strand" off \
	32	"Mate reverse strand" off \
	64	"First in pair" off \
	128	"Second in pair" off \
	256	"Not primary alignment" off \
	512	"Read fails platform/vendor quality checks" off \
	1024	"Read is PCR or optical duplicate" off \
| awk 'BEGIN { i=0;} {i+=int($1);} END {printf("%s\n",i);}'
