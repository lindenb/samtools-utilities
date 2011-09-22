CC=gcc
CFLAGS=-Wall -O2
BIN=../bin
all:$(BIN)/ttview $(BIN)/bamsorted $(BIN)/bam2wig 

checksamenv:
	if [ -z "${SAMDIR}" ]; then echo '###ERROR: the environment variable $${SAMDIR} is not defined.'; exit -1; fi
	echo "Compiling with SAMDIR=${SAMDIR}"
checktabixenv:
	if [ -z "${TABIXDIR}" ]; then echo '###ERROR: the environment variable $${TABIXDIR} is not defined.'; exit -1; fi
	echo "Compiling with TABIXDIR=${TABIXDIR}"

$(BIN)/ttview: ttview.c $(BIN) checksamenv
	$(CC) -o $@ -DSTANDALONE_VERSION  -I${SAMDIR} -L${SAMDIR} -L${SAMDIR}/bcftools  $<  ${SAMDIR}/bam2bcf.o  ${SAMDIR}/errmod.o  ${SAMDIR}/bam_color.o ${SAMDIR}/libbam.a -lbcf  -lm -lz
$(BIN)/jointabix:jointabix.c $(BIN) checktabixenv
	$(CC) -o $@ ${CFLAGS} -I${TABIXDIR} -L${TABIXDIR} $<  -ltabix -lz


$(BIN)/bamsorted:bamsorted.c $(BIN) checksamenv
	$(CC) ${CFLAGS} -o $@ -I ${SAMDIR} -L ${SAMDIR} $< -lbam -lz
$(BIN)/bam2wig:bam2wig.c $(BIN) checksamenv
	$(CC) ${CFLAGS} -o $@ -I ${SAMDIR} -L ${SAMDIR} $< -lbam -lz
$(BIN)/faidx.cgi:faidxcgi.c cgi.c $(BIN) checksamenv
	$(CC) ${CFLAGS}  -o $@ -I ${SAMDIR} -L ${SAMDIR} -DBUILD=$(BUILD) -DGENOME_PATH=$(GENOME_PATH)  $< cgi.c ${SAMDIR}/faidx.o ${SAMDIR}/razf.o ${SAMDIR}/knetfile.o -lz




	
test:test-samtools $(BIN)/bam2wig
	$(BIN)/bam2wig ${SAMDIR}/examples/toy.bam
	$(BIN)/bam2wig ${SAMDIR}/examples/toy.bam "ref2:10-20"
	
test-samtools:
	${SAMDIR}/samtools view -b ${SAMDIR}/examples/toy.sam -t ${SAMDIR}/examples/toy.fa -o ${SAMDIR}/examples/toy.bam
	${SAMDIR}/samtools index  ${SAMDIR}/examples/toy.bam

test-tabix:$(BIN)/jointabix
	wget -O cytoBand.txt.gz "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz"
	gunzip -q cytoBand.txt.gz
	${TABIXDIR}/bgzip cytoBand.txt
	${TABIXDIR}/tabix -f -s 1 -b 2 -e 3 -0 cytoBand.txt.gz
	echo "chr2	16500900	18600000" | $(BIN)/jointabix -c 1 -s 2 -e 3 -f cytoBand.txt.gz
	
clean:
	rm -f $(BIN)/battview $(BIN)/bamsorted $(BIN)/selectflag $(BIN)/bam2wig
