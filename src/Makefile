CC=gcc
CFLAGS=-Wall -O2
BIN=../bin
all:$(BIN)/ttview  $(BIN)/bamsorted $(BIN)/bam2wig 
#$(BIN)/faidx.cgi


checkenv:
	echo "Compiling with SAMDIR=${SAMDIR} and TABIXDIR=${TABIXDIR}"

$(BIN)/ttview: ttview.c $(BIN) checkenv
	$(CC) -o $@ -DSTANDALONE_VERSION  -I${SAMDIR} -L${SAMDIR} -L${SAMDIR}/bcftools  $<  ${SAMDIR}/bam2bcf.o  ${SAMDIR}/errmod.o  ${SAMDIR}/bam_color.o ${SAMDIR}/libbam.a -lbcf  -lm -lz



$(BIN)/bamsorted:bamsorted.c $(BIN) checkenv
	$(CC) ${CFLAGS} -o $@ -I ${SAMDIR} -L ${SAMDIR} $< -lbam -lz
$(BIN)/bam2wig:bam2wig.c $(BIN) checkenv
	$(CC) ${CFLAGS} -o $@ -I ${SAMDIR} -L ${SAMDIR} $< -lbam -lz
$(BIN)/faidx.cgi:faidxcgi.c cgi.c $(BIN) checkenv
	$(CC) ${CFLAGS}  -o $@ -I ${SAMDIR} -L ${SAMDIR} -DBUILD=$(BUILD) -DGENOME_PATH=$(GENOME_PATH)  $< cgi.c ${SAMDIR}/faidx.o ${SAMDIR}/razf.o ${SAMDIR}/knetfile.o -lz




	
test:test-samtools $(BIN)/bam2wig
	$(BIN)/bam2wig ${SAMDIR}/examples/toy.bam
	$(BIN)/bam2wig ${SAMDIR}/examples/toy.bam "ref2:10-20"
	
test-samtools:
	${SAMDIR}/samtools view -b ${SAMDIR}/examples/toy.sam -t ${SAMDIR}/examples/toy.fa -o ${SAMDIR}/examples/toy.bam
	${SAMDIR}/samtools index  ${SAMDIR}/examples/toy.bam

clean:
	rm -f $(BIN)/battview $(BIN)/bamsorted $(BIN)/selectflag
