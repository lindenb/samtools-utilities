CC=gcc
CFLAGS=-Wall -O2

all:bin/ttview bin/selectflag bin/bamsorted bin/bam2wig


checkenv:
	echo "Compiling with SAMDIR=${SAMDIR}"

bin/ttview: src/ttview.c bin checkenv
	$(CC) -o $@ -DSTANDALONE_VERSION  -I${SAMDIR} -L${SAMDIR} -L${SAMDIR}/bcftools  $<  ${SAMDIR}/bam2bcf.o  ${SAMDIR}/errmod.o  ${SAMDIR}/bam_color.o ${SAMDIR}/libbam.a -lbcf  -lm -lz

bin/selectflag:script/selectflag.sh bin
	cp $< $@
	chmod +x $@

bin/bamsorted:src/bamsorted.c bin checkenv
	$(CC) ${CFLAGS} -o $@ -I ${SAMDIR} -L ${SAMDIR} $< -lbam -lz
bin/bam2wig:src/bam2wig.c bin checkenv
	$(CC) ${CFLAGS} -o $@ -I ${SAMDIR} -L ${SAMDIR} $< -lbam -lz


bin:
	mkdir -p bin

test:test-samtools bin/bam2wig
	bin/bam2wig ${SAMDIR}/examples/toy.bam
	bin/bam2wig ${SAMDIR}/examples/toy.bam "ref2:10-20"
	
test-samtools:
	${SAMDIR}/samtools view -b ${SAMDIR}/examples/toy.sam -t ${SAMDIR}/examples/toy.fa -o ${SAMDIR}/examples/toy.bam
	${SAMDIR}/samtools index  ${SAMDIR}/examples/toy.bam

clean:
	rm -f bin/battview bin/bamsorted bin/selectflag
