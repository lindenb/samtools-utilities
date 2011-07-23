CC=gcc
CFLAGS=-Wall -O2

all:bin/ttview bin/selectflag bin/bamsorted


checkenv:
	echo "Compiling with SAMDIR=${SAMDIR}"

bin/ttview: src/ttview.c bin checkenv
	$(CC) -o $@ -DSTANDALONE_VERSION  -I${SAMDIR} -L${SAMDIR} -L${SAMDIR}/bcftools  $<  ${SAMDIR}/bam2bcf.o  ${SAMDIR}/errmod.o  ${SAMDIR}/bam_color.o ${SAMDIR}/libbam.a -lbcf  -lm -lz

bin/selectflag:script/selectflag.sh bin
	cp $< $@
	chmod +x $@

bin/bamsorted:src/bamsorted.c bin checkenv
	$(CC) ${CFLAGS} -o $@ -I ${SAMDIR} -L ${SAMDIR} $< -lbam -lz

bin:
	mkdir -p bin

clean:
	rm -f bin/battview bin/bamsorted bin/selectflag
