all:bin bin/selectflag
	(cd src; make )
bin/selectflag:script/selectflag.sh bin
	cp $< $@
	chmod +x $@
bin:
	mkdir $@
lib:
	mkdir $@
clean:
	(cd src; make clean)
	rm -f script/selectflag.sh
