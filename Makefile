CC=gcc
CFLAGS=-O3
LDFLAGS=-static -lm
CFILES= applythresh.c chainelim.c distributionanalyzer.c fastexp2.c gettimer.c mainexactpscore.c mainmisalign.c mainscreening.c mainshevek.c misalignnrpang.c numrecpang.c numrecutilities.c prcerr.c randnumgen.c rcount2pang.c redundanteliminator.c support.c
HEADERFILES=definitions.h fastexp2.h fastexp2.inc nrutilp.h rcont2p.h shevek.h support.h

shevek: ${CFILES} ${HEADERFILES}
	${CC} ${CFLAGS} ${CFILES} -o $@ ${LDFLAGS}

clean:
	rm shevek
