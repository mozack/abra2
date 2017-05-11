# Make file for ABRA
# libAbra is invoked from the ABRA java code

SRCDIR=src/main/c

all: clean native ssw java

java:
	mvn package

mktargetdir:
	mkdir target

native: mktargetdir
	g++ -g -O2 -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -shared -fPIC $(SRCDIR)/assembler.cpp $(SRCDIR)/sg_aligner.cpp -o target/libAbra.so

clean_ssw:
	cd src/main/c/ssw && make clean

ssw:
	cd src/main/c/ssw && make libssw.so && make libsswjni.so && cp libssw*.so ../../../../target

standalone:
	g++ -g -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux $(SRCDIR)/assembler.c -o abra

sga:
	g++ -g -O2 -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux $(SRCDIR)/sg_aligner.cpp -o sga

clean:	clean_ssw
	rm -rf target
	mvn clean

# TODO: Parameterize version
javah: java
	javah -classpath target/abra-0.53-SNAPSHOT.jar abra.NativeAssembler
