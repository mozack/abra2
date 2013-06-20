# Make file for ABRA
# libAbra is invoked from the ABRA java code

SRCDIR=src/main/c

all: clean native java

java:
	mvn package

mktargetdir:
	mkdir target

native: mktargetdir
	g++ -g -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -shared -fPIC $(SRCDIR)/assembler.c -o target/libAbra.so

standalone:
	g++ -g -I. -I$(JAVA_HOME)/include assembler.c -o abra

clean:
	rm -rf target
	mvn clean
