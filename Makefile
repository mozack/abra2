# Make file for ABRA
# libAbra is invoked from the ABRA java code

SRCDIR=/home/lmose/code/abra/src/main/c

all: clean java native

java:
	mvn package

native:
	g++ -g -I$(SRCDIR) -I$(JAVA_HOME)/include -shared $(SRCDIR)/assembler.c -o target/libAbra.so

standalone:
	g++ -g -I. -I$(JAVA_HOME)/include assembler.c -o abra

clean:
	rm -rf target
