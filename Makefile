# Make file for ABRA
# libAbra is invoked from the ABRA java code

SRCDIR=src/main/c

all: clean java native

java:
	mvn clean package

native:
	g++ -g -I$(SRCDIR) -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -shared -fPIC $(SRCDIR)/assembler.c -o target/libAbra.so

standalone:
	g++ -g -I. -I$(JAVA_HOME)/include assembler.c -o abra

clean:
	mvn clean
	rm -rf target
