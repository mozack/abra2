# Make file for ABRA
# libAbra is invoked from the ABRA java code

all:
	g++ -g -I. -I$(JAVA_HOME)/include -shared assembler.c -o libAbra.so
#	g++ -g -I. -I/usr/lib/jvm/java-6-openjdk/include -shared assembler.c -o libAbra.so

standalone:
	g++ -g -I. -I$(JAVA_HOME)/include assembler.c -o abra

clean:
	rm libAbra.so
