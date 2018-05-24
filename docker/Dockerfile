FROM  alpine
LABEL MAINTAINER "Lisle Mose"

RUN apk -U add \
     g++ \ 
     git \
     make \
     maven \
     openjdk8

ARG ABRA2_VERSION=2.15
ENV ABRA2_VERSION ${ABRA2_VERSION}

RUN git clone https://github.com/mozack/abra2.git && \
       cd abra2 && \
       git checkout tags/v${ABRA2_VERSION} && \ 
       make

ENV JAVA_MEM_REQ 5G

ENTRYPOINT ["sh", "-c", "echo For this Dockerized version of ABRA2-v${ABRA2_VERSION}, add -e JAVA_MEM_REQ=$JAVA_MEM_REQ to \n echo the 'docker run ...' to change the amount of memory for Java to request ; /usr/bin/java -Xmx${JAVA_MEM_REQ} -jar /abra2/target/abra2-${ABRA2_VERSION}-jar-with-dependencies.jar "]

