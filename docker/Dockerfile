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
       make && \
       ln -s /abra2/target/abra2-${ABRA2_VERSION}-jar-with-dependencies.jar /abra2.jar


ENTRYPOINT [ "/usr/bin/java", "-Xmx16G","-jar","/abra2.jar"]

