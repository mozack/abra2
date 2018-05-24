FROM  alpine
LABEL MAINTAINER "Lisle Mose"

RUN apk -U add \
     g++ \ 
     git \
     make \
     maven \
     openjdk8

RUN git clone https://github.com/mozack/abra2.git && \
       cd abra2 && \
       make
       
ENTRYPOINT ["/usr/bin/java", "-Xmx5G", "-jar", "/abra2/target/abra2-2.15-jar-with-dependencies.jar"]
