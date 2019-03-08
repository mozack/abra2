FROM  alpine
LABEL MAINTAINER "Lisle Mose <lmose@email.unc.edu>"
LABEL MAINTAINER "Alan Hoyle <alanh@unc.edu>"

RUN apk -U add \
     libc6-compat \ 
     openjdk8

ARG ABRA2_VERSION=2.20
ENV ABRA2_VERSION ${ABRA2_VERSION}
ENV JAVA_OPTS "-Xmx16G"

ADD https://github.com/mozack/abra2/releases/download/v${ABRA2_VERSION}/abra2-${ABRA2_VERSION}.jar /
ADD run-abra2.sh /

RUN chmod 755 /abra2-${ABRA2_VERSION}.jar /run-abra2.sh && \
    ln -s /abra2-${ABRA2_VERSION}.jar /abra2.jar 


# CMD java $JAVA_OPTS -jar /abra2.jar
CMD /run-abra2.sh
# CMD [ --help ]
