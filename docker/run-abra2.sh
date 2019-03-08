#! /bin/sh

if [ -z $JAVA_OPTS ] ; then
    JAVA_OPTS="-Xmx16g"
fi

if [ -z $1 ] ; then 

    java $JAVA_OPTS -jar /abra2.jar 

    echo ""
    echo "To use this Docker image, structure your command as follows:"
    echo ""
    echo "   docker run abra:$ABRA2_VERSION java $JAVA_OPTS -jar /abra2.jar [options for ABRA2...]"
    echo ""
    echo "or pass in your JAVA_OPTS this way"
    echo ""
    echo "   docker run --env JAVA_OPTS abra:$ABRA2_VERSION sh /run-abra2.sh [options for ABRA2...]"
    echo ""
    echo "Pay special care to set the --tmpdir option to a good location"

else

    java $JAVA_OPTS -jar /abra2.jar $@

fi
