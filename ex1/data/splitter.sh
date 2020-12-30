if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters"
	exit 1
fi

linenumber=`grep -n ">id:${1} " human_data.fasta | cut -f1 -d:`

linenumber="$(($linenumber-1))"

head -n $linenumber $2