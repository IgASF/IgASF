#! /bin/bash

methods=(global macroS macroN macroR element)

type=TYPE

if [ "$#" -lt 1 ]; then
    echo "$0 test_file[s]"
    echo ""
    echo "   test_file[s]  The test file[s] to be executed."
    echo ""
    echo "The timings are written to log.txt"
    echo "The errors between the matrices are written to diffs.txt"
    exit 1
fi

for test in "$@"; do
    for i in ${methods[*]}; do
        echo "./runTest-$type $1 -m $i -l log.txt -o tmp-$i.dat"
        ./runTest-$type $1 -m $i -l log.txt -o tmp-$i.dat
    done
    echo >> diffs.txt
    echo $test >> diffs.txt
    echo >> diffs.txt
    for i in ${methods[*]}; do
        DIFF=
        for j in ${methods[*]}; do
            DIFF0=`./compareMatrices-$type tmp-$i.dat tmp-$j.dat`
            DIFF="$DIFF $DIFF0"
        done
        echo $DIFF >> diffs.txt
    done
    for i in ${methods[*]}; do
        rm tmp-$i.dat
    done
    echo "Differences of assembled matrices are written to diffs.txt."
    echo
done

