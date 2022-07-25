#!/bin/bash
set -e

##
## HiC-Pro
## Internal function
## Merge valid interactions files and remove duplicates
##

rmDup=0
prefix=""
while getopts ":dp:" opt; do
    case "$opt" in
        d) rmDup=1 ;;
        p) prefix=$OPTARG ;;
    esac
done
shift $(( OPTIND - 1 ))

vpairs="$@"
vpairs_sorted=$(echo $vpairs | sed -e 's/validPairs/sorted.validPairs/g')

mkdir -p ./tmp/

if [[ ${rmDup} == 1 ]]; then
    ## Sort individual validPairs files
    fcounts=0
    for vfile in ${vpairs}
    do
        echo "Sorting ${vfile} ..."
        fcounts=$((fcounts+1))
        ofile=$(echo ${vfile} | sed -e 's/validPairs/sorted.validPairs/')
        #sort -k2,2V -k3,3n -k5,5V -k6,6n -T ./tmp/ -o ${ofile} ${vfile}
        sort -k2,2 -k5,5 -k3,3n -k6,6n -T ./tmp/ -o ${ofile} ${vfile}
    done

    if [[ $fcounts -gt 1 ]]
    then
        echo "Merging and removing the duplicates ..."
        ## Sort valid pairs and remove read pairs with same starts (i.e duplicated read pairs)
        #sort -k2,2V -k3,3n -k5,5V -k6,6n -T ./tmp/ -m ${vpairs_sorted} | \
        sort -k2,2 -k5,5 -k3,3n -k6,6n -T ./tmp/ -m ${vpairs_sorted} | \
            awk -F"\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=$2 || c2!=$5 || s1!=$3 || s2!=$6){print;c1=$2;c2=$5;s1=$3;s2=$6}' > ${prefix}.allValidPairs
    else
        echo "Removing the duplicates ..."
        cat ${vpairs_sorted} | awk -F"\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=$2 || c2!=$5 || s1!=$3 || s2!=$6){print;c1=$2;c2=$5;s1=$3;s2=$6}' > ${prefix}.allValidPairs
    fi

    ## clean
    /bin/rm -rf ${vpairs_sorted}
else
    cat ${vpairs} > ${prefix}.allValidPairs
fi

echo -e -n "valid_interaction\t" > ${prefix}_allValidPairs.mergestat
cat ${vpairs} | wc -l >> ${prefix}_allValidPairs.mergestat
echo -e -n "valid_interaction_rmdup\t" >> ${prefix}_allValidPairs.mergestat
cat ${prefix}.allValidPairs | wc -l >> ${prefix}_allValidPairs.mergestat

## Count short range (<20000) vs long range contacts
awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} $2 == $5{cis=cis+1; d=$6>$3?$6-$3:$3-$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} $2!=$5{trans=trans+1}END{print "trans_interaction\t"trans"\ncis_interaction\t"cis"\ncis_shortRange\t"sr"\ncis_longRange\t"lr}' ${prefix}.allValidPairs >> ${prefix}_allValidPairs.mergestat

## clean
/bin/rm -rf ./tmp/
