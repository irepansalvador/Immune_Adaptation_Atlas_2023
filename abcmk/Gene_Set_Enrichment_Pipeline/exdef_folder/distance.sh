#!/bin/bash
# set -e


usage() { echo "Usage: $0 [-c <case file>] [-a <annotation file>]" 1>&2; exit 1; }

while getopts ":c:a:" o; do
    case "${o}" in
        c)
            c=${OPTARG}
            ;;
        a)
            a=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${c}" ] || [ -z "${a}" ]; then
    usage
fi


####################

RND=$RANDOM

#fgrep -f ${c} ${a} | awk -F'\t' 'BEGIN {OFS = FS} {print $1,int(($2+$3)/2),int(($2+$3)/2+1),$4,$5}' | sort -k1,2V >  tmp_${RND}
fgrep -f <(fgrep yes ${c} | cut -d' ' -f1) ${a} | awk -F'\t' 'BEGIN {OFS = FS} {print $1,int(($2+$3)/2),int(($2+$3)/2+1),$4,$5}' | sort -k1,2V >  tmp_${RND}

for i in `cut -f1 ${a} | sort -uV`
do

    genes=$(grep -P "^${i}\t" tmp_${RND})

    if [[ -z ${genes} ]];then
        continue
    else
        # echo ${i}
        closestBed -d -a <(grep -P "^${i}\t" ${a} | awk -F'\t' 'BEGIN {OFS = FS} {print $1,int(($2+$3)/2),int(($2+$3)/2+1),$4,$5}' | sort -k2,3n) -b <(grep -P "^${i}\t" tmp_${RND}) | awk '{print $4" "$11}' >> distance_${RND}
    fi
done

# echo "here"
sort distance_${RND} > params/distance_`basename ${c}`

rm tmp_${RND} distance_${RND}
