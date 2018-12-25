#!/bin/bash

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    --from)
    FROM="$2"
    shift # past argument
    ;;
    --to)
    TO="$2"
    shift # past argument
    ;;
    *)
          # unknown option
    ;;
esac
shift # past argument or value
done

if [ -z "${FROM}" ]; then
   FROM=$(git tag --sort v:refname | grep -v rc | tail -1)
fi

if [ -z "${TO}" ]; then
   TO=""
fi

START=$(git log "${FROM}".. --simplify-by-decoration --pretty="format:%ai %d" | tail -1 | gawk '{ print $1 }')
END=$(git log "${TO}".. --simplify-by-decoration --pretty="format:%ai %d" | head -1 | gawk '{ print $1 }')

git log "${FROM}".. --format='%an#%s' | grep -v Merge > commits

# Include a summary by contributor in release notes:
# gawk -F '#' '{ print "- " $1 }' commits | sort | uniq

echo -e "Stats since <${FROM}> [${START} - ${END}] \n"

AUTHORS=$(gawk -F '#' '{ print $1 }' commits | sort | uniq | wc -l)
echo -e "Number of authors: ${AUTHORS} \n"

TCOMMITS=$(gawk -F '#' '{ print $1 }' commits | wc -l)
echo -e "Total commits    : ${TCOMMITS} \n"

# Include a summary count of commits included in the release by contributor:
# cat commits | gawk -F '#' '{ print $1 }' | sort | uniq -c | sort -nr

/bin/rm commits
