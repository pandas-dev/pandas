#!/bin/bash

while [[ $# > 1 ]]
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

if [ -z "$FROM" ]; then
   FROM=`git tag --sort v:refname | grep -v rc | tail -1`
fi

if [ -z "$TO" ]; then
   TO=""
fi

START=`git log $FROM.. --simplify-by-decoration --pretty="format:%ai %d" | tail -1 | gawk '{ print $1 }'`
END=`git log $TO.. --simplify-by-decoration --pretty="format:%ai %d" | head -1 | gawk '{ print $1 }'`

git log $FROM.. --format='%an#%s' | grep -v Merge > commits

# Include a summary by contributor in release notes:
# cat commits | gawk -F '#' '{ print "- " $1 }' | sort | uniq

echo "Stats since <$FROM> [$START - $END]"
echo ""

AUTHORS=`cat commits | gawk -F '#' '{ print $1 }' | sort | uniq | wc -l`
echo "Number of authors: $AUTHORS"

TCOMMITS=`cat commits | gawk -F '#' '{ print $1 }'| wc -l`
echo "Total commits    : $TCOMMITS"

# Include a summary count of commits included in the release by contributor:
# cat commits | gawk -F '#' '{ print $1 }' | sort | uniq -c | sort -nr

/bin/rm commits
