#!/bin/sh

LAST=`git tag --sort v:refname | grep -v rc | tail -1`
START=`git log $LAST.. --simplify-by-decoration --pretty="format:%ai %d" | tail -1 | gawk '{ print $1 }'`
END=`git log $LAST.. --simplify-by-decoration --pretty="format:%ai %d" | head -1 | gawk '{ print $1 }'`

git log $LAST.. --format='%an#%s' | grep -v Merge > commits

# Include a summary by contributor in release notes:
# cat commits | gawk -F '#' '{ print "- " $1 }' | sort | uniq

echo "Stats since <$LAST> [$START - $END]"
echo ""

AUTHORS=`cat commits | gawk -F '#' '{ print $1 }' | sort | uniq | wc -l`
echo "Number of authors: $AUTHORS"

TCOMMITS=`cat commits | gawk -F '#' '{ print $1 }'| wc -l`
echo "Total commits    : $TCOMMITS"

# Include a summary count of commits included in the release by contributor:
# cat commits | gawk -F '#' '{ print $1 }' | sort | uniq -c | sort -nr

/bin/rm commits
