set -e

if [[ ! $(git tag) ]]; then
    echo "No git tags in clone, please sync your git tags with upstream using:"
    echo "    git fetch --tags upstream"
    echo "    git push --tags origin"
    echo ""
    echo "If the issue persists, the clone depth needs to be increased in .travis.yml"
    exit 1
fi

# This will error if there are no tags and we omit --always
DESCRIPTION=$(git describe --long --tags)
echo "$DESCRIPTION"

if [[ "$DESCRIPTION" == *"untagged"* ]]; then
    echo "Unable to determine most recent tag, aborting build"
    exit 1
else
    if [[ "$DESCRIPTION" != *"g"* ]]; then
	# A good description will have the hash prefixed by g, a bad one will be
	# just the hash
	echo "Unable to determine most recent tag, aborting build"
	exit 1
    else
	echo "$(git tag)"
    fi
fi
