#!/bin/bash

GBQ_JSON_FILE=$1

if [[ $# -ne 1 ]]; then
    echo -e "Too few arguments.\nUsage: ./travis_encrypt_gbq.sh "\
    "<gbq-json-credentials-file>"
    exit 1
fi

if [[ $GBQ_JSON_FILE != *.json ]]; then
    echo "ERROR: Expected *.json file"
    exit 1
fi

if [[ ! -f $GBQ_JSON_FILE ]]; then
    echo "ERROR: File $GBQ_JSON_FILE does not exist"
    exit 1
fi

echo "Encrypting $GBQ_JSON_FILE..."
read -d "\n" TRAVIS_KEY TRAVIS_IV <<<$(travis encrypt-file $GBQ_JSON_FILE \
travis_gbq.json.enc -f | grep -o "\w*_iv\|\w*_key");

echo "Adding your secure key to travis_gbq_config.txt ..."
echo -e "TRAVIS_IV_ENV=$TRAVIS_IV\nTRAVIS_KEY_ENV=$TRAVIS_KEY"\
> travis_gbq_config.txt

echo "Done. Removing file $GBQ_JSON_FILE"
rm $GBQ_JSON_FILE

echo -e "Created encrypted credentials file travis_gbq.json.enc.\n"\
     "NOTE: Do NOT commit the *.json file containing your unencrypted" \
     "private key"
