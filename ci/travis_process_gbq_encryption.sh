#!/bin/bash

source ci/travis_gbq_config.txt

if [[ -n ${!TRAVIS_IV_ENV} ]]; then
    openssl aes-256-cbc -K ${!TRAVIS_KEY_ENV} -iv ${!TRAVIS_IV_ENV} \
    -in ci/travis_gbq.json.enc -out ci/travis_gbq.json -d;
    export GBQ_PROJECT_ID=$GBQ_PROJECT_ID;
    echo 'Successfully decrypted gbq credentials'
fi

