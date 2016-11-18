#!/bin/bash

source ci/travis_gbq_config.txt

if [[ -n ${SERVICE_ACCOUNT_KEY} ]]; then
    echo "${SERVICE_ACCOUNT_KEY}" > ci/travis_gbq.json;
elif [[ -n ${!TRAVIS_IV_ENV} ]]; then
    openssl aes-256-cbc -K ${!TRAVIS_KEY_ENV} -iv ${!TRAVIS_IV_ENV} \
    -in ci/travis_gbq.json.enc -out ci/travis_gbq.json -d;
    export GBQ_PROJECT_ID='pandas-travis';
    echo 'Successfully decrypted gbq credentials'
fi

