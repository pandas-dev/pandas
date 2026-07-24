Regenerate minified files with::

    npm install -g uglify-js
    for f in $(ls sphinx/search/non-minified-js/); \
    do echo $f && \
    npx uglifyjs sphinx/search/non-minified-js/$f --compress --mangle --output sphinx/search/minified-js/$f; \
    done
