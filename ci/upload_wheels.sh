# Modified from numpy's https://github.com/numpy/numpy/blob/main/tools/wheels/upload_wheels.sh

set_upload_vars() {
    echo "IS_PUSH is $IS_PUSH"
    echo "IS_SCHEDULE_DISPATCH is $IS_SCHEDULE_DISPATCH"
    if [[ "$IS_PUSH" == "true" ]]; then
        echo push and tag event
        export ANACONDA_ORG="multibuild-wheels-staging"
        export TOKEN="$PANDAS_STAGING_UPLOAD_TOKEN"
        export ANACONDA_UPLOAD="true"
    elif [[ "$IS_SCHEDULE_DISPATCH" == "true" ]]; then
        echo scheduled or dispatched event
        export ANACONDA_ORG="scipy-wheels-nightly"
        export TOKEN="$PANDAS_NIGHTLY_UPLOAD_TOKEN"
        export ANACONDA_UPLOAD="true"
    else
        echo non-dispatch event
        export ANACONDA_UPLOAD="false"
    fi
}
upload_wheels() {
    echo ${PWD}
    if [[ ${ANACONDA_UPLOAD} == true ]]; then
        if [ -z ${TOKEN} ]; then
            echo no token set, not uploading
        else
            # sdists are located under dist folder when built through setup.py
            if compgen -G "./dist/*.gz"; then
                echo "Found sdist"
                anaconda -q -t ${TOKEN} upload --skip -u ${ANACONDA_ORG} ./dist/*.gz
            elif compgen -G "./wheelhouse/*.whl"; then
                echo "Found wheel"
                anaconda -q -t ${TOKEN} upload --skip -u ${ANACONDA_ORG} ./wheelhouse/*.whl
            else
                echo "Files do not exist"
                return 1
            fi
            echo "PyPI-style index: https://pypi.anaconda.org/$ANACONDA_ORG/simple"
        fi
    fi
}
