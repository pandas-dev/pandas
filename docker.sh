#!/bin/bash
#
# Set up a development environment with Docker.

run_container () {
    # Run a container and bind the local forked repo to it
    docker run -it --rm -v "$(pwd)":/home/pandas-"$USER" pandas-"$USER"-env
}

# Check if pandas image already exists
docker image ls | grep "pandas-$USER-env" &> /dev/null

if [[ $? == 0 ]]; then

    run_container

else

    # Pass the Github username in the Dockerfile
    read -rp "Github username: " gh_username
    sed -i "s/gh_username=pandas-dev/gh_username=$gh_username/" Dockerfile

    docker build --tag pandas-"$USER"-env .

    # Revert change made to Dockerfile. This will prevent the Dockerfile
    # from being commited with the contributor's username in it if
    # the contributor eventually runs "git commit -am".
    sed -i "s/$gh_username/pandas-dev/" Dockerfile

    run_container
fi
