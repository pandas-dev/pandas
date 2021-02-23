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

    # Pass the Github username as the build variable
    read -rp "Github username: " gh_username
    docker build --tag pandas-"$USER"-env --build-arg gh_username=$gh_username .

    run_container

fi
