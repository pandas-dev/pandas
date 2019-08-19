# Working with git

Now that you have an issue you want to fix, enhancement to add, or
documentation to improve, you need to learn how to work with GitHub and
the *pandas* code base.

## Version control, Git, and GitHub

To the new user, working with Git is one of the more daunting aspects of
contributing to *pandas*. It can very quickly become overwhelming, but
sticking to the guidelines below will help keep the process
straightforward and mostly trouble free. As always, if you are having
difficulties please feel free to ask for help.

The code is hosted on
[GitHub](https://www.github.com/pandas-dev/pandas). To contribute you
will need to sign up for a [free GitHub account](https://github.com/signup/free).
We use [Git](http://git-scm.com/) for version control to allow many people to
work together on the project.

Some great resources for learning Git:

-   the [GitHub help pages](http://help.github.com/).
-   the [NumPy documentation](http://docs.scipy.org/doc/numpy/dev/index.html).
-   Matthew Brett's [Pydagogue](http://matthew-brett.github.com/pydagogue/).

## Getting started with Git

[GitHub has instructions](http://help.github.com/set-up-git-redirect)
for installing git, setting up your SSH key, and configuring git. All
these steps need to be completed before you can work seamlessly between
your local repository and GitHub.

## Forking

You will need your own fork to work on the code. Go to the [pandas
project page](https://github.com/pandas-dev/pandas) and hit the
`Fork` button. You will want to clone your fork to your machine:

    git clone https://github.com/your-user-name/pandas.git pandas-yourname
    cd pandas-yourname
    git remote add upstream https://github.com/pandas-dev/pandas.git

This creates the directory `pandas-yourname` and connects
your repository to the upstream (main project) *pandas* repository.

## Creating a branch

You want your master branch to reflect only production-ready code, so
create a feature branch for making your changes. For example:

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to:

    git checkout -b shiny-new-feature

This changes your working directory to the shiny-new-feature branch.
Keep any changes in this branch specific to one bug or feature so it is
clear what the branch brings to *pandas*. You can have many
shiny-new-features and switch in between them using the git checkout
command.

When creating this branch, make sure your master branch is up to date
with the latest upstream master version. To update your local master
branch, you can do:

    git checkout master
    git pull upstream master --ff-only

When you want to update the feature branch with changes in master after
you created the branch, check the section on
updating a PR.
