# Submitting your changes to _pandas_

## Committing your code

Keep style fixes to a separate commit to make your pull request more
readable.

Once you've made changes, you can see them by typing:

    git status

If you have created a new file, it is not being tracked by git. Add it
by typing:

    git add path/to/file-to-be-added.py

Doing \'git status\' again should give something like:

    # On branch shiny-new-feature
    #
    #       modified:   /relative/path/to/file-you-added.py
    #

Finally, commit your changes to your local repository with an
explanatory message. *pandas* uses a convention for commit message
prefixes and layout. Here are some common prefixes along with general
guidelines for when to use them:

-   ENH: Enhancement, new functionality
-   BUG: Bug fix
-   DOC: Additions/updates to documentation
-   TST: Additions/updates to tests
-   BLD: Updates to the build process/scripts
-   PERF: Performance improvement
-   CLN: Code cleanup

The following defines how a commit message should be structured. Please
reference the relevant GitHub issues in your commit message using GH1234
or \#1234. Either style is fine, but the former is generally preferred:

-   a subject line with \< 80 chars.
-   One blank line.
-   Optionally, a commit message body.

Now you can commit your changes in your local repository:

    git commit -m

## Pushing your changes

When you want your changes to appear publicly on your GitHub page, push
your forked feature branch\'s commits:

    git push origin shiny-new-feature

Here `origin` is the default name given to your remote repository on
GitHub. You can see the remote repositories:

    git remote -v

If you added the upstream repository as described above you will see
something like:

    origin  git@github.com:yourname/pandas.git (fetch)
    origin  git@github.com:yourname/pandas.git (push)
    upstream        git://github.com/pandas-dev/pandas.git (fetch)
    upstream        git://github.com/pandas-dev/pandas.git (push)

Now your code is on GitHub, but it is not yet a part of the *pandas*
project. For that to happen, a pull request needs to be submitted on
GitHub.

## Review your code

When you\'re ready to ask for a code review, file a pull request. Before
you do, once again make sure that you have followed all the guidelines
outlined in this document regarding code style, tests, performance
tests, and documentation. You should also double check your branch
changes against the branch it was based on:

1.  Navigate to your repository on GitHub: <https://github.com/your-user-name/pandas>
2.  Click on `Branches`
3.  Click on the `Compare` button for your feature branch
4.  Select the `base` and `compare` branches, if necessary. This will be
    `master` and `shiny-new-feature`, respectively.

## Finally, make the pull request

If everything looks good, you are ready to make a pull request. A pull
request is how code from a local repository becomes available to the
GitHub community and can be looked at and eventually merged into the
master version. This pull request and its associated changes will
eventually be committed to the master branch and available in the next
release. To submit a pull request:

1.  Navigate to your repository on GitHub
2.  Click on the `Pull Request` button
3.  You can then click on `Commits` and `Files Changed` to make sure
    everything looks okay one last time
4.  Write a description of your changes in the `Preview Discussion` tab
5.  Click `Send Pull Request`.

This request then goes to the repository maintainers, and they will
review the code.

## Updating your pull request

Based on the review you get on your pull request, you will probably need
to make some changes to the code. In that case, you can make them in
your branch, add a new commit to that branch, push it to GitHub, and the
pull request will be automatically updated. Pushing them to GitHub again
is done by:

    git push origin shiny-new-feature

This will automatically update your pull request with the latest code
and restart the Continuous Integration tests.

Another reason you might need to update your pull request is to solve
conflicts with changes that have been merged into the master branch
since you opened your pull request.

To do this, you need to \"merge upstream master\" in your branch:

    git checkout shiny-new-feature
    git fetch upstream
    git merge upstream/master

If there are no conflicts (or they could be fixed automatically), a file
with a default commit message will open, and you can simply save and
quit this file.

If there are merge conflicts, you need to solve those conflicts. See for
example at
<https://help.github.com/articles/resolving-a-merge-conflict-using-the-command-line/>
for an explanation on how to do this. Once the conflicts are merged and
the files where the conflicts were solved are added, you can run
`git commit` to save those fixes.

If you have uncommitted changes at the moment you want to update the
branch with master, you will need to `stash` them prior to updating (see
the [stash docs](https://git-scm.com/book/en/v2/Git-Tools-Stashing-and-Cleaning)).
This will effectively store your changes and they can be reapplied after
updating.

After the feature branch has been update locally, you can now update
your pull request by pushing to the branch on GitHub:

    git push origin shiny-new-feature

## Delete your merged branch (optional)

Once your feature branch is accepted into upstream, you\'ll probably
want to get rid of the branch. First, merge upstream master into your
branch so git knows it is safe to delete your branch:

    git fetch upstream
    git checkout master
    git merge upstream/master

Then you can do:

    git branch -d shiny-new-feature

Make sure you use a lower-case `-d`, or else git won\'t warn you if your
feature branch has not actually been merged.

The branch will still exist on GitHub, so to delete it there do:

    git push origin --delete shiny-new-feature
