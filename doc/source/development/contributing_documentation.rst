.. _contributing:

{{ header }}

**********************
Contributing to pandas
**********************

.. contents:: Table of contents:
   :local:

Where to start?
===============

All contributions, bug reports, bug fixes, documentation improvements,
enhancements, and ideas are welcome.

If you are brand new to pandas or open-source development, we recommend going
through the `GitHub "issues" tab <https://github.com/pandas-dev/pandas/issues>`_
to find issues that interest you. There are a number of issues listed under `Docs
<https://github.com/pandas-dev/pandas/issues?labels=Docs&sort=updated&state=open>`_
and `good first issue
<https://github.com/pandas-dev/pandas/issues?labels=good+first+issue&sort=updated&state=open>`_
where you could start out. Once you've found an interesting issue, you can
return here to get your development environment setup.

When you start working on an issue, it's a good idea to assign the issue to yourself,
so nobody else duplicates the work on it. GitHub restricts assigning issues to maintainers
of the project only. In most projects, and until recently in pandas, contributors added a
comment letting others know they are working on an issue. While this is ok, you need to
check each issue individually, and it's not possible to find the unassigned ones.

For this reason, we implemented a workaround consisting of adding a comment with the exact
text ``take``. When you do it, a GitHub action will automatically assign you the issue
(this will take seconds, and may require refreshing the page to see it).
By doing this, it's possible to filter the list of issues and find only the unassigned ones.

So, a good way to find an issue to start contributing to pandas is to check the list of
`unassigned good first issues <https://github.com/pandas-dev/pandas/issues?q=is%3Aopen+is%3Aissue+label%3A%22good+first+issue%22+no%3Aassignee>`_
and assign yourself one you like by writing a comment with the exact text ``take``.

If for whatever reason you are not able to continue working with the issue, please try to
unassign it, so other people know it's available again. You can check the list of
assigned issues, since people may not be working in them anymore. If you want to work on one
that is assigned, feel free to kindly ask the current assignee if you can take it
(please allow at least a week of inactivity before considering work in the issue discontinued).

Feel free to ask questions on the `mailing list
<https://groups.google.com/forum/?fromgroups#!forum/pydata>`_ or on `Gitter <https://gitter.im/pydata/pandas/>`_.

.. _contributing.bug_reports:

Bug reports and enhancement requests
====================================

Bug reports are an important part of making pandas more stable. Having a complete bug report
will allow others to reproduce the bug and provide insight into fixing. See
`this stackoverflow article <https://stackoverflow.com/help/mcve>`_ and
`this blogpost <https://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports>`_
for tips on writing a good bug report.

Trying the bug-producing code out on the *master* branch is often a worthwhile exercise
to confirm the bug still exists. It is also worth searching existing bug reports and pull requests
to see if the issue has already been reported and/or fixed.

Bug reports must:

#. Include a short, self-contained Python snippet reproducing the problem.
   You can format the code nicely by using `GitHub Flavored Markdown
   <https://github.github.com/github-flavored-markdown/>`_::

      ```python
      >>> from pandas import DataFrame
      >>> df = DataFrame(...)
      ...
      ```

#. Include the full version string of pandas and its dependencies. You can use the built-in function::

      >>> import pandas as pd
      >>> pd.show_versions()

#. Explain why the current behavior is wrong/not desired and what you expect instead.

The issue will then show up to the pandas community and be open to comments/ideas from others.

.. _contributing.github:

Working with the code
=====================

Now that you have an issue you want to fix, enhancement to add, or documentation to improve,
you need to learn how to work with GitHub and the pandas code base.

.. _contributing.version_control:

Version control, Git, and GitHub
--------------------------------

To the new user, working with Git is one of the more daunting aspects of contributing to pandas.
It can very quickly become overwhelming, but sticking to the guidelines below will help keep the process
straightforward and mostly trouble free.  As always, if you are having difficulties please
feel free to ask for help.

The code is hosted on `GitHub <https://www.github.com/pandas-dev/pandas>`_. To
contribute you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_. We use `Git <https://git-scm.com/>`_ for
version control to allow many people to work together on the project.

Some great resources for learning Git:

* the `GitHub help pages <https://help.github.com/>`_.
* the `NumPy's documentation <https://numpy.org/doc/stable/dev/index.html>`_.
* Matthew Brett's `Pydagogue <https://matthew-brett.github.com/pydagogue/>`_.

Getting started with Git
------------------------

`GitHub has instructions <https://help.github.com/set-up-git-redirect>`__ for installing git,
setting up your SSH key, and configuring git.  All these steps need to be completed before
you can work seamlessly between your local repository and GitHub.

.. _contributing.forking:

Forking
-------

You will need your own fork to work on the code. Go to the `pandas project
page <https://github.com/pandas-dev/pandas>`_ and hit the ``Fork`` button. You will
want to clone your fork to your machine::

    git clone https://github.com/your-user-name/pandas.git pandas-yourname
    cd pandas-yourname
    git remote add upstream https://github.com/pandas-dev/pandas.git

This creates the directory ``pandas-yourname`` and connects your repository to
the upstream (main project) *pandas* repository.

Note that performing a shallow clone (with ``--depth==N``, for some ``N`` greater
or equal to 1) might break some tests and features as ``pd.show_versions()``
as the version number cannot be computed anymore.

.. _contributing.dev_env:

Creating a development environment
----------------------------------

To test out code changes, you'll need to build pandas from source, which
requires a C/C++ compiler and Python environment. If you're making documentation
changes, you can skip to :ref:`contributing.documentation` but if you skip
creating the development environment you won't be able to build the documentation
locally before pushing your changes.

Using a Docker container
~~~~~~~~~~~~~~~~~~~~~~~~

Instead of manually setting up a development environment, you can use `Docker
<https://docs.docker.com/get-docker/>`_ to automatically create the environment with just several
commands. pandas provides a ``DockerFile`` in the root directory to build a Docker image
with a full pandas development environment.

**Docker Commands**

Pass your GitHub username in the ``DockerFile`` to use your own fork::

    # Build the image pandas-yourname-env
    docker build --tag pandas-yourname-env .
    # Run a container and bind your local forked repo, pandas-yourname, to the container
    docker run -it --rm -v path-to-pandas-yourname:/home/pandas-yourname pandas-yourname-env

Even easier, you can integrate Docker with the following IDEs:

**Visual Studio Code**

You can use the DockerFile to launch a remote session with Visual Studio Code,
a popular free IDE, using the ``.devcontainer.json`` file.
See https://code.visualstudio.com/docs/remote/containers for details.

**PyCharm (Professional)**

Enable Docker support and use the Services tool window to build and manage images as well as
run and interact with containers.
See https://www.jetbrains.com/help/pycharm/docker.html for details.

Note that you might need to rebuild the C extensions if/when you merge with upstream/master using::

    python setup.py build_ext -j 4
    
Creating a branch
-----------------

You want your master branch to reflect only production-ready code, so create a
feature branch for making your changes. For example::

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to::

    git checkout -b shiny-new-feature

This changes your working directory to the shiny-new-feature branch.  Keep any
changes in this branch specific to one bug or feature so it is clear
what the branch brings to pandas. You can have many shiny-new-features
and switch in between them using the git checkout command.

When creating this branch, make sure your master branch is up to date with
the latest upstream master version. To update your local master branch, you
can do::

    git checkout master
    git pull upstream master --ff-only

When you want to update the feature branch with changes in master after
you created the branch, check the section on
:ref:`updating a PR <contributing.update-pr>`.


.. _contributing.documentation:

Contributing to the documentation
=================================

Contributing to the documentation benefits everyone who uses pandas.
We encourage you to help us improve the documentation, and
you don't have to be an expert on pandas to do so! In fact,
there are sections of the docs that are worse off after being written by
experts. If something in the docs doesn't make sense to you, updating the
relevant section after you figure it out is a great way to ensure it will help
the next person.

.. contents:: Documentation:
   :local:


About the pandas documentation
--------------------------------

The documentation is written in **reStructuredText**, which is almost like writing
in plain English, and built using `Sphinx <https://www.sphinx-doc.org/en/master/>`__. The
Sphinx Documentation has an excellent `introduction to reST
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`__. Review the Sphinx docs to perform more
complex changes to the documentation as well.

Some other important things to know about the docs:

* The pandas documentation consists of two parts: the docstrings in the code
  itself and the docs in this folder ``doc/``.

  The docstrings provide a clear explanation of the usage of the individual
  functions, while the documentation in this folder consists of tutorial-like
  overviews per topic together with some other information (what's new,
  installation, etc).

* The docstrings follow a pandas convention, based on the **Numpy Docstring
  Standard**. Follow the :ref:`pandas docstring guide <docstring>` for detailed
  instructions on how to write a correct docstring.

  .. toctree::
     :maxdepth: 2

     contributing_docstring.rst

* The tutorials make heavy use of the `IPython directive
  <https://matplotlib.org/sampledoc/ipython_directive.html>`_ sphinx extension.
  This directive lets you put code in the documentation which will be run
  during the doc build. For example::

      .. ipython:: python

          x = 2
          x**3

  will be rendered as::

      In [1]: x = 2

      In [2]: x**3
      Out[2]: 8

  Almost all code examples in the docs are run (and the output saved) during the
  doc build. This approach means that code examples will always be up to date,
  but it does make the doc building a bit more complex.

* Our API documentation files in ``doc/source/reference`` house the auto-generated
  documentation from the docstrings. For classes, there are a few subtleties
  around controlling which methods and attributes have pages auto-generated.

  We have two autosummary templates for classes.

  1. ``_templates/autosummary/class.rst``. Use this when you want to
     automatically generate a page for every public method and attribute on the
     class. The ``Attributes`` and ``Methods`` sections will be automatically
     added to the class' rendered documentation by numpydoc. See ``DataFrame``
     for an example.

  2. ``_templates/autosummary/class_without_autosummary``. Use this when you
     want to pick a subset of methods / attributes to auto-generate pages for.
     When using this template, you should include an ``Attributes`` and
     ``Methods`` section in the class docstring. See ``CategoricalIndex`` for an
     example.

  Every method should be included in a ``toctree`` in one of the documentation files in
  ``doc/source/reference``, else Sphinx
  will emit a warning.

.. note::

    The ``.rst`` files are used to automatically generate Markdown and HTML versions
    of the docs. For this reason, please do not edit ``CONTRIBUTING.md`` directly,
    but instead make any changes to ``doc/source/development/contributing.rst``. Then, to
    generate ``CONTRIBUTING.md``, use `pandoc <https://johnmacfarlane.net/pandoc/>`_
    with the following command::

      pandoc doc/source/development/contributing.rst -t markdown_github > CONTRIBUTING.md

The utility script ``scripts/validate_docstrings.py`` can be used to get a csv
summary of the API documentation. And also validate common errors in the docstring
of a specific class, function or method. The summary also compares the list of
methods documented in the files in ``doc/source/reference`` (which is used to generate
the `API Reference <https://pandas.pydata.org/pandas-docs/stable/api.html>`_ page)
and the actual public methods.
This will identify methods documented in ``doc/source/reference`` that are not actually
class methods, and existing methods that are not documented in ``doc/source/reference``.


Updating a pandas docstring
-----------------------------

When improving a single function or method's docstring, it is not necessarily
needed to build the full documentation (see next section).
However, there is a script that checks a docstring (for example for the ``DataFrame.mean`` method)::

    python scripts/validate_docstrings.py pandas.DataFrame.mean

This script will indicate some formatting errors if present, and will also
run and test the examples included in the docstring.
Check the :ref:`pandas docstring guide <docstring>` for a detailed guide
on how to format the docstring.

The examples in the docstring ('doctests') must be valid Python code,
that in a deterministic way returns the presented output, and that can be
copied and run by users. This can be checked with the script above, and is
also tested on Travis. A failing doctest will be a blocker for merging a PR.
Check the :ref:`examples <docstring.examples>` section in the docstring guide
for some tips and tricks to get the doctests passing.

When doing a PR with a docstring update, it is good to post the
output of the validation script in a comment on github.


How to build the pandas documentation
---------------------------------------

Requirements
~~~~~~~~~~~~

First, you need to have a development environment to be able to build pandas
(see the docs on :ref:`creating a development environment above <contributing.dev_env>`).

Building the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

So how do you build the docs? Navigate to your local
``doc/`` directory in the console and run::

    python make.py html

Then you can find the HTML output in the folder ``doc/build/html/``.

The first time you build the docs, it will take quite a while because it has to run
all the code examples and build all the generated docstring pages. In subsequent
evocations, sphinx will try to only build the pages that have been modified.

If you want to do a full clean build, do::

    python make.py clean
    python make.py html

You can tell ``make.py`` to compile only a single section of the docs, greatly
reducing the turn-around time for checking your changes.

::

    # omit autosummary and API section
    python make.py clean
    python make.py --no-api

    # compile the docs with only a single section, relative to the "source" folder.
    # For example, compiling only this guide (doc/source/development/contributing.rst)
    python make.py clean
    python make.py --single development/contributing.rst

    # compile the reference docs for a single function
    python make.py clean
    python make.py --single pandas.DataFrame.join

For comparison, a full documentation build may take 15 minutes, but a single
section may take 15 seconds. Subsequent builds, which only process portions
you have changed, will be faster.

You can also specify to use multiple cores to speed up the documentation build::

    python make.py html --num-jobs 4

Open the following file in a web browser to see the full documentation you
just built::

    doc/build/html/index.html

And you'll have the satisfaction of seeing your new and improved documentation!

.. _contributing.dev_docs:

Building master branch documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When pull requests are merged into the pandas ``master`` branch, the main parts of
the documentation are also built by Travis-CI. These docs are then hosted `here
<https://pandas.pydata.org/docs/dev/>`__, see also
the :ref:`Continuous Integration <contributing.ci>` section.

Contributing your changes to pandas
=====================================

.. _contributing.commit-code:

Committing your code
--------------------

Keep style fixes to a separate commit to make your pull request more readable.

Once you've made changes, you can see them by typing::

    git status

If you have created a new file, it is not being tracked by git. Add it by typing::

    git add path/to/file-to-be-added.py

Doing 'git status' again should give something like::

    # On branch shiny-new-feature
    #
    #       modified:   /relative/path/to/file-you-added.py
    #

Finally, commit your changes to your local repository with an explanatory message. pandas
uses a convention for commit message prefixes and layout.  Here are
some common prefixes along with general guidelines for when to use them:

* ENH: Enhancement, new functionality
* BUG: Bug fix
* DOC: Additions/updates to documentation
* TST: Additions/updates to tests
* BLD: Updates to the build process/scripts
* PERF: Performance improvement
* TYP: Type annotations
* CLN: Code cleanup

The following defines how a commit message should be structured.  Please reference the
relevant GitHub issues in your commit message using GH1234 or #1234.  Either style
is fine, but the former is generally preferred:

* a subject line with ``< 80`` chars.
* One blank line.
* Optionally, a commit message body.

Now you can commit your changes in your local repository::

    git commit -m

.. _contributing.push-code:

Pushing your changes
--------------------

When you want your changes to appear publicly on your GitHub page, push your
forked feature branch's commits::

    git push origin shiny-new-feature

Here ``origin`` is the default name given to your remote repository on GitHub.
You can see the remote repositories::

    git remote -v

If you added the upstream repository as described above you will see something
like::

    origin  git@github.com:yourname/pandas.git (fetch)
    origin  git@github.com:yourname/pandas.git (push)
    upstream        git://github.com/pandas-dev/pandas.git (fetch)
    upstream        git://github.com/pandas-dev/pandas.git (push)

Now your code is on GitHub, but it is not yet a part of the pandas project. For that to
happen, a pull request needs to be submitted on GitHub.

Review your code
----------------

When you're ready to ask for a code review, file a pull request. Before you do, once
again make sure that you have followed all the guidelines outlined in this document
regarding code style, tests, performance tests, and documentation. You should also
double check your branch changes against the branch it was based on:

#. Navigate to your repository on GitHub -- https://github.com/your-user-name/pandas
#. Click on ``Branches``
#. Click on the ``Compare`` button for your feature branch
#. Select the ``base`` and ``compare`` branches, if necessary. This will be ``master`` and
   ``shiny-new-feature``, respectively.

Finally, make the pull request
------------------------------

If everything looks good, you are ready to make a pull request.  A pull request is how
code from a local repository becomes available to the GitHub community and can be looked
at and eventually merged into the master version.  This pull request and its associated
changes will eventually be committed to the master branch and available in the next
release.  To submit a pull request:

#. Navigate to your repository on GitHub
#. Click on the ``Pull Request`` button
#. You can then click on ``Commits`` and ``Files Changed`` to make sure everything looks
   okay one last time
#. Write a description of your changes in the ``Preview Discussion`` tab
#. Click ``Send Pull Request``.

This request then goes to the repository maintainers, and they will review
the code.

.. _contributing.update-pr:

Updating your pull request
--------------------------

Based on the review you get on your pull request, you will probably need to make
some changes to the code. In that case, you can make them in your branch,
add a new commit to that branch, push it to GitHub, and the pull request will be
automatically updated.  Pushing them to GitHub again is done by::

    git push origin shiny-new-feature

This will automatically update your pull request with the latest code and restart the
:ref:`Continuous Integration <contributing.ci>` tests.

Another reason you might need to update your pull request is to solve conflicts
with changes that have been merged into the master branch since you opened your
pull request.

To do this, you need to "merge upstream master" in your branch::

    git checkout shiny-new-feature
    git fetch upstream
    git merge upstream/master

If there are no conflicts (or they could be fixed automatically), a file with a
default commit message will open, and you can simply save and quit this file.

If there are merge conflicts, you need to solve those conflicts. See for
example at https://help.github.com/articles/resolving-a-merge-conflict-using-the-command-line/
for an explanation on how to do this.
Once the conflicts are merged and the files where the conflicts were solved are
added, you can run ``git commit`` to save those fixes.

If you have uncommitted changes at the moment you want to update the branch with
master, you will need to ``stash`` them prior to updating (see the
`stash docs <https://git-scm.com/book/en/v2/Git-Tools-Stashing-and-Cleaning>`__).
This will effectively store your changes and they can be reapplied after updating.

After the feature branch has been update locally, you can now update your pull
request by pushing to the branch on GitHub::

    git push origin shiny-new-feature

Autofixing formatting errors
----------------------------

We use several styling checks (e.g. ``black``, ``flake8``, ``isort``) which are run after
you make a pull request. If there is a scenario where any of these checks fail then you
can comment::

    @github-actions pre-commit

on that pull request. This will trigger a workflow which will autofix formatting errors.

Delete your merged branch (optional)
------------------------------------

Once your feature branch is accepted into upstream, you'll probably want to get rid of
the branch. First, merge upstream master into your branch so git knows it is safe to
delete your branch::

    git fetch upstream
    git checkout master
    git merge upstream/master

Then you can do::

    git branch -d shiny-new-feature

Make sure you use a lower-case ``-d``, or else git won't warn you if your feature
branch has not actually been merged.

The branch will still exist on GitHub, so to delete it there do::

    git push origin --delete shiny-new-feature

.. _Gitter: https://gitter.im/pydata/pandas


Tips for a successful pull request
==================================

If you have made it to the `Review your code`_ phase, one of the core contributors may
take a look. Please note however that a handful of people are responsible for reviewing
all of the contributions, which can often lead to bottlenecks.

To improve the chances of your pull request being reviewed, you should:

- **Reference an open issue** for non-trivial changes to clarify the PR's purpose
- **Ensure you have appropriate tests**. These should be the first part of any PR
- **Keep your pull requests as simple as possible**. Larger PRs take longer to review
- **Ensure that CI is in a green state**. Reviewers may not even look otherwise
- **Keep** `Updating your pull request`_, either by request or every few days

