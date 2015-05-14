.. _contributing:

**********************
Contributing to pandas
**********************

.. contents:: Table of contents:
   :local:

Where to start?
===============

All contributions, bug reports, bug fixes, documentation improvements,
enhancements and ideas are welcome.

If you are simply looking to start working with the *pandas* codebase, navigate to the
`GitHub "issues" tab <https://github.com/pydata/pandas/issues>`_ and start looking through
interesting issues.  There are a number of issues listed under `Docs
<https://github.com/pydata/pandas/issues?labels=Docs&sort=updated&state=open>`_
and `Difficulty Novice
<https://github.com/pydata/pandas/issues?q=is%3Aopen+is%3Aissue+label%3A%22Difficulty+Novice%22>`_
where you could start out.

Or maybe through using *pandas* you have an idea of you own or are looking for something
in the documentation and thinking 'this can be improved'...you can do something
about it!

Feel free to ask questions on `mailing list
<https://groups.google.com/forum/?fromgroups#!forum/pydata>`_

Bug Reports/Enhancement Requests
================================

Bug reports are an important part of making *pandas* more stable.  Having a complete bug report
will allow others to reproduce the bug and provide insight into fixing.  Since many versions of
*pandas* are supported, knowing version information will also identify improvements made since
previous versions.  Often trying the bug-producing code out on the *master* branch is a worthwhile exercise
to confirm the bug still exists.  It is also worth searching existing bug reports and pull requests
to see if the issue has already been reported and/or fixed.

Bug reports must:

#. Include a short, self-contained Python snippet reproducing the problem.
   You can have the code formatted nicely by using `GitHub Flavored Markdown
   <http://github.github.com/github-flavored-markdown/>`_: ::

      ```python
      >>> from pandas import DataFrame
      >>> df = DataFrame(...)
      ...
      ```

#. Include the full version string of *pandas* and its dependencies. In recent (>0.12) versions
   of *pandas* you can use a built in function: ::

      >>> from pandas.util.print_versions import show_versions
      >>> show_versions()

   and in 0.13.1 onwards: ::

      >>> pd.show_versions()

#. Explain why the current behavior is wrong/not desired and what you expect instead.

The issue will then show up to the *pandas* community and be open to comments/ideas from others.

Working with the code
=====================

Now that you have an issue you want to fix, enhancement to add, or documentation to improve,
you need to learn how to work with GitHub and the *pandas* code base.

Version Control, Git, and GitHub
--------------------------------

To the new user, working with Git is one of the more daunting aspects of contributing to *pandas*.
It can very quickly become overwhelming, but sticking to the guidelines below will make the process
straightforward and will work without much trouble.  As always, if you are having difficulties please
feel free to ask for help.

The code is hosted on `GitHub <https://www.github.com/pydata/pandas>`_. To
contribute you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_. We use `Git <http://git-scm.com/>`_ for
version control to allow many people to work together on the project.

Some great resources for learning git:

* the `GitHub help pages <http://help.github.com/>`_.
* the `NumPy's documentation <http://docs.scipy.org/doc/numpy/dev/index.html>`_.
* Matthew Brett's `Pydagogue <http://matthew-brett.github.com/pydagogue/>`_.

Getting Started with Git
------------------------

`GitHub has instructions <http://help.github.com/set-up-git-redirect>`__ for installing git,
setting up your SSH key, and configuring git.  All these steps need to be completed before
working seamlessly with your local repository and GitHub.

.. _contributing.forking:

Forking
-------

You will need your own fork to work on the code. Go to the `pandas project
page <https://github.com/pydata/pandas>`_ and hit the *fork* button. You will
want to clone your fork to your machine: ::

    git clone git@github.com:your-user-name/pandas.git pandas-yourname
    cd pandas-yourname
    git remote add upstream git://github.com/pydata/pandas.git

This creates the directory `pandas-yourname` and connects your repository to
the upstream (main project) *pandas* repository.

The testing suite will run automatically on Travis-CI once your Pull Request is
submitted.  However, if you wish to run the test suite on a branch prior to
submitting the Pull Request, then Travis-CI needs to be hooked up to your
GitHub repository.  Instructions are for doing so are `here
<http://about.travis-ci.org/docs/user/getting-started/>`__.

Creating a Branch
-----------------

You want your master branch to reflect only production-ready code, so create a
feature branch for making your changes. For example::

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to::

    git checkout -b shiny-new-feature

This changes your working directory to the shiny-new-feature branch.  Keep any
changes in this branch specific to one bug or feature so it is clear
what the branch brings to *pandas*. You can have many shiny-new-features
and switch in between them using the git checkout command.

To update this branch, you need to retrieve the changes from the master branch::

    git fetch upstream
    git rebase upstream/master

This will replay your commits on top of the lastest pandas git master.  If this
leads to merge conflicts, you must resolve these before submitting your Pull
Request.  If you have uncommitted changes, you will need to `stash` them prior
to updating.  This will effectively store your changes and they can be reapplied
after updating.

.. _contributing.dev_env:

Creating a Development Environment
----------------------------------

An easy way to create a *pandas* development environment is as follows.

- Install either :ref:`Install Anaconda <install-anaconda>` or :ref:`Install miniconda <install-miniconda>`
- Make sure that you have :ref:`cloned the repository <contributing-forking>`
- ``cd`` to the pandas source directory

Tell ``conda`` to create a new environment, named ``pandas_dev``, or any name you would like for this environment by running:

::

      conda create -n pandas_dev --file ci/requirements_dev.txt


For a python 3 environment

::

      conda create -n pandas_dev python=3 --file ci/requirements_dev.txt


If you are on ``windows``, then you will need to install the compiler linkages:

::

      conda install -n pandas_dev libpython

This will create the new environment, and not touch any of your existing environments, nor any existing python installation. It will install all of the basic dependencies of *pandas*, as well as the development and testing tools. If you would like to install other dependencies, you can install them as follows:

::

      conda install -n pandas_dev -c pandas pytables scipy

To install *all* pandas dependencies you can do the following:

::

      conda install -n pandas_dev -c pandas --file ci/requirements_all.txt

To work in this environment, ``activate`` it as follows:

::

      activate pandas_dev

At which point, the prompt will change to indicate you are in the new development environment.

.. note::

   The above syntax is for ``windows`` environments. To work on ``macosx/linux``, use:

   ::

       source activate pandas_dev

To view your environments:

::

      conda info -e

To return to you home root environment:

::

      deactivate

See the full ``conda`` docs `here
<http://conda.pydata.org/docs>`__.

At this point you can easily do an *in-place* install, as detailed in the next section.

.. _contributing.getting_source:

Making changes
--------------

Before making your code changes, it is often necessary to build the code that was
just checked out.  There are two primary methods of doing this.

#. The best way to develop *pandas* is to build the C extensions in-place by
   running::

      python setup.py build_ext --inplace

   If you startup the Python interpreter in the *pandas* source directory you
   will call the built C extensions

#. Another very common option is to do a ``develop`` install of *pandas*::

      python setup.py develop

   This makes a symbolic link that tells the Python interpreter to import *pandas*
   from your development directory. Thus, you can always be using the development
   version on your system without being inside the clone directory.

Contributing to the documentation
=================================

If you're not the developer type, contributing to the documentation is still
of huge value. You don't even have to be an expert on
*pandas* to do so! Something as simple as rewriting small passages for clarity
as you reference the docs is a simple but effective way to contribute. The
next person to read that passage will be in your debt!

Actually, there are sections of the docs that are worse off by being written
by experts. If something in the docs doesn't make sense to you, updating the
relevant section after you figure it out is a simple way to ensure it will
help the next person.

.. contents:: Documentation:
   :local:


About the pandas documentation
------------------------------

The documentation is written in **reStructuredText**, which is almost like writing
in plain English, and built using `Sphinx <http://sphinx.pocoo.org/>`__. The
Sphinx Documentation has an excellent `introduction to reST
<http://sphinx.pocoo.org/rest.html>`__. Review the Sphinx docs to perform more
complex changes to the documentation as well.

Some other important things to know about the docs:

- The *pandas* documentation consists of two parts: the docstrings in the code
  itself and the docs in this folder ``pandas/doc/``.

  The docstrings provide a clear explanation of the usage of the individual
  functions, while the documentation in this folder consists of tutorial-like
  overviews per topic together with some other information (what's new,
  installation, etc).

- The docstrings follow the **Numpy Docstring Standard** which is used widely
  in the Scientific Python community. This standard specifies the format of
  the different sections of the docstring. See `this document
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
  for a detailed explanation, or look at some of the existing functions to
  extend it in a similar manner.

- The tutorials make heavy use of the `ipython directive
  <http://matplotlib.org/sampledoc/ipython_directive.html>`_ sphinx extension.
  This directive lets you put code in the documentation which will be run
  during the doc build. For example:

  ::

      .. ipython:: python

          x = 2
          x**3

  will be rendered as

  ::

      In [1]: x = 2

      In [2]: x**3
      Out[2]: 8

  This means that almost all code examples in the docs are always run (and the
  output saved) during the doc build. This way, they will always be up to date,
  but it makes the doc building a bit more complex.


How to build the pandas documentation
-------------------------------------

Requirements
~~~~~~~~~~~~

To build the *pandas* docs there are some extra requirements: you will need to
have ``sphinx`` and ``ipython`` installed. `numpydoc
<https://github.com/numpy/numpydoc>`_ is used to parse the docstrings that
follow the Numpy Docstring Standard (see above), but you don't need to install
this because a local copy of ``numpydoc`` is included in the *pandas* source
code.

It is easiest to :ref:`create a development environment <contributing-dev_env>`, then install:

::

      conda install -n pandas_dev sphinx ipython

Furthermore, it is recommended to have all `optional dependencies
<http://pandas.pydata.org/pandas-docs/dev/install.html#optional-dependencies>`_
installed. This is not strictly necessary, but be aware that you will see some error
messages. Because all the code in the documentation is executed during the doc
build, the examples using this optional dependencies will generate errors.
Run ``pd.show_versions()`` to get an overview of the installed version of all
dependencies.

.. warning::

   Sphinx version >= 1.2.2 or the older 1.1.3 is required.

Building the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

So how do you build the docs? Navigate to your local the folder
``pandas/doc/`` directory in the console and run::

    python make.py html

And then you can find the html output in the folder ``pandas/doc/build/html/``.

The first time it will take quite a while, because it has to run all the code
examples in the documentation and build all generated docstring pages.
In subsequent evocations, sphinx will try to only build the pages that have
been modified.

If you want to do a full clean build, do::

    python make.py clean
    python make.py build


Starting with 0.13.1 you can tell ``make.py`` to compile only a single section
of the docs, greatly reducing the turn-around time for checking your changes.
You will be prompted to delete `.rst` files that aren't required.  This is okay
since the prior version can be checked out from git, but make sure to
not commit the file deletions.

::

    #omit autosummary and API section
    python make.py clean
    python make.py --no-api

    # compile the docs with only a single
    # section, that which is in indexing.rst
    python make.py clean
    python make.py --single indexing

For comparison, a full documentation build may take 10 minutes. a ``-no-api`` build
may take 3 minutes and a single section may take 15 seconds.  However, subsequent
builds only process portions you changed.  Now, open the following file in a web
browser to see the full documentation you just built::

    pandas/docs/build/html/index.html

And you'll have the satisfaction of seeing your new and improved documentation!

.. _contributing.dev_docs:

Built Master Branch Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When pull-requests are merged into the pandas *master* branch, the main parts of the documentation are
also built by Travis-CI. These docs are then hosted `here <http://pandas-docs.github.io/pandas-docs-travis>`__.

Contributing to the code base
=============================

.. contents:: Code Base:
   :local:

Code Standards
--------------

*pandas* uses the `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ standard.
There are several tools to ensure you abide by this standard.

We've written a tool to check that your commits are PEP8 great, `pip install pep8radius <https://github.com/hayd/pep8radius>`_.
Look at PEP8 fixes in your branch vs master with::

    pep8radius master --diff` and make these changes with `pep8radius master --diff --in-place`

Alternatively, use `flake8 <http://pypi.python.org/pypi/flake8>`_ tool for checking the style of your code.
Additional standards are outlined on the `code style wiki page <https://github.com/pydata/pandas/wiki/Code-Style-and-Conventions>`_.

Please try to maintain backward-compatibility. *Pandas* has lots of users with lots of existing code, so
don't break it if at all possible.  If you think breakage is required clearly state why
as part of the Pull Request.  Also, be careful when changing method signatures and add
deprecation warnings where needed.

Test-driven Development/Writing Code
------------------------------------

*Pandas* is serious about testing and strongly encourages individuals to embrace `Test-driven Development (TDD)
<http://en.wikipedia.org/wiki/Test-driven_development>`_.
This development process "relies on the repetition of a very short development cycle:
first the developer writes an (initially failing) automated test case that defines a desired
improvement or new function, then produces the minimum amount of code to pass that test."
So, before actually writing any code, you should write your tests.  Often the test can be
taken from the original GitHub issue.  However, it is always worth considering additional
use cases and writing corresponding tests.

Adding tests is one of the most common requests after code is pushed to *pandas*.  It is worth getting
in the habit of writing tests ahead of time so this is never an issue.

Like many packages, *pandas* uses the `Nose testing system
<http://somethingaboutorange.com/mrl/projects/nose/>`_ and the convenient
extensions in `numpy.testing
<http://docs.scipy.org/doc/numpy/reference/routines.testing.html>`_.

Writing tests
~~~~~~~~~~~~~

All tests should go into the *tests* subdirectory of the specific package.
There are probably many examples already there and looking to these for
inspiration is suggested.  If you test requires working with files or
network connectivity there is more information on the `testing page
<https://github.com/pydata/pandas/wiki/Testing>`_ of the wiki.

The ``pandas.util.testing`` module has many special ``assert`` functions that
make it easier to make statements about whether Series or DataFrame objects are
equivalent. The easiest way to verify that your code is correct is to
explicitly construct the result you expect, then compare the actual result to
the expected correct result:

::

    def test_pivot(self):
        data = {
            'index' : ['A', 'B', 'C', 'C', 'B', 'A'],
            'columns' : ['One', 'One', 'One', 'Two', 'Two', 'Two'],
            'values' : [1., 2., 3., 3., 2., 1.]
        }

        frame = DataFrame(data)
        pivoted = frame.pivot(index='index', columns='columns', values='values')

        expected = DataFrame({
            'One' : {'A' : 1., 'B' : 2., 'C' : 3.},
            'Two' : {'A' : 1., 'B' : 2., 'C' : 3.}
        })

        assert_frame_equal(pivoted, expected)

Running the test suite
~~~~~~~~~~~~~~~~~~~~~~

The tests can then be run directly inside your git clone (without having to
install *pandas*) by typing:::

    nosetests pandas

The tests suite is exhaustive and takes around 20 minutes to run.  Often it is
worth running only a subset of tests first around your changes before running the
entire suite.  This is done using one of the following constructs:

::

    nosetests pandas/tests/[test-module].py
    nosetests pandas/tests/[test-module].py:[TestClass]
    nosetests pandas/tests/[test-module].py:[TestClass].[test_method]


Running the performance test suite
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Performance matters and it is worth considering that your code has not introduced
performance regressions.  Currently *pandas* uses the `vbench library <https://github.com/pydata/vbench>`__
to enable easy monitoring of the performance of critical *pandas* operations.
These benchmarks are all found in the ``pandas/vb_suite`` directory.  vbench
currently only works on python2.

To install vbench::

    pip install git+https://github.com/pydata/vbench

Vbench also requires sqlalchemy, gitpython, and psutil which can all be installed
using pip.  If you need to run a benchmark, change your directory to the *pandas* root and run::

    ./test_perf.sh -b master -t HEAD

This will checkout the master revision and run the suite on both master and
your commit.  Running the full test suite can take up to one hour and use up
to 3GB of RAM.  Usually it is sufficient to past a subset of the results in
to the Pull Request to show that the committed changes do not cause unexpected
performance regressions.

You can run specific benchmarks using the *-r* flag which takes a regular expression.

See the `performance testing wiki <https://github.com/pydata/pandas/wiki/Performance-Testing>`_ for information
on how to write a benchmark.

Documenting your code
---------------------

Changes should be reflected in the release notes located in `doc/source/whatsnew/vx.y.z.txt`.
This file contains an ongoing change log for each release.  Add an entry to this file to
document your fix, enhancement or (unavoidable) breaking change.  Make sure to include the
GitHub issue number when adding your entry.

If your code is an enhancement, it is most likely necessary to add usage examples to the
existing documentation.  This can be done following the section regarding documentation.

Contributing your changes to *pandas*
=====================================

Committing your code
--------------------

Keep style fixes to a separate commit to make your PR more readable.

Once you've made changes, you can see them by typing::

    git status

If you've created a new file, it is not being tracked by git. Add it by typing ::

    git add path/to/file-to-be-added.py

Doing 'git status' again should give something like ::

    # On branch shiny-new-feature
    #
    #       modified:   /relative/path/to/file-you-added.py
    #

Finally, commit your changes to your local repository with an explanatory message.  *Pandas*
uses a convention for commit message prefixes and layout.  Here are
some common prefixes along with general guidelines for when to use them:

    * ENH: Enhancement, new functionality
    * BUG: Bug fix
    * DOC: Additions/updates to documentation
    * TST: Additions/updates to tests
    * BLD: Updates to the build process/scripts
    * PERF: Performance improvement
    * CLN: Code cleanup

The following defines how a commit message should be structured.  Please reference the
relevant GitHub issues in your commit message using `GH1234` or `#1234`.  Either style
is fine, but the former is generally preferred:

    * a subject line with `< 80` chars.
    * One blank line.
    * Optionally, a commit message body.

Now you can commit your changes in your local repository::

    git commit -m

If you have multiple commits, it is common to want to combine them into one commit, often
referred to as "squashing" or "rebasing".  This is a common request by package maintainers
when submitting a Pull Request as it maintains a more compact commit history.  To rebase your commits::

    git rebase -i HEAD~#

Where # is the number of commits you want to combine.  Then you can pick the relevant
commit message and discard others.

Pushing your changes
--------------------

When you want your changes to appear publicly on your GitHub page, push your
forked feature branch's commits ::

    git push origin shiny-new-feature

Here `origin` is the default name given to your remote repository on GitHub.
You can see the remote repositories ::

    git remote -v

If you added the upstream repository as described above you will see something
like ::

    origin  git@github.com:yourname/pandas.git (fetch)
    origin  git@github.com:yourname/pandas.git (push)
    upstream        git://github.com/pydata/pandas.git (fetch)
    upstream        git://github.com/pydata/pandas.git (push)

Now your code is on GitHub, but it is not yet a part of the *pandas* project.  For that to
happen, a Pull Request needs to be submitted on GitHub.

Review your code
----------------

When you're ready to ask for a code review, you will file a Pull Request. Before you do,
again make sure you've followed all the guidelines outlined in this document regarding
code style, tests, performance tests, and documentation. You should also double check
your branch changes against the branch it was based off of:

#. Navigate to your repository on GitHub--https://github.com/your-user-name/pandas.
#. Click on `Branches`.
#. Click on the `Compare` button for your feature branch.
#. Select the `base` and `compare` branches, if necessary. This will be `master` and
   `shiny-new-feature`, respectively.

Finally, make the Pull Request
------------------------------

If everything looks good you are ready to make a Pull Request.  A Pull Request is how
code from a local repository becomes available to the GitHub community and can be looked
at and eventually merged into the master version.  This Pull Request and its associated
changes will eventually be committed to the master branch and available in the next
release.  To submit a Pull Request:

#. Navigate to your repository on GitHub.
#. Click on the `Pull Request` button.
#. You can then click on `Commits` and `Files Changed` to make sure everything looks okay one last time.
#. Write a description of your changes in the `Preview Discussion` tab.
#. Click `Send Pull Request`.

This request then appears to the repository maintainers, and they will review
the code. If you need to make more changes, you can make them in
your branch, push them to GitHub, and the pull request will be automatically
updated.  Pushing them to GitHub again is done by::

    git push -f origin shiny-new-feature

This will automatically update your Pull Request with the latest code and restart the Travis-CI tests.

Delete your merged branch (optional)
------------------------------------

Once your feature branch is accepted into upstream, you'll probably want to get rid of
the branch. First, merge upstream master into your branch so git knows it is safe to delete your branch ::

    git fetch upstream
    git checkout master
    git merge upstream/master

Then you can just do::

    git branch -d shiny-new-feature

Make sure you use a lower-case -d, or else git won't warn you if your feature
branch has not actually been merged.

The branch will still exist on GitHub, so to delete it there do ::

    git push origin --delete shiny-new-feature
