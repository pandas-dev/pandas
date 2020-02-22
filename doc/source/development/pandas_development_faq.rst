.. _pandas_development_faq:

{{ header }}

======================
Pandas Development FAQ
======================

.. contents:: Table of contents:
   :local:

Purpose
=======

Based on https://github.com/pydata/pandas/pull/4404#issuecomment-22864665 this
wiki page gathers oft-asked questions/comments from contributors to make the
contribution process a bit less painful.

The aim is to make it easier for

* Core developers to give advice & accept new code contributions.
* New contributors to find an easier way in for quick and efficient bug-fixes
  or feature additions

While some questions/comments/advice may be applicable to general programming,
these are things that directly relate to ``pandas`` development.

* `**PR == pull request** <https://help.github.com/articles/using-pull-requests>`_
* **core developer:** A person contributing on very high frequency & who is
  familiar with the code base and development process of ``pandas``.
* **contributors:** The occasional contributor, maybe from a specific domain,
  contributes bug fixes, features or documentation with low frequency, may not
  be an every day programmer (e.g. programming scientists or engineer using
  pandas for data processing) and looks at things from an end-user perspective.

Pandas Development & Release Process
====================================

Testing
-------

**Q:** What are some recommendations for writing unit tests?

**A:** Your test should be self-contained. That is, it should test preferably a
single thing, e.g., a method that you've added to the ``DataFrame`` class. Your
test function/method should start with ``test_`` and the rest of the name should
be related to whatever functionality you're testing, like
``test_replace_with_dict_regex``.

**Q:** Help! I can't get the tests to run!

**A:** You probably either have multiple Python versions installed and there's
an ABI (application binary interface) issue or you forgot to build the extension
modules in place. The latter can be done with

.. code-block:: shell

    python setup.py build_ext --inplace

from the ``pandas`` directory.

Travis
------

**Q:** Where do I need to change the settings in my GitHub configuration and/or
Travis configuration for the Travis to start builds from my fork?

**A:** To be filled out.

**Q:** Why do I need a Travis file in my repo if it's already in the head
repository?

**A:** Because we're not using subversion. Okay, seriously, it's because as far
as ``git`` is concerned *your* repository is the *only* one that exists. There's
really no such thing as a "head" repository in the eyes of ``git``, those are
concepts that we impose on it to make collaboration more effective and easier.
This is one of the nice aspects of
`distributed version control <http://en.wikipedia.org/wiki/Distributed_revision_control>`_.

Documentation
-------------

**Q:** Does Travis build documentation?

**A:** Currently, no. There are some issues surrounding Sphinx error reporting.
We are investigating ways to solve this problem.

Workflow
--------

* What is a typical workflow on my local fork?
* Shall I work in a virtualenvironment?
* Shall I work in a virtualenvironment and then copy my changes over into a
  clean local fork of my own repo?

**Q:** Who will be responsible for evaluating my PR?

**A:** Technically, anyone with push rights to the ``pydata/pandas`` can
evaluate it. In practice, there are a handful of people who are constantly
watching the ``pandas`` repo for new PRs, so most likely it'll be one of them
that evaluates it. I'm not going to list names, but it's not that hard to figure
out...

Criteria for PR
---------------

**Q:** What are the criteria for acceptance of a PR?

**A:** First and foremost, your fix **must not break any existing
functionality**, one indicator of this is that your Travis build passes. Second,
just give it some time. Everyone is busy and @wesm has not (yet?) amassed a
``pandas`` development army.

**Q:** Do I need to open an issue first?

**A:** Not necessarily. If you want to submit a documentation change, e.g., a
typo fix, then opening an issue is not necessary.

Coding Style
------------

**Q:** What level of commenting is accepted?

**A:** The common sense level. Don't overdo it on the comments, and make sure
if you *do* comment that your comments explain *what* your code is doing, not
*how* it is doing it (that's what code is for).

Obligatory example:

BAD:

.. code-block:: python

    # increment i
    i = int

    i += 1

GOOD:

.. code-block:: python

    # add a person to the person count
    i = int

    i += 1

Debugging
---------

**Q:** How can I debug without adding loads of ``print`` statements/calls
everywhere?

**A:** You can use the Python standard library's ``pdb`` and set a breakpoint.
Put ``import pdb; pdb.set_trace()`` at the line where you want to stop.
`ipdb <https://github.com/gotcha/ipdb>`_ is ``pdb`` with tab-completion and a
few other bells and whistles, making debugging less painful. There's also
`ipdbplugin <https://github.com/flavioamieiro/nose-ipdb>`_ which allows you to
drop into ``ipdb`` from `nose <https://github.com/nose-devs/nose>`_ when a test
fails via

.. code-block:: shell

    nosetests --ipdb # or --ipdb-failures

**Q:** Would a logging hook be a solution?

**A:** That's probably a bit overkill. See the suggestions above.

Pandas Library
==============

Source Comments
---------------

* It would be nice to add more source comments to quickly understand the context
  when chiming in to fix an issue

Testing
-------

**Q:** Why don't test functions have a docstring?

**A:** If your tests are self-contained and aren't
`sprawling ecosystems of spaghetti <http://cdn.memegenerator.net/instances/250x250/26336623.jpg>`_
then having a docstring is redundant. Also, the test name is usually (and
should be!) very descriptive. Remember there's no character limit for variable
names. We're not using FORTRAN.

**Q:** ``DataFrame`` and other ``pandas`` objects often many properties/methods.
What is the level of detail that I should consider when I'm writing my test(s)?

**A:** See the previous question/answer. Strive to test one and only one thing.
You could even separate out your tests by their formal parameters if you want
things to be *really* self-contained.

**Q:** Should I consider possible corner cases of my implementation?

**A:** The answer is a resounding **YES**! In some cases you may come across
something that is very pathological. In those cases you should ask a core
developer.

Complexity
----------

* Some modules (e.g. io/parsers.py) seem to have grown into very high complexity.
  It is very time consuming to find out what is done where just for fixing a
  small bug.
* a splitting into several modules would be good
* more in-code comments telling why something is done and under which condition
  and for what expected result.


Docstrings
----------

* even internal functions shall have a simple 1-line docstring
