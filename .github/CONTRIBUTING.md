# Contributing to pandas

Whether you are a novice or experienced software developer, all contributions and suggestions are welcome!

Our main contributing guide can be found [in this repo](https://github.com/pandas-dev/pandas/blob/master/doc/source/development/contributing.rst) or [on the website](https://pandas.pydata.org/docs/dev/development/contributing.html). If you do not want to read it in its entirety, we will summarize the main ways in which you can contribute and point to relevant sections of that document for further information.

## Getting Started

If you are looking to contribute to the *pandas* codebase, the best place to start is the [GitHub "issues" tab](https://github.com/pandas-dev/pandas/issues). This is also a great place for filing bug reports and making suggestions for ways in which we can improve the code and documentation.

If you have additional questions, feel free to ask them on the [mailing list](https://groups.google.com/forum/?fromgroups#!forum/pydata) or on [Gitter](https://gitter.im/pydata/pandas). Further information can also be found in the "[Where to start?](https://github.com/pandas-dev/pandas/blob/master/doc/source/development/contributing.rst#where-to-start)" section.

## Filing Issues

If you notice a bug in the code or documentation, or have suggestions for how we can improve either, feel free to create an issue on the [GitHub "issues" tab](https://github.com/pandas-dev/pandas/issues) using [GitHub's "issue" form](https://github.com/pandas-dev/pandas/issues/new). The form contains some questions that will help us best address your issue. For more information regarding how to file issues against *pandas*, please refer to the "[Bug reports and enhancement requests](https://github.com/pandas-dev/pandas/blob/master/doc/source/development/contributing.rst#bug-reports-and-enhancement-requests)" section.

## Contributing to the Codebase

This is an intermediate guide for making contributions, which is meant to summarize key points in the [Contributing Guide](https://pandas.pydata.org/docs/dev/development/contributing.html).  References are provided, but it may be beneficial to read the whole guide.

###	Clone Pandas to Make Changes

o	The code is hosted on GitHub, so you will need clone the repo to have a copy of the code
o	To get started, you will need to create your own fork to work on the code.  This should be done in a separate development environment separate from your existing Python environment.

Forking commands:
git clone https://github.com/your-user-name/pandas.git pandas-yourname
cd pandas-yourname
git remote add upstream https://github.com/pandas-dev/pandas.git

o	Please see the [Working with the code](https://pandas.pydata.org/docs/dev/development/contributing.html#working-with-the-code) section of the Contributing Guide for further information.

### Follow the Code Standards and Style Guide

o	Before submitting your changes, make sure the align to the pandas [Code Standards](https://pandas.pydata.org/docs/dev/development/contributing_codebase.html#code-standards)
o	Also, see the [Pandas Code Style Guide](https://pandas.pydata.org/docs/dev/development/code_style.html)

### Ensure your changes do not break any tests

o	You can find the testing information [Running the Test Suite](https://pandas.pydata.org/docs/dev/development/contributing_codebase.html#running-the-test-suite) section of the guide.
o	Tests can be run directly in your Git clone with the following command: 'pytest pandas'

### Push your changes

o	Once your changes are ready to be submitted, push your changes with the following command (replacing “shiny-new-feature” with "your branch name”: 'git push origin shiny-new-feature'
o	More information can be found in [Contributing your changes to pandas](https://pandas.pydata.org/docs/dev/development/contributing.html#contributing-your-changes-to-pandas) section

###	Create a pull request
o	This is done on the GitHub repo via the “Pull Request button”
o	If you have not done this before, please see the [Finally, make a pull request](https://pandas.pydata.org/docs/dev/development/contributing.html#finally-make-the-pull-request) section of the guide

### Review changes

o	We will review your changes, and you will most likely be asked to make additional changes before the final merge
o	Make any necessary changes to your contribution

###	Congratulations!

o	After the final review, we will merge your changes
o	You have now successfully contributed to the codebase
