<h1 id="contributing-to-pandas">Contributing to pandas</h1>
<h2 id="where-to-start">Where to start?</h2>
<p>All contributions, bug reports, bug fixes, documentation improvements, enhancements and ideas are welcome.</p>
<p>If you are simply looking to start working with the <em>pandas</em> codebase, navigate to the <a href="https://github.com/pydata/pandas/issues">GitHub &quot;issues&quot; tab</a> and start looking through interesting issues. There are a number of issues listed under <a href="https://github.com/pydata/pandas/issues?labels=Docs&amp;sort=updated&amp;state=open">Docs</a> and <a href="https://github.com/pydata/pandas/issues?labels=Good+as+first+PR&amp;sort=updated&amp;state=open">Good as first PR</a> where you could start out.</p>
<p>Or maybe through using <em>pandas</em> you have an idea of you own or are looking for something in the documentation and thinking 'this can be improved'...you can do something about it!</p>
<p>Feel free to ask questions on <a href="https://groups.google.com/forum/?fromgroups#!forum/pydata">mailing list</a></p>
<h2 id="bug-reportsenhancement-requests">Bug Reports/Enhancement Requests</h2>
<p>Bug reports are an important part of making <em>pandas</em> more stable. Having a complete bug report will allow others to reproduce the bug and provide insight into fixing. Since many versions of <em>pandas</em> are supported, knowing version information will also identify improvements made since previous versions. Often trying the bug-producing code out on the <em>master</em> branch is a worthwhile exercise to confirm the bug still exists. It is also worth searching existing bug reports and pull requests to see if the issue has already been reported and/or fixed.</p>
<p>Bug reports must:</p>
<ol>
<li><p>Include a short, self-contained Python snippet reproducing the problem. You can have the code formatted nicely by using <a href="http://github.github.com/github-flavored-markdown/">GitHub Flavored Markdown</a>: :</p>
<pre><code>```python
&gt;&gt;&gt; from pandas import DataFrame
&gt;&gt;&gt; df = DataFrame(...)
...
```</code></pre></li>
<li><p>Include the full version string of <em>pandas</em> and its dependencies. In recent (&gt;0.12) versions of <em>pandas</em> you can use a built in function: :</p>
<pre><code>&gt;&gt;&gt; from pandas.util.print_versions import show_versions
&gt;&gt;&gt; show_versions()</code></pre>
<p>and in 0.13.1 onwards: :</p>
<pre><code>&gt;&gt;&gt; pd.show_versions()</code></pre></li>
<li>Explain why the current behavior is wrong/not desired and what you expect instead.</li>
</ol>
<p>The issue will then show up to the <em>pandas</em> community and be open to comments/ideas from others.</p>
<h2 id="working-with-the-code">Working with the code</h2>
<p>Now that you have an issue you want to fix, enhancement to add, or documentation to improve, you need to learn how to work with GitHub and the <em>pandas</em> code base.</p>
<h3 id="version-control-git-and-github">Version Control, Git, and GitHub</h3>
<p>To the new user, working with Git is one of the more daunting aspects of contributing to <em>pandas</em>. It can very quickly become overwhelming, but sticking to the guidelines below will make the process straightforward and will work without much trouble. As always, if you are having difficulties please feel free to ask for help.</p>
<p>The code is hosted on <a href="https://www.github.com/pydata/pandas">GitHub</a>. To contribute you will need to sign up for a <a href="https://github.com/signup/free">free GitHub account</a>. We use <a href="http://git-scm.com/">Git</a> for version control to allow many people to work together on the project.</p>
<p>Some great resources for learning git:</p>
<ul>
<li>the <a href="http://help.github.com/">GitHub help pages</a>.</li>
<li>the <a href="http://docs.scipy.org/doc/numpy/dev/index.html">NumPy's documentation</a>.</li>
<li>Matthew Brett's <a href="http://matthew-brett.github.com/pydagogue/">Pydagogue</a>.</li>
</ul>
<h3 id="getting-started-with-git">Getting Started with Git</h3>
<p><a href="http://help.github.com/set-up-git-redirect">GitHub has instructions</a> for installing git, setting up your SSH key, and configuring git. All these steps need to be completed before working seamlessly with your local repository and GitHub.</p>
<h3 id="forking">Forking</h3>
<p>You will need your own fork to work on the code. Go to the <a href="https://github.com/pydata/pandas">pandas project page</a> and hit the <em>fork</em> button. You will want to clone your fork to your machine: :</p>
<pre><code>git clone git@github.com:your-user-name/pandas.git pandas-yourname
cd pandas-yourname
git remote add upstream git://github.com/pydata/pandas.git</code></pre>
<p>This creates the directory pandas-yourname and connects your repository to the upstream (main project) <em>pandas</em> repository.</p>
<p>You will also need to hook up Travis-CI to your GitHub repository so the suite is automatically run when a Pull Request is submitted. Instructions are <a href="http://about.travis-ci.org/docs/user/getting-started/">here</a>.</p>
<h3 id="creating-a-branch">Creating a Branch</h3>
<p>You want your master branch to reflect only production-ready code, so create a feature branch for making your changes. For example:</p>
<pre><code>git branch shiny-new-feature
git checkout shiny-new-feature</code></pre>
<p>The above can be simplified to:</p>
<pre><code>git checkout -b shiny-new-feature</code></pre>
<p>This changes your working directory to the shiny-new-feature branch. Keep any changes in this branch specific to one bug or feature so it is clear what the branch brings to <em>pandas</em>. You can have many shiny-new-features and switch in between them using the git checkout command.</p>
<h3 id="making-changes">Making changes</h3>
<p>Before making your code changes, it is often necessary to build the code that was just checked out. There are two primary methods of doing this.</p>
<ol>
<li><p>The best way to develop <em>pandas</em> is to build the C extensions in-place by running:</p>
<pre><code>python setup.py build_ext --inplace</code></pre>
<p>If you startup the Python interpreter in the <em>pandas</em> source directory you will call the built C extensions</p></li>
<li><p>Another very common option is to do a <code>develop</code> install of <em>pandas</em>:</p>
<pre><code>python setup.py develop</code></pre>
<p>This makes a symbolic link that tells the Python interpreter to import <em>pandas</em> from your development directory. Thus, you can always be using the development version on your system without being inside the clone directory.</p></li>
</ol>
<h2 id="contributing-to-the-documentation">Contributing to the documentation</h2>
<p>If you're not the developer type, contributing to the documentation is still of huge value. You don't even have to be an expert on <em>pandas</em> to do so! Something as simple as rewriting small passages for clarity as you reference the docs is a simple but effective way to contribute. The next person to read that passage will be in your debt!</p>
<p>Actually, there are sections of the docs that are worse off by being written by experts. If something in the docs doesn't make sense to you, updating the relevant section after you figure it out is a simple way to ensure it will help the next person.</p>
<h3 id="about-the-pandas-documentation">About the pandas documentation</h3>
<p>The documentation is written in <strong>reStructuredText</strong>, which is almost like writing in plain English, and built using <a href="http://sphinx.pocoo.org/">Sphinx</a>. The Sphinx Documentation has an excellent <a href="http://sphinx.pocoo.org/rest.html">introduction to reST</a>. Review the Sphinx docs to perform more complex changes to the documentation as well.</p>
<p>Some other important things to know about the docs:</p>
<ul>
<li><p>The <em>pandas</em> documentation consists of two parts: the docstrings in the code itself and the docs in this folder <code>pandas/doc/</code>.</p>
<p>The docstrings provide a clear explanation of the usage of the individual functions, while the documentation in this folder consists of tutorial-like overviews per topic together with some other information (what's new, installation, etc).</p></li>
<li>The docstrings follow the <strong>Numpy Docstring Standard</strong> which is used widely in the Scientific Python community. This standard specifies the format of the different sections of the docstring. See <a href="https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt">this document</a> for a detailed explanation, or look at some of the existing functions to extend it in a similar manner.</li>
<li><p>The tutorials make heavy use of the <a href="http://matplotlib.org/sampledoc/ipython_directive.html">ipython directive</a> sphinx extension. This directive lets you put code in the documentation which will be run during the doc build. For example:</p>
<pre><code>.. ipython:: python

    x = 2
    x**3</code></pre>
<p>will be rendered as</p>
<pre><code>In [1]: x = 2

In [2]: x**3
Out[2]: 8</code></pre>
<p>This means that almost all code examples in the docs are always run (and the output saved) during the doc build. This way, they will always be up to date, but it makes the doc building a bit more complex.</p></li>
</ul>
<h3 id="how-to-build-the-pandas-documentation">How to build the pandas documentation</h3>
<h4 id="requirements">Requirements</h4>
<p>To build the <em>pandas</em> docs there are some extra requirements: you will need to have <code>sphinx</code> and <code>ipython</code> installed. <a href="https://github.com/numpy/numpydoc">numpydoc</a> is used to parse the docstrings that follow the Numpy Docstring Standard (see above), but you don't need to install this because a local copy of <code>numpydoc</code> is included in the <em>pandas</em> source code.</p>
<p>Furthermore, it is recommended to have all <a href="http://pandas.pydata.org/pandas-docs/dev/install.html#optional-dependencies">optional dependencies</a> installed. This is not needed, but be aware that you will see some error messages. Because all the code in the documentation is executed during the doc build, the examples using this optional dependencies will generate errors. Run <code>pd.show_versions()</code> to get an overview of the installed version of all dependencies.</p>
<blockquote>
<p><strong>warning</strong></p>
<p>Sphinx version &gt;= 1.2.2 or the older 1.1.3 is required.</p>
</blockquote>
<h4 id="building-the-documentation">Building the documentation</h4>
<p>So how do you build the docs? Navigate to your local the folder <code>pandas/doc/</code> directory in the console and run:</p>
<pre><code>python make.py html</code></pre>
<p>And then you can find the html output in the folder <code>pandas/doc/build/html/</code>.</p>
<p>The first time it will take quite a while, because it has to run all the code examples in the documentation and build all generated docstring pages. In subsequent evocations, sphinx will try to only build the pages that have been modified.</p>
<p>If you want to do a full clean build, do:</p>
<pre><code>python make.py clean
python make.py build</code></pre>
<p>Starting with 0.13.1 you can tell <code>make.py</code> to compile only a single section of the docs, greatly reducing the turn-around time for checking your changes. You will be prompted to delete .rst files that aren't required, since the last committed version can always be restored from git.</p>
<pre><code>#omit autosummary and API section
python make.py clean
python make.py --no-api

# compile the docs with only a single
# section, that which is in indexing.rst
python make.py clean
python make.py --single indexing</code></pre>
<p>For comparison, a full documentation build may take 10 minutes. a <code>-no-api</code> build may take 3 minutes and a single section may take 15 seconds. However, subsequent builds only process portions you changed. Now, open the following file in a web browser to see the full documentation you just built:</p>
<pre><code>pandas/docs/build/html/index.html</code></pre>
<p>And you'll have the satisfaction of seeing your new and improved documentation!</p>
<h2 id="contributing-to-the-code-base">Contributing to the code base</h2>
<h3 id="code-standards">Code Standards</h3>
<p><em>pandas</em> uses the <a href="http://www.python.org/dev/peps/pep-0008/">PEP8</a> standard. There are several tools to ensure you abide by this standard.</p>
<p>We've written a tool to check that your commits are PEP8 great, <a href="https://github.com/hayd/pep8radius">pip install pep8radius</a>. Look at PEP8 fixes in your branch vs master with:</p>
<pre><code>pep8radius master --diff` and make these changes with `pep8radius master --diff --in-place`</code></pre>
<p>Alternatively, use <a href="http://pypi.python.org/pypi/flake8">flake8</a> tool for checking the style of your code. Additional standards are outlined on the <a href="https://github.com/pydata/pandas/wiki/Code-Style-and-Conventions">code style wiki page</a>.</p>
<p>Please try to maintain backward-compatibility. <em>Pandas</em> has lots of users with lots of existing code, so don't break it if at all possible. If you think breakage is required clearly state why as part of the Pull Request. Also, be careful when changing method signatures and add deprecation warnings where needed.</p>
<h3 id="test-driven-developmentwriting-code">Test-driven Development/Writing Code</h3>
<p><em>Pandas</em> is serious about <a href="http://en.wikipedia.org/wiki/Test-driven_development">Test-driven Development (TDD)</a>. This development process &quot;relies on the repetition of a very short development cycle: first the developer writes an (initially failing) automated test case that defines a desired improvement or new function, then produces the minimum amount of code to pass that test.&quot; So, before actually writing any code, you should write your tests. Often the test can be taken from the original GitHub issue. However, it is always worth considering additional use cases and writing corresponding tests.</p>
<p>Adding tests is one of the most common requests after code is pushed to <em>pandas</em>. It is worth getting in the habit of writing tests ahead of time so this is never an issue.</p>
<p>Like many packages, <em>pandas</em> uses the <a href="http://somethingaboutorange.com/mrl/projects/nose/">Nose testing system</a> and the convenient extensions in <a href="http://docs.scipy.org/doc/numpy/reference/routines.testing.html">numpy.testing</a>.</p>
<h4 id="writing-tests">Writing tests</h4>
<p>All tests should go into the <em>tests</em> subdirectory of the specific package. There are probably many examples already there and looking to these for inspiration is suggested. If you test requires working with files or network connectivity there is more information on the <a href="https://github.com/pydata/pandas/wiki/Testing">testing page</a> of the wiki.</p>
<p>The <code>pandas.util.testing</code> module has many special <code>assert</code> functions that make it easier to make statements about whether Series or DataFrame objects are equivalent. The easiest way to verify that your code is correct is to explicitly construct the result you expect, then compare the actual result to the expected correct result:</p>
<pre><code>def test_pivot(self):
    data = {
        &#39;index&#39; : [&#39;A&#39;, &#39;B&#39;, &#39;C&#39;, &#39;C&#39;, &#39;B&#39;, &#39;A&#39;],
        &#39;columns&#39; : [&#39;One&#39;, &#39;One&#39;, &#39;One&#39;, &#39;Two&#39;, &#39;Two&#39;, &#39;Two&#39;],
        &#39;values&#39; : [1., 2., 3., 3., 2., 1.]
    }

    frame = DataFrame(data)
    pivoted = frame.pivot(index=&#39;index&#39;, columns=&#39;columns&#39;, values=&#39;values&#39;)

    expected = DataFrame({
        &#39;One&#39; : {&#39;A&#39; : 1., &#39;B&#39; : 2., &#39;C&#39; : 3.},
        &#39;Two&#39; : {&#39;A&#39; : 1., &#39;B&#39; : 2., &#39;C&#39; : 3.}
    })

    assert_frame_equal(pivoted, expected)</code></pre>
<h4 id="running-the-test-suite">Running the test suite</h4>
<p>The tests can then be run directly inside your git clone (without having to install <em>pandas</em>) by typing::</p>
<pre><code>nosetests pandas</code></pre>
<p>The tests suite is exhaustive and takes around 20 minutes to run. Often it is worth running only a subset of tests first around your changes before running the entire suite. This is done using one of the following constructs:</p>
<pre><code>nosetests pandas/tests/[test-module].py
nosetests pandas/tests/[test-module].py:[TestClass]
nosetests pandas/tests/[test-module].py:[TestClass].[test_method]</code></pre>
<h4 id="running-the-performance-test-suite">Running the performance test suite</h4>
<p>Performance matters and it is worth considering that your code has not introduced performance regressions. Currently <em>pandas</em> uses the <a href="https://github.com/pydata/vbench">vbench library</a> to enable easy monitoring of the performance of critical <em>pandas</em> operations. These benchmarks are all found in the <code>pandas/vb_suite</code> directory. vbench currently only works on python2.</p>
<p>To install vbench:</p>
<pre><code>pip install git+https://github.com/pydata/vbench</code></pre>
<p>Vbench also requires sqlalchemy, gitpython, and psutil which can all be installed using pip. If you need to run a benchmark, change your directory to the <em>pandas</em> root and run:</p>
<pre><code>./test_perf.sh -b master -t HEAD</code></pre>
<p>This will checkout the master revision and run the suite on both master and your commit. Running the full test suite can take up to one hour and use up to 3GB of RAM. Usually it is sufficient to past a subset of the results in to the Pull Request to show that the committed changes do not cause unexpected performance regressions.</p>
<p>You can run specific benchmarks using the <em>-r</em> flag which takes a regular expression.</p>
<p>See the <a href="https://github.com/pydata/pandas/wiki/Performance-Testing">performance testing wiki</a> for information on how to write a benchmark.</p>
<h3 id="documenting-your-code">Documenting your code</h3>
<p>Changes should be reflected in the release notes located in doc/source/whatsnew/vx.y.z.txt. This file contains an ongoing change log for each release. Add an entry to this file to document your fix, enhancement or (unavoidable) breaking change. Make sure to include the GitHub issue number when adding your entry.</p>
<p>If your code is an enhancement, it is most likely necessary to add usage examples to the existing documentation. This can be done following the section regarding documentation.</p>
<h2 id="contributing-your-changes-to-pandas">Contributing your changes to <em>pandas</em></h2>
<h3 id="committing-your-code">Committing your code</h3>
<p>Keep style fixes to a separate commit to make your PR more readable.</p>
<p>Once you've made changes, you can see them by typing:</p>
<pre><code>git status</code></pre>
<p>If you've created a new file, it is not being tracked by git. Add it by typing :</p>
<pre><code>git add path/to/file-to-be-added.py</code></pre>
<p>Doing 'git status' again should give something like :</p>
<pre><code># On branch shiny-new-feature
#
#       modified:   /relative/path/to/file-you-added.py
#</code></pre>
<p>Finally, commit your changes to your local repository with an explanatory message. An informal commit message format is in effect for the project. Please try to adhere to it. Here are some common prefixes along with general guidelines for when to use them:</p>
<blockquote>
<ul>
<li>ENH: Enhancement, new functionality</li>
<li>BUG: Bug fix</li>
<li>DOC: Additions/updates to documentation</li>
<li>TST: Additions/updates to tests</li>
<li>BLD: Updates to the build process/scripts</li>
<li>PERF: Performance improvement</li>
<li>CLN: Code cleanup</li>
</ul>
</blockquote>
<p>The following defines how a commit message should be structured. Please reference the relevant GitHub issues in your commit message using GH1234 or #1234. Either style is fine, but the former is generally preferred:</p>
<blockquote>
<ul>
<li>a subject line with &lt; 80 chars.</li>
<li>One blank line.</li>
<li>Optionally, a commit message body.</li>
</ul>
</blockquote>
<p>Now you can commit your changes in your local repository:</p>
<pre><code>git commit -m</code></pre>
<p>If you have multiple commits, it is common to want to combine them into one commit, often referred to as &quot;squashing&quot; or &quot;rebasing&quot;. This is a common request by package maintainers when submitting a Pull Request as it maintains a more compact commit history. To rebase your commits:</p>
<pre><code>git rebase -i HEAD~#</code></pre>
<p>Where # is the number of commits you want to combine. Then you can pick the relevant commit message and discard others.</p>
<h3 id="pushing-your-changes">Pushing your changes</h3>
<p>When you want your changes to appear publicly on your GitHub page, push your forked feature branch's commits :</p>
<pre><code>git push origin shiny-new-feature</code></pre>
<p>Here origin is the default name given to your remote repository on GitHub. You can see the remote repositories :</p>
<pre><code>git remote -v</code></pre>
<p>If you added the upstream repository as described above you will see something like :</p>
<pre><code>origin  git@github.com:yourname/pandas.git (fetch)
origin  git@github.com:yourname/pandas.git (push)
upstream        git://github.com/pydata/pandas.git (fetch)
upstream        git://github.com/pydata/pandas.git (push)</code></pre>
<p>Now your code is on GitHub, but it is not yet a part of the <em>pandas</em> project. For that to happen, a Pull Request needs to be submitted on GitHub.</p>
<h3 id="review-your-code">Review your code</h3>
<p>When you're ready to ask for a code review, you will file a Pull Request. Before you do, again make sure you've followed all the guidelines outlined in this document regarding code style, tests, performance tests, and documentation. You should also double check your branch changes against the branch it was based off of:</p>
<ol>
<li>Navigate to your repository on GitHub--<a href="https://github.com/your-user-name/pandas">https://github.com/your-user-name/pandas</a>.</li>
<li>Click on Branches.</li>
<li>Click on the Compare button for your feature branch.</li>
<li>Select the base and compare branches, if necessary. This will be master and shiny-new-feature, respectively.</li>
</ol>
<h3 id="finally-make-the-pull-request">Finally, make the Pull Request</h3>
<p>If everything looks good you are ready to make a Pull Request. A Pull Request is how code from a local repository becomes available to the GitHub community and can be looked at and eventually merged into the master version. This Pull Request and its associated changes will eventually be committed to the master branch and available in the next release. To submit a Pull Request:</p>
<ol>
<li>Navigate to your repository on GitHub.</li>
<li>Click on the Pull Request button.</li>
<li>You can then click on Commits and Files Changed to make sure everything looks okay one last time.</li>
<li>Write a description of your changes in the Preview Discussion tab.</li>
<li>Click Send Pull Request.</li>
</ol>
<p>This request then appears to the repository maintainers, and they will review the code. If you need to make more changes, you can make them in your branch, push them to GitHub, and the pull request will be automatically updated. Pushing them to GitHub again is done by:</p>
<pre><code>git push -f origin shiny-new-feature</code></pre>
<p>This will automatically update your Pull Request with the latest code and restart the Travis-CI tests.</p>
<h3 id="delete-your-merged-branch-optional">Delete your merged branch (optional)</h3>
<p>Once your feature branch is accepted into upstream, you'll probably want to get rid of the branch. First, merge upstream master into your branch so git knows it is safe to delete your branch :</p>
<pre><code>git fetch upstream
git checkout master
git merge upstream/master</code></pre>
<p>Then you can just do:</p>
<pre><code>git branch -d shiny-new-feature</code></pre>
<p>Make sure you use a lower-case -d, or else git won't warn you if your feature branch has not actually been merged.</p>
<p>The branch will still exist on GitHub, so to delete it there do :</p>
<pre><code>git push origin --delete shiny-new-feature</code></pre>
