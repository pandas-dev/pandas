# Getting started

## Try it in your browser

You can try `pandas` in your browser with the following interactive shell
without installing anything on your computer.

*Note it can take up to 30 seconds before the shell finishes loading and is ready to run commands.*

<iframe
  src="./lite/repl/index.html?toolbar=1&kernel=python&code=import%20pandas%20as%20pd&code=df%20=%20pd.DataFrame(%7B'num_legs':%20%5B2,%204%5D,%20'num_wings':%20%5B2,%200%5D%7D,%20index=%5B'falcon',%20'dog'%5D)"
  width="100%"
  height="500px"
></iframe>

## Installation instructions

The next steps provides the easiest and recommended way to set up your
environment to use pandas. Other installation options can be found in
the [advanced installation page]({{ base_url}}/docs/getting_started/install.html).

1. Download [Anaconda](https://www.anaconda.com/distribution/) for your operating system and
   the latest Python version, run the installer, and follow the steps. Please note:

    - It is not needed (and discouraged) to install Anaconda as root or administrator.
    - When asked if you wish to initialize Anaconda3, answer yes.
    - Restart the terminal after completing the installation.

    Detailed instructions on how to install Anaconda can be found in the
    [Anaconda documentation](https://docs.anaconda.com/anaconda/install/).

2. In the Anaconda prompt (or terminal in Linux or MacOS), start JupyterLab:

    <img class="img-fluid" alt="" src="{{ base_url }}/static/img/install/anaconda_prompt.png"/>

3. In JupyterLab, create a new (Python 3) notebook:

    <img class="img-fluid" alt="" src="{{ base_url }}/static/img/install/jupyterlab_home.png"/>

4. In the first cell of the notebook, you can import pandas and check the version with:

    <img class="img-fluid" alt="" src="{{ base_url }}/static/img/install/pandas_import_and_version.png"/>

5. Now you are ready to use pandas, and you can write your code in the next cells.

## Tutorials

You can learn more about pandas in the [tutorials]({{ base_url }}/docs/getting_started/intro_tutorials/),
and more about JupyterLab in the
[JupyterLab documentation](https://jupyterlab.readthedocs.io/en/stable/user/interface.html).

## Books

The book we recommend to learn pandas is [Python for Data Analysis](https://amzn.to/2KI5JJw),
by [Wes McKinney](https://wesmckinney.com/), creator of pandas.

<a href="https://amzn.to/2KI5JJw">
    <img alt="Python for Data Analysis" src="{{ base_url }}/static/img/pydata_book.gif"/>
</a>

## Videos

<iframe width="560" height="315" frameborder="0"
src="https://www.youtube.com/embed/_T8LGqJtuGc"
allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
allowfullscreen></iframe>

## Cheat sheet

[pandas cheat sheet](https://pandas.pydata.org/Pandas_Cheat_Sheet.pdf)
