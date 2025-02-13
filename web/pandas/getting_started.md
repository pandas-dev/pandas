# Getting started

## Installation instructions

To install pandas, please reference the [installation page]({{ base_url}}docs/getting_started/install.html)
from the pandas documentation.

## Tutorials

You can learn more about pandas in the [tutorials]({{ base_url }}docs/getting_started/intro_tutorials/),
and more about JupyterLab in the
[JupyterLab documentation](https://jupyterlab.readthedocs.io/en/stable/user/interface.html).

## Books

The book we recommend to learn pandas is [Python for Data Analysis](https://amzn.to/3DyLaJc),
by [Wes McKinney](https://wesmckinney.com/), creator of pandas.

<a href="https://amzn.to/3DyLaJc">
    <img alt="Python for Data Analysis" src="{{ base_url }}static/img/books/pydata_book.gif"/>
</a>

## Videos

<iframe
  src="https://www.youtube.com/embed/_T8LGqJtuGc"
  style="width: 100%; max-width: 560px; height: 315px;"
  frameborder="0"
  allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
  allowfullscreen
></iframe>

## Cheat sheet

[pandas cheat sheet](https://pandas.pydata.org/Pandas_Cheat_Sheet.pdf)

## Try pandas in your browser (experimental)

You can try pandas in your browser with the following interactive shell
without needing to install anything on your system.

**Please note it can take a while (>30 seconds) before the shell is initialised and ready to run commands.**

**Running it requires a reasonable amount of bandwidth and resources, so it may not work properly on all devices or networks.**

<iframe
  src="./lite/repl/index.html?toolbar=1&kernel=python&execute=0&code=import%20pandas%20as%20pd%0Adf%20%3D%20pd.DataFrame%28%7B%22num_legs%22%3A%20%5B2%2C%204%5D%2C%20%22num_wings%22%3A%20%5B2%2C%200%5D%7D%2C%20index%3D%5B%22falcon%22%2C%20%22dog%22%5D%29%0Adf"
  style="width: 100%; max-width: 650px; height: 600px; border: 1px solid #130753;"
></iframe>
