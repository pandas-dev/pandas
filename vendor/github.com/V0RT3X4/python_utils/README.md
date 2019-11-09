# python_utils [![CircleCI](https://circleci.com/gh/V0RT3X4/python_utils.svg?style=svg&circle-token=30fa8fb22fa45a521a5d728e9accde63c242c2b4)](https://circleci.com/gh/V0RT3X4/python_utils)
Python utilities and helper functions/classes/modules

## Sub Packages

- [AWS](#aws)
- [Database](#database)
- [Deployment](#deployment)

## Installation

Installation is done by using [submodule vendoring](#vendoring).
Vendor the package into your project as [below](#vendoring) then you can install
with
```
pip install vendor/github.com/V0RT3X4/python_utils/<sub_package>
```
or
```
echo vendor/github.com/V0RT3X4/python_utils/<sub_package> >> requirements.txt
pip install -r requirements.txt
```

## Aws

Helper modules for `s3` client side encryption. `ses` email processing
(s3 as an inbox). `lambda` function handeler types.

## Database

Data base connection helpers to get you a
[`SQLAlchemy`](https://www.sqlalchemy.org/) connection [`Engine`](https://docs.sqlalchemy.org/en/latest/core/engines_connections.html)
to an RDS or RedShift database using
`aws secretsmanager` for managing connection credentials and rotation, and with
SSL encryption.

## Deployment

Custom Deployment Jazz

## Installation - Vendoring the subtree
To install the scripts into your project it is recommended to vendor this module as a `git subtree` as opposed to a `git submodule`. You will have a version of this code in your repo, and you can easily update and push changes back upstream.

To make your life easier install [git-vendor](https://github.com/brettlangdon/git-vendor)

Then you can vendor the module into your repo and run installation scripts:
```
git vendor add python_utils git@github.com:V0RT3X4/python_utils.git master
```

finally you can install the modules you want
```
pip install vendor/github.com/V0RT3X4/python_utils/<sub_package>
```

to update the reference
```
git vendor update python_utils master
```

## AS Submodule

In the project directory
```
git submodule add \
    --name github.com/V0RT3X4/python_utils \
    git@github.com:V0RT3X4/python_utils.git \
    vendor/github.com/V0RT3X4/python_utils
```

Subsequently when you check out the source code (say in
[circleCI](https://circleci.com) or locally).
```
git clone git@github.com:/V0RT3X4/<your_project_name>.git
cd <your_project_name>
git submodule init
git submodule update --remote
```

finally you can install the modules you want
```
pip install vendor/github.com/V0RT3X4/python_utils/<sub_package>
```

## Contributing
To contribute and push changes back upstream add this repo as a remote.
```
git remote add -f python_utils git@github.com:V0RT3X4/python_utils.git
```
Push changes in the sub tree
```
git subtree push --prefix=vendor/github.com/V0RT3X4/python_utils python_utils some_branch
```

## [git-vendor](https://github.com/brettlangdon/git-vendor) installation

```
cd $(mktemp -d) && \
git clone https://github.com/brettlangdon/git-vendor &> /dev/null && \
cd git-vendor && \
sudo make install
```

or

```
brew install git-vendor
```
