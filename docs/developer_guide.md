# Developer Guide

## Avoid pushing to main

Always make a development branch for developing code. 
Once happy with the code, and assuming it passes the required tests, make a pull request.


## Pre-commit

Pre-commit hooks are set up to perform linting and formatting. 
Before committing you should check that the code passes the pre-commit checks. 
If you haven't already, install [pre-commit](https://pre-commit.com/):
```
pip install pre-commit
```

To check files before committing run:
```
pre-commit run --all-files
```


## Conventional Commits

You are encouraged to use [commitizen](https://commitizen-tools.github.io/commitizen/) to generate [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/). 
To install commitizen:
```
pip install --user -U Commitizen
```

Modify some file(s) and then add the file to the staging area:
```
git add path/to/modified_file.py
```

Use commitizen to generate the commit comment:
```
cz commit
# or
cz c
```
Fill out the form. 

When you're happy with the changes to be committed:

```
git push
```