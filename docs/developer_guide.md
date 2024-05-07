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

## Python conventions

### Style Guide

See [PEP 8 – Style Guide for Python Code](https://peps.python.org/pep-0008/) for details on the style guide.

### Function Signatures

To help with documentation, you're encouraged to use type hinting and full docstrings for functions and classes. Consider the following example:

```python
def kl_divergence(data1: np.ndarray, data2: np.ndarray) -> float:
    """Calculate the Kullback-Leibler Divergence.
    This provides a metric for comparing two 2D distributions

    Parameters
    ----------
        data1 (numpy.ndarray)                   - Array to compare against.
        data2 (numpy.ndarray)                   - Array to compare.

    Returns
    ----------
        kl  (float)                             - KL Divergence score.
    """
    kl = data1 * np.log(data1 / data2)
    kl = kl.sum()
    return kl
```


Here we have the function `kl_divergence` which takes two arguments `data1` and `data2` and returns a `float`. In Python, we can use "[type hinting](https://docs.python.org/3/library/typing.html)"  to help the user know what to expect from the function. The typing tells the user that the function expects two `np.ndarray` and will return a `float`.

When we write functions we can use the following format:

```python
def function_name(argument: type = default) -> return_type:
    ...
```


Where `argument` is the name of the argument. We specify the "type" of the argument with  `: type` after the argument. We can also specify a default value using `= default`.
Finally, we specify the return type using `-> return_type`. For example, if we have a function that takes an argument "apples" which is a `float` and an argument "normalize" which is a `bool` (which defaults to `False`) and returns an `int` we would write:

```python
def my_function( apples : float, normalize : bool = False) -> int:
    ...
```

Docstrings should also be used:

* To describe, in a line, what the function does, with an optional longer description.
* List the parameter, their type and what the parameters do within the function and default values.
* List the returns, their type and what the parameter represents.

Consider the following example:
```python
def get_pi(t_sleep : float, get_tau : bool = False) -> float:
    '''Get the value for pi or tau

    Function to sleep for t_sleep seconds and return either pi or tau.

    Parameters
    ----------
        t_sleep (float)                     - Number of seconds to sleep for
        get_tau (bool)                      - If true function will return 2*pi
                                              Defaults to False (return pi)


    Returns
    ----------
        pi (float)                          - Pi (or Tau if get_tau == True).
    '''

    sleep(t_sleep)

    pi = np.pi
    if get_tau:
        pi *= 2
    return pi
```
in the above, we specify the types of the variables that are getting passed using typing. In the docstring we have:

* a short message (`Get the value for pi or tau`) to give a short description of the function.
* a longer message (`Function to sleep for t_sleep seconds and return either pi or tau.`) to give a more detailed message about the function.
* the `Parameters` passed, what they're used for and what default values are.
* the `Returns` which tell what is being returned, their type and the expected behavior when various parameters are applied.

When a user types `help(get_pi)`, they will be given the message:
```
    '''Get the value for pi or tau

    Function to sleep for t_sleep seconds and return either pi or tau.

    Parameters
    ----------
        t_sleep (float)                     - Number of seconds to sleep for
        get_tau (bool)                      - If true function will return 2*pi
                                              Defaults to False (return pi)


    Returns
    ----------
        pi (float)                          - Pi (or Tau if get_tau == True).
    '''
```
Additionally this message will be used to generate the documentation for the function.
