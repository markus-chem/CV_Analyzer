#!/usr/bin/env python
import json
from setuptools import setup


if __name__ == '__main__':
    """
    The first part compiles the broad package, the necessary
    information is given in the setup.json file.
    """
    with open('setup.json', 'r') as info:
        kwargs = json.load(info)
    setup(
        **kwargs
        )
