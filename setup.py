import os
from setuptools import setup,find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "workflowV2",
    version = "2.2219",
    author = "Patrick Neal",
    author_email = "prnmac12@gmail.com",
    description = ("An automation tool for QM calculations"),
    url = "https://github.com/neal-p/workflowV2",
    packages=find_packages(),
    long_description=read('README'),
)


