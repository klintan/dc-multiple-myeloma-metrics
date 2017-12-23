import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "dc-metrics",
    version = "0.0.1",
    author = "Andreas Klintberg",
    author_email = "ankl@kth.se",
    description = ("A range of medical research related metrics"),
    license = "MIT",
    keywords = "metrics accuracy",
    url = "https://github.com/klintan/dc-multiple-myeloma-metrics",
    packages=['metrics', 'tests'],
    long_description=read('README')
)