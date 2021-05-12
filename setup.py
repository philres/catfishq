"""
CatfishQ

"""
import os
from setuptools import setup, find_packages

os.environ['GIT_SSL_NO_VERIFY'] = 'true'

__version__ = '1.1.3'

setup(
    name='catfishq',
    version=__version__,
    author='philres',
    description='Cat FASTQ files',
    zip_safe=False,
    install_requires=[
        'pysam'
    ],
    packages=find_packages(exclude=("tests",)),
    entry_points={
        "console_scripts": [
            'catfishq = catfishq.cat_fastq:main',
        ]
    },
)
