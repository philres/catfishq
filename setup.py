"""
CatfishQ

"""
import os
from setuptools import setup, find_packages

from catfishq import __version__

os.environ['GIT_SSL_NO_VERIFY'] = 'true'

setup(
    name='catfishq',
    version=__version__,
    author='philres',
    description='Cat FASTQ files',
    license='MIT',
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
