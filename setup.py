import setuptools
from vfdbQuery.__init__ import __version__, __author__, __email__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="vfdbQuery",
    install_requires=['click', 'pandas', 'matplotlib', 'numpy'],
    python_requires='~=3.6',
    description="Package for quickly querying an assembly against the VFDB",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bfssi-forest-dussault/vfdbQuery",
    packages=setuptools.find_packages(),
    version=__version__,
    author=__author__,
    author_email=__email__,
    entry_points={
        'console_scripts': [
            'vfdbQuery=vfdbQuery.vfdbQuery:cli'
        ]
    }
)
