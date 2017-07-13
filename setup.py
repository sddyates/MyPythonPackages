
# coding: utf-8 

from setuptools import setup

setup(
    name='sample',

    version='1.3.0',

    description='A package for analysing simulation results from
                 the MHD PLUTO.',

    # The project's main homepage.
    url='https://github.com/sddyates/plutopy',

    # Author details
    author='Simon Daley-Yates',
    author_email='sddyates@gmail.com',

    # License
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Data Analysis :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
    ],

    # What does your project relate to?
    keywords='Simulation analysis',

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['vtk', 'yt'],

)
