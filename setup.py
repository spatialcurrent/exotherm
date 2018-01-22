#!/usr/bin/env python


from setuptools import setup

setup(
    name='exotherm',
    version='1.0.0',
    install_requires=[],
    author='Spatial Current Developers',
    author_email='opensource@spatialcurrent.io',
    license='BSD License',
    url='https://github.com/spatialcurrent/exotherm/',
    keywords='python',
    description='A Python library for transforming common geospatial objects.',
    long_description=open('README.rst').read(),
    download_url="https://github.com/spatialcurrent/exotherm/zipball/master",
    packages=["exotherm"],
    package_data={'': ['LICENSE', 'NOTICE'], 'exotherm': ['*.pem']},
    package_dir={'exotherm': 'exotherm'},
    include_package_data=True,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ]
)
