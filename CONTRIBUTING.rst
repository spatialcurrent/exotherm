Contributing to exotherm
==============

Contributor License Agreement
-----------------

Thank you for your interest in contributing.  You will first need to agree to the license.  Simply post the following paragraph, with "your name" and "your github account" substituted as needed, to the first issue `#1`_.  Otherwise, you can email us `here`_.

.. _#1: https://github.com/spatialcurrent/exotherm/issues/1
.. _here: mailto:opensource@spatialcurrent.io

I, **< YOUR NAME > (< YOUR GITHUB ACCOUNT >)**, agree to the license terms.  My contributions to this repo are granted to **Spatial Current, Inc.** under a perpetual, worldwide, non-exclusive, no-charge, royalty-free, irrevocable license and/or copyright is transferred.

Versioning
-----------------

This application is still in alpha and not under semantic versioning yet.

Authors
-----------------

See `AUTHORS`_ for a list of contributors.

.. _AUTHORS: https://github.com/spatialcurrent/exotherm/blob/master/AUTHORS

Documentation
-----------------

To build documentation first install sphinx via PyPI_ with:

.. _PyPI: https://pypi.python.org/pypi

.. code-block:: bash

    pip install sphinx sphinx_rtd_theme

Change directory to the **docs** directory then build with:

.. code-block:: bash

    ./make-docs.sh

Tests
-----------------

Run tests with:

.. code-block:: bash

    python -m unittest exotherm.test
