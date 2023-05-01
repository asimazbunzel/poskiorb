===========
Development
===========

.. warning::

   These page is incomplete: it is, still, being under active development

This section gives more information for the development of the module and the steps to follow to
contribute to it

Installation
------------

To start contributing, install it with all the development tools needed

1. Clone the repository onto your machine:

::

    git clone git@github.com:asimazbunzel/poskiorb.git
    
or

::

    git clone https://github.com/asimazbunzel/poskiorb.git

2. Create a conda environment from the ``environment.yml`` file. It will install necessary modules
   for development such as ``black``, ``pyupgrade``, etc:

::
  
    conda activate poskiorb-dev

3. Activate the environment

::

    conda activate poskiorb-dev

4. Install ``poskiorb`` by running

::

    pip install .
