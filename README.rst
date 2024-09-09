Tissue-Simulation-Toolkit 2.0
==============================

Welcome to Tissue Simulation Toolkit (TST) 2.0, a library for
two-dimensional simulations of Glazier and Graner's Cellular Potts
model (Glazier and Graner, 1993).

TST 2.0 is an efficient C++ library for two-dimensional Cellular Potts Simulations. It is suitable for simulations with live visualization as well as batch simulations on clusters.

TST 2.0 provides many recent extensions to the CPM, including

* Efficient edgelist algorithm
* Infinite number of PDE layers (forward Euler)
* A reacion diffusion solver on the CPU and on CUDA
* Interaction of CPM cells and PDE (secretion, absorption)
* Chemotaxis
* Length and connectivity constraints
* Act-CPM model (Niculescu et al., PLOS Comput Biol 2015)
* Discrete fibrous extracellular matrix (molecular dynamics)

The current version of the TST includes example programs for the
following published simulations:

* Differential adhesion driven cell sorting (Glazier and Graner, 1993)
* Cell elongation dependent vasculogenesis (blood vessel growth) (Merks et al., 2006) 
* Contact-inhibition dependent vasculogenesis and angiogenesis (Merks and Glazier, 2005; Merks and Glazier, 2006; Merks et al, PLoS Comput Biol 2008)
* Hybrid cellular Potts and bead-spring modeling of cells in fibrous extracellular matrix (Tsingos and Bakker et al, Biophys J. 2023)


and visualization of:

* Cells, according to cell type or anything you wish
* Chemical fields, using color ramps and contour lines (level sets)

Documentation
-------------

Documentation is avaiable `here <https://www.mathbioleiden.nl/tst-docs/docs/html/index.html>`_. Note that this is work in process and not all functions are fully documented yet!

Downloading and installing
--------------------------

TST 2.0 is available from GitHub at https://github.com/ComputationalAnimalDevelopment/Tissue-Simulation-Toolkit. It can be built and run on Windows, macOS and Linux using the instructions below.

Windows
~~~~~~~

The easiest way to install and work with the TST on Windows is via the Windows Subsystem for Linux. This provides an Ubuntu Linux-like environment within Windows, within which you can install TST. Opening a WSL2 terminal and following the Linux instructions should get you there.

macOS
~~~~~~~

On macOS, you need to install the XCode development environment from Apple to get the required tools, including the command line tools. You will need to specifically select the command line tools in the installer.

To install the dependencies, we recommend installing [Homebrew](https://brew.sh). Once you have that installed, you can install QT5, libpng and zlib using (see note on Qt below)

.. code-block:: bash

    brew install qt@5 libpng zlib

You may have to edit `src/Tissue-Simulation-Toolkit.pri` for qmake to be able to find them. 

Note on Qt: If you have an existing Qt installation (e.g. the open source installation through qt.io)  do not install Qt again through homebrew. Instead, ensure that qmake is in the path or edit the Makefile such that the full path for qmake is given. 

Next, you can get the source by cloning the repository from GitHub. You can use the following commands in a Terminal:

.. code-block:: bash

    git clone --recursive -b TST2.0 git@github.com:ComputationalAnimalDevelopment/Tissue-Simulation-Toolkit.git

If you are on a Mac then you will have to modify the file :code:`lib/MultiCellDS/v1.0/v1.0.0/Makefile` to get the TST to compile. Find the line

.. code-block:: bash

    export COMPILE_CFLAGS := -O3 -s -mfpmath=both -m64 -std=c++11

And remove the :code:`-s -mfpmath=both` so that it reads

.. code-block:: bash

    export COMPILE_CFLAGS := -O3 -m64 -std=c++11

The TST can then be built using

.. code-block:: bash

    Tissue-Simulation-Toolkit$ make

See below for how to run a simple simulation to test if it's all working.

Linux
~~~~~

To compile the TST, C and C++ compilers are needed, as well as the usual helper tools like :code:`ar` and :code:`ranlib`, and :code:`make` for the build system. The TST also requires the zlib, libpng, OpenCL and QT5 libraries. On a recent Ubuntu or another Debian-based distribution (we tested Ubuntu 22.04), you can install the requirements using

.. code-block:: bash

    apt install gcc g++ binutils make zlib1g-dev libpng-dev ocl-icd-opencl-dev libqt5opengl5-dev

To get the source, clone the repository from GitHub:

.. code-block:: bash

    git clone --recursive -b TST2.0 git@github.com:rmerks/Tissue-Simulation-Toolkit.git

The TST can then be built using

.. code-block:: bash

    Tissue-Simulation-Toolkit$ make

See below for how to run a simple simulation to test if it's all working.

Test the Tissue Simulation Toolkit
----------------------------------

If compilation was successful, then the 'bin/' folder contains an executable called 'vessel'. This executable needs to be run from the `bin/` folder, and passed the location of a parameter file. You can run a test simulation like this:

.. code-block:: bash

    Tissue-Simulation-Toolkit$ cd bin
    Tissue-Simulation-Toolkit/bin$ ./vessel ../data/chemotaxis.par

Troubleshooting
---------------

MultiCellDS not found
~~~~~~~~~~~~~~~~~~~~~

If you get the error

.. code-block:: bash

    lib/MultiCellDS/v1.0/v1.0.0/libMCDS/xsde: No such file or directory

you probablly forgot to specify the '--recursive' keyword when cloning from github. You can solve this with

.. code-block:: bash

    cd Tissue-Simulation-Toolkit
    Tissue-Simulation-Toolkit$ git submodule init
    Tissue-Simulation-Toolkit$ git submodule update

Unkown FP unit
~~~~~~~~~~~~~~

If you get the error:

.. code-block:: bash

    error: unknown FP unit 'both'
    make[2]: *** [MultiCellDS.o] Error 1
    make[1]: *** [objects] Error 2
    make: *** [MCDS] Error 2

Find the line

.. code-block:: bash

    export COMPILE_CFLAGS := -O3 -s -mfpmath=both -m64 -std=c++11

And remove the :code:`-s -mfpmath=both` so that it reads

.. code-block:: bash

    export COMPILE_CFLAGS := -O3 -m64 -std=c++11

Contributing
-------------

Contributions are very welcome! You can contribute by forking this repository and creating a pull request.

Formatting
~~~~~~~~~~~~~~~
Code is easier to read if it formatted in a nice way. To format your code you can use :code:`clang-format` with the :code:`.clang_format` style file. 
To do this run 

.. code-block:: bash
    
    clang-format --style=file:./clang_format YOUR_FILE


We also created a script that automatically formats any file that you commit using :code:`git commit`. To activate this script you have to let git know where it is located. You can do this by running:

.. code-block:: bash
    git config --local core.hooksPath .githooks/


Who do I talk to?
-----------------

* Roeland Merks
