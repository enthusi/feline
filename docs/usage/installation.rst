Installation
============

FELINE requires specific software dependencies for installation and operation. Please follow the instructions below to set up the environment for FELINE:

Prerequisites
-------------
Ensure the following software is installed on your system:

- ``python 3.x`` (minimum version 3.8)
- ``python3.x-dev`` package
- ``python3.x-venv`` package
- A C/C++ compiler, such as ``clang`` (recommended due to a significant performance boost compared to gcc) or ``gcc``
- ``SDL2`` (Optional: Needed for graphical output during runtime)

Setup Instructions
------------------

.. note::
   | **macOS users**: If you use ``clang``, make sure native ``clang`` is installed (``brew install clang``) since Apple's Clang doesn't support ``openmp``. Secondly, you need ``libomp``, e.g., ``brew install libomp``.
   | For users who want to use ``gcc``, only need to adjust the following Makefile lines:
   |
   | ``[1] CC = gcc-<version>``
   | ``[2] CFLAGS = -O3 -ffast-math -fopenmp -g -std=c99``
   |
   | **Linux users (Debian/Ubuntu)**: If you use ``clang``, you only need to install ``libomp-dev``, e.g., ``apt install libomp-dev``.
   | For users who want to use ``gcc``, only need to adjust the following Makefile lines:
   |
   | ``[1] CC = gcc``
   | ``[2] CFLAGS = -O3 -ffast-math -fopenmp -g -std=c99``

1. **Clone the repository**:
   Clone the FELINE repository to your local machine:

   .. code-block:: bash

      git clone git@github.com:enthusi/feline.git
      cd feline

2. **Set up a virtual environment**:
   Create and activate a virtual environment for Python:

   .. code-block:: bash

      python3 -m venv venv
      source venv/bin/activate

3. **Install dependencies**:
   Install the required Python packages listed in ``requirements.txt``:

   .. code-block:: bash

      pip install -r requirements.txt

4. **Build the program**:
   Compile the main FELINE binary:

   .. code-block:: bash

      make

   .. note::
      | **For SDL Usage**: By default the Makefile activates SDL if it is available.
      | If you want to disable SDL you need to add the following macro to the make command:
      | ``make SDLavailable=0``