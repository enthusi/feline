Installation
============

FELINE requires specific software dependencies for installation and operation. Please follow the instructions below to set up the environment for FELINE:

Prerequisites
-------------
Ensure the following software is installed on your system:

- ``python 3.x`` (minimum version 3.8)
- ``python3.x-dev`` package
- ``python3.x-venv`` package
- A C/C++ compiler, such as ``gcc`` or ``clang``
- ``SDL2`` (Optional: Needed for graphical output during runtime)

Setup Instructions
------------------

1. **Clone the repository**:
   Clone the FELINE repository to your local machine:

   .. code-block:: bash

      git clone <repository-url>
      cd <repository-directory>

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

Optional GPU Acceleration
-------------------------
If you have an NVIDIA GPU and the required CUDA tools installed, you can enable GPU acceleration by running:

.. code-block:: bash

   make cuda

Testing the Installation
------------------------
To verify your installation, execute the workflow provided in the ``Makefile``:

1. Edit the ``CUBENAME`` and ``CUBELINK`` parameters in the ``Makefile``:

   - If the cube file is stored locally, place it in the project root directory and update ``CUBENAME`` accordingly.
   - Alternatively, provide the cube file URL and Name in ``CUBELINK`` and ``CUBENAME``.

2. Run the workflow:

   .. code-block:: bash

      make run

3. (Optional) For GPU acceleration, execute:

   .. code-block:: bash

      make cuda

Output
------
The final results will be available as PDF files in the ``data/pdf_files`` directory:

.. code-block:: bash

   data/pdf_files/result_*.pdf

Clean Up
--------
To remove temporary files and reset the project directory:

.. code-block:: bash

   make clean