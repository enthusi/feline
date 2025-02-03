Automatic Usage
---------------
To use FELINE in an automated manner, you can use the pre-configured Makefile workflow.
Herein, a 3D data cube observed with VLT/MUSE is being fetched from the wonderful
public `AMUSED <https://amused.univ-lyon1.fr>`_ data base.

1. **Edit the Makefile**:
   Edit the `CUBELINK` and `CUBENAME` parameters inside the `Makefile`.
   If the cube file is stored locally, copy it into the /data/raw directory and update `CUBENAME`.
   If you wish to use a URL, update `CUBELINK`.

   .. code-block:: bash

       [27] CUBELINK := <link>
       [28] CUBENAME := <name>

2. **Run the workflow**:
   To run the workflow, use the following command:

   .. code-block:: bash

      make run

   .. note::
      | **For SDL Usage**: By default the Makefile activates SDL if it is available.
      | If you want to disable SDL you need to add the following macro to the make command:
      | ``make run SDLavailable=0``

3. **Optional GPU acceleration**:
   If a compatible NVIDIA GPU is available, and the necessary developer tools are installed, you can enable GPU acceleration:

   .. code-block:: bash

      make cuda

4. **Find results**:
   Once the workflow completes successfully and the PDF files are merged, you can find the results in the `data/pdf_files/` directory:

   .. code-block:: bash

      data/pdf_files/result_*.pdf

   PDF Format: ``result_YYYY_MM_DD_HH:MM:SS``

5. **Clean up**:
   To clean up temporary files:

   .. code-block:: bash

      make clean