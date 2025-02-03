Manual Usage
-------------
Alternatively, you can manually execute each step of the FELINE workflow. Follow the instructions below to process the data manually. Each command must be executed from the project's root directory.

1. **Copy the cube file**:
   Place the cube file into the `data/raw/` directory of the project.

2. **Create a virtual environment**:
   Create a Python virtual environment to manage dependencies:

   .. code-block:: bash

      python3.x -m venv venv

3. **Activate the environment**:
   Activate the virtual environment:

   .. code-block:: bash

      source venv/bin/activate

4. **Install dependencies**:
   Install the required Python packages:

   .. code-block:: bash

      pip install -r requirements.txt

5. **Preprocessing**:
   Apply a median filter to the data cube to remove continuum sources:

   .. code-block:: bash

      python src/preprocessing/lsdcat/median-filter-cube.py data/raw/<CUBENAME>.fits --signalHDU=1 --varHDU=2 --num_cpu=<num_cores> --width=151 --output=data/processed/med_filt.fits

6. **Apply spatial template matching**:
   Filter the data cube with a spatial template:

   .. code-block:: bash

      python src/preprocessing/lsdcat/lsd_cc_spatial.py --input=data/processed/med_filt.fits --SHDU=1 --NHDU=2 --threads=<num_cores> --gaussian --lambda0=7050 -pc 0.7 --classic --output=data/processed/spatial_cc.fits --overwrite

7. **Apply spectral template matching**:
   Filter the data cube with a spectral template:

   .. code-block:: bash

      python src/preprocessing/lsdcat/lsd_cc_spectral.py --input=data/processed/spatial_cc.fits --threads=<num_cores> --FWHM=250 --SHDU=1 --NHDU=2 --classic --output=data/processed/spectral_cc.fits --overwrite

8. **Construct signal-to-noise cube**:
   Build the S/N cube:

   .. code-block:: bash

      python src/preprocessing/lsdcat/s2n-cube.py --input=data/processed/spectral_cc.fits --output=data/processed/s2n_v250.fits --clobber --NHDU=2 --SHDU=1

9. **Transpose the cube**:
   For better cache access, transpose the cube:

   .. code-block:: bash

      python -m src.preprocessing.masking_and_transpose <CUBENAME>.fits s2n_v250.fits

10. **Compile and run the main program**:
    Compile and execute the FELINE binary with the following commands:

    .. code-block:: bash

      make
      ./feline.bin <ZLOW> <ZHIGH> <MAX_MATCH> <IGNORE_BELOW>

    .. note::
      | **For SDL Usage**: By default the Makefile activates SDL if it is available.
      | If you want to disable SDL you need to add the following macro to the make command:
      | ``make SDLavailable=0``

11. **Postprocessing**:
    Run scripts to detect objects, generate plots and create the PDF file:

    .. code-block:: bash

      python -m src.postprocessing.detect_objects s2n_v250.fits
      python -m src.postprocessing.create_final_plots <CUBENAME>.fits s2n_v250.fits sorted_catalog.txt med_filt.fits J0014m0028
      python -m src.postprocessing.create_pdf

12. **Find results**:
    After postprocessing, the results will be saved in the `data/pdf_files/` directory:

    .. code-block:: bash

      data/pdf_files/result_*.pdf

    PDF Format: ``result_YYYY_MM_DD_HH:MM:SS``

13. **Clean up**:
    To clean up temporary files after processing:

    .. code-block:: bash

      make clean