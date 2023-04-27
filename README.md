# Mouse/Candida Peptide Overlaps

## Identification of intersection of tryptic peptides between Candida albicans and M. musculus

This repository contains Jupyter notebooks used to determine which peptides are non-discriminatory in co-MS experiments involving Mouse and Candida.

There are two notebooks available:
* **Candida_Mouse_Peptide_Intersections.ipynb** identifies tryptic peptides which are shared between proteins present in both Mouse and Candida. This _must_ be run first to generate the outputs for the second notebook.

* **MS_Output_Screening.ipynb** takes outputs from Spectronaut and identifies which of the peptides from an MS experiment are non-discriminatory.

This repository is available on MyBinder from where the notebooks can be run - just click the button below:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/bartongroup/Mouse-Candida-Peptide-Overlaps/HEAD)


### Standalone Setup

Installation of prerequisite packages is most easily carried out using [Conda](https://docs.conda.io/en/latest/).

1. Clone the git repository:

    ```git clone https://github.com/bartongroup/Mouse-Candida-Peptide-Overlaps.git```

2. Change into the directory containing the clone of the repository

    ```cd Mouse-Candida-Peptide-Overlaps```

3. Create a conda environment. The prerequisite packages are defined in `environment.yml`, and will be installed into an environment named `PepOverlap` 

    *MacOS Users* The EMBOSS package is not available from bioconda, so the environment creation will fail. Edit the environment.yml file to remove the 'EMBOSS=6.6.0' line prior to creating the environment. EMBOSS can instead be installed from homebrew with `brew install emboss`

    ```conda env create -f environment.yml```

4. Activate the conda environment:
    ```conda activate PepOverlap```

4. Start the Jupyter instance. A browser window will open displaying the JupyterLab interface

    ```jupyter lab```

5. Open the notebooks by double-clicking on their names in the left hand pane. The notebook can be run by selecting `Run All Cells` from the Jupyter `Run` menu

