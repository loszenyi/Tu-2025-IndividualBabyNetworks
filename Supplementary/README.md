## Parcellation comparison notebooks

Two Jupyter notebooks were involved in the application of the CNC toolbox:

- `parcellation_comparison.ipynb`
- `parcellation_comparison_customizedAtlases.ipynb`

The main difference between them is whether customized atlases beyond those included in the CNC toolbox are used.

The standard notebook, `parcellation_comparison.ipynb`, uses the default atlas resources available through the CNC toolbox and does not require any additional setup.

The customized notebook, `parcellation_comparison_customizedAtlases.ipynb`, uses additional customized atlas files. To run this notebook, you will also need the files provided in the `customized_atlas_files` folder. These files should be added to the corresponding folders under the `data` path inside the CNC Python library.

Specifically, make sure to add the customized atlas files, configuration files, network name files, and updated list files to the appropriate locations in the CNC Python library.

The CNC toolbox is available here:
https://github.com/rubykong/cbig_network_correspondence