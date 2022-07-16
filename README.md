# PSP
The proposed approach consists of 2 main components, a fragment generator and a protein predictor.
First run the fragment generator, then run the predictor.
If you already have the fragments, you don't need to run the fragment generator. Just set the PSP project's protein_loader.py file to the directory where the fragments are.


The protein-data-fetcher project is used to generate the protein fragments.

Use the frag_picker.sh file to configure the protocol used to generate the fragments.

Use the protein_fetcher_config.yaml file to configure the directory of secondary structure predictors used.

To run the project use:
"python3 protein_fetcher.py protein_name"
ex: python3 protein_fetcher.py 1crn

-----------------------------

The PSP project is used to make protein predictions.

The protein_data folder inside the project contains all the input needed to perform the prediction, including the fragments generated in the previous step. So before running this project, make sure this folder is up to date with all the information.

The protein_loader.py file is responsible for collecting this information. Here you can configure the path of the protein_data folder.

In the main.py file you configure the directory where the result will be stored.

In the file config.yaml you can configure some parameters of the algorithm, use the parameter proteinName to identify the protein that will be executed.

Note: the project runs one protein at a time.

To run this project use the command:
python3.6 main.py
