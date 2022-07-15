# PSP
O projeto completo é composto por 2 projetos, um gerador de fragmentos e um preditor de proteínas.
Primeiro execute o gerador de fragmentos, em seguida execute o preditor.
Caso você já tenha os fragmentos, não é necessário executar o gerador de fragmentos. Basta configurar o arquivo protein_loader.py do projeto PSP para o diretório onde os fragmentos se encontram.


O projeto protein-data-fetcher é utilizado para gerar os fragmentos de proteínas.

Utilize o arquivo frag_picker.sh para configurar o protocolo utilizado para gerar os fragmentos.

Utilize o arquivo protein_fetcher_config.yaml para configurar o diretório dos preditores de estrutura secundária utilizados.

Para executar o projeto utilize:
"python3 protein_fetcher.py nome_proteina"
ex: python3 protein_fetcher.py 1crn

-----------------------------

O projeto PSP é utilizado para fazer as previsões das proteínas.

A pasta protein_data dentro do projeto contem toda a entrada necessária para realizar a predição, inclusive os fragmentos gerados na etapa anterior. Então antes de executar este projeto, certifique-se de que esta pasta está atualizada com todas as informações.

O arquivo protein_loader.py é responsável por colher estas informações. Aqui você pode configurar o caminho da pasta protein_data.

No arquivo main.py você configura o diretório onde o resultado será armazenado.

No arquivo config.yaml você pode configurar alguns parâmetros do algoritmo, utilize o parâmetro proteinName para identificar a proteína que será executada.

Obs: o projeto executa uma proteína por vez.

Para executar este projeto utilize o comando:
python3.6 main.py

-------------------------------

# PSP
The complete project consists of 2 projects, a fragment generator and a protein predictor.
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
