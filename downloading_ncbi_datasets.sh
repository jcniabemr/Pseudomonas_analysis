####When downloading datasets from ncbi to commandline first install ncbi_datasets 
#Crete environment and install ncbi_datasets  
conda create -n ncbi_datasets
#Activate environment 
conda activate ncbi_datasets
#Install ncbi_datasets 
conda install -c conda-forge ncbi-datasets-cli
#alternitively in one step create and install 
conda create -n ncbi_datasets ncbi-datasets-cli -y 
#Then move to the desired download location
#Here i download 2700 genome with the following command 
datasets download genome taxon 'Pseudomonas syringae group' 
#For individual genome use 
datasets download genome accession "GCA_example"
#Download only CDS for eg transcript counting
datasets download genome accession "GCA_example" --include cds
#If this will take a considerable amount of time you may want to first run a screen session and request resources as below  
screen -S download_ncbi_data
#then
srun --partition=medium --cpus-per-task=4 --mem=8G --pty bash 
