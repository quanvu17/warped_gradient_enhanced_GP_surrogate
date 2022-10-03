# Instructions to install TensorFlow v1 

These instructions are for Linux system. If you use a different operating system then the steps need some slight modifications.

Step 1: Install miniconda for your correct system Python version. Download it from here https://docs.conda.io/en/latest/miniconda.html and then run 

`bash -u Miniconda3-latest-Linux-x86_64.sh`

Step 2: Add miniconda3 to the .bashrc file (PATH). For example, if you installed miniconda in your home directory, add the following to your .bashrc file

`export PATH=$PATH:/home/<USERNAME>/miniconda3/bin`

Step 3: Initialise conda for shell interaction

`conda init`

Step 4: Create a new conda environment with all the required software. A full dependency list is given in the table here: https://www.tensorflow.org/install/source#gpu
                      
`conda create -n TFv1 -c conda-forge cudnn=7.6 cudatoolkit=10.0 gcc_linux-64=7.3 tensorflow-gpu=1.15`

Step 5: Install tensorflow in R

`install.packages("tensorflow")`

Step 6: Restart R and modify the script to load the installed conda environment

`library("reticulate")`

`use_condaenv("~/miniconda3/envs/TFv1/")`

Step 7: Run the scripts