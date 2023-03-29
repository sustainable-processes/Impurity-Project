# Impurity Prediction

Contains developing code base for impurity prediction of reactions. At the moment, users need to have installed conda, anaconda or miniconda. Please refer to [Reaction Impurity Prediction using a Data Mining Approach](https://doi.org/10.1002/cmtd.202200062) for more details.

![image](https://user-images.githubusercontent.com/45038622/228471794-195ba7fa-7ec3-4454-b6f1-081c94a2f92a.png)


## Installation

Navigate to local project directory and:

`git clone https://github.com/sustainable-processes/Impurity-Project.git` 

OR

`git clone git@github.com:sustainable-processes/Impurity-Project.git` if SSH is configured

To install dependencies, type the following in anaconda prompt or a compatible IDE:

`conda env create -f environment_win.yml -n <your-env-name>` if OS is windows

OR

`conda env create -f environment_linux.yml -n <your-env-name>` if OS is linux/server based

Remember to activate the environment via `conda activate <your-env-name>`

## Viewing Results

Two tutorial files are given under the `Tutorial_Final` folder. `Tutorial-Jupyter Notebook_final.html` shows the workflow for the paracetamol case study (predicting impurities in paracetamol synthesis) with widgets and results visualized in an html file. It is advised to download this github repository as a .zip, and open `Tutorial-Jupyter Notebook_final.html` in Google Chrome for best results. `Tutorial-Jupyter Notebook_final.ipynb` is a jupyter notebook with the same workflow, but with no widget results shown (due to storage space limitations).
