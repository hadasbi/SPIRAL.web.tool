
This repo contains the code to the SPIRAL website and software, as described in:

SPIRAL: Significant Process InfeRence ALgorithm for single cell RNA-sequencing and spatial transcriptomics\ Hadas Biran, Tamar Hashimshony, Tamar Lahav, Or Efrat, Yael Mandel-Gutfreund and Zohar Yakhini.

Full details are available at: https://spiral.technion.ac.il/

# How can you run SPIRAL on your data?

Choose one of the options below to run SPIRAL on your data. When done, you can upload the "spiral_results.zip" file to the website to view your results in the SPIRAL results panel, as explained at https://spiral.technion.ac.il/how_to_view.

## Option 1: download the packaged software 

A detailed explanation on this can be found at: https://spiral.technion.ac.il/how_to_run

## Option 2: clone the repository

Step-by-step explanation (for Windows and linux):

1. SPIRAL requires Python>=3.10, so verify your version of Python:
    ```sh 
    python --version
    ```
    or:
    ```sh 
    python3 --version
    ```
    If your version of Python is lower than 3.10, follow the instructions at [https://www.python.org/downloads/](https://www.python.org/downloads/) to upgrade it.

2. In the terminal, navigate to where you want the repository to be saved and clone it:
    ```sh 
    git clone https://github.com/hadasbi/SPIRAL.web.tool.git
    ```
3. Enter the repository folder:
    ```sh 
    cd SPIRAL.web.tool
    ```
4. Create a virtual environment:
    ```sh 
    python -m venv spiral_venv
    ```
    or:
    ```sh 
    python3 -m venv spiral_venv
    ```
5. Activate the virtual environment.
    
    In linux:
    ```sh 
    source ./spiral_venv/bin/activate
    ```
    In Windows:
    ```sh 
    spiral_venv\Scripts\activate
    ```
6. Install required Python packages:
    ```sh 
    pip3 install -r ./requirements.txt
    ```
7. Create a new directory called "analysis" inside the static directory:
    ```sh 
    cd static
    mkdir analysis
    cd ..
    ```
8. Run the script:
    ```sh 
    python SPIRAL.py
    ```
9. Follow the prompts to input necessary parameters and data.

10. When done, deactivate the virtual environment:
    ```sh 
    deactivate
    ```
_______________________________________________________________________________________________________________________________________
Please contact us by [email](mailto:spiral.web.tool@gmail.com) if you have any questions.
