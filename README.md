# Important note:
This repo contains code to a website that enables data upload by the user, then runs SPIRAL on out server and displays the results panel. We ultimately came to the conclusion that this is not the right solution. Instead, we developed a SPIRAL software that the user can download and run on any server. The result.zip file can then be uploaded to our website to view the results panel. 

Full details are available at: https://spiral.technion.ac.il/

________________________________________________________________________________________________________________________________________

# SPIRAL
SPIRAL is an algorithm that relies on a Gaussian statistical model to produce a comprehensive overview of significant
processes in single cell RNA-seq, spatial transcriptomics or bulk RNA-seq. SPIRAL detects biological processes by
identifying the subset of genes involved and the subset of cells, spots or samples.

SPIRAL is available as a web interface at https://spiral.technion.ac.il/.
You can also run SPIRAL locally by downloading the repository, and then run SPIRAL_local_run.py.

## How to run the SPIRAL algorithm on your data locally
To execute the SPIRAL algorithm on your data locally, bypassing the need to utilize the webtool, simply employ the SPIRAL_Local_run.py script. 
This convenient tool facilitates running the algorithm directly on your local data.

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
    python SPIRAL_Local_run.py
    ```
9. Follow the prompts to input necessary parameters and data.
10. When done, deactivate the virtual environment:
    ```sh 
    deactivate
    ```

Please contact us by [email](mailto:spiral.web.tool@gmail.com) if you have any questions.
