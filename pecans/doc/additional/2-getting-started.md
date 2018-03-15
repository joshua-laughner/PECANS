# Getting started

## Quick start

Assuming you want to use one of the built in mechanisms:

1. (optional, but recommended) Set up a virtual environment (see below)
2. Clone/download the PECANS model
3. Install the required Python packages with `pip3 install -r requirements.txt`
4. Build a mechanism with the `build_pecans.py` script (run it to select options interactively).
5. (_not yet implemented_) Execute the model with the `run_pecans.py` script.

## Detail start
### 1) Set up a virtual environment (optional)

Although this is optional, we do recommend that you set up PECANS in a Python 3 virtual environment
if at all possible. This ensures that you can install the necessary packages without changing your
system wide Python installation.  

With a standard Python 3 installation, you create a virtual environment with the command `python3 -m venv MyVirtEnv`.
This will create a folder `MyVirtEnv` in the current folder (you can change that name to whatever 
you wish).

Users with an Anaconda python installation may wish to use `conda create --name MyVirtEnv` instead.

After creating the virtual environment, `cd` into that directory and type `source bin/activate` to
activate that virtual environment. You will need to perform this step each time you wish to
build or run the model, _if_ you build it in a virtual environment.

The first time you activate the virtual environment, it's a good idea to make sure `pip` (the Python
package installer) is up to date by running the command `pip3 install --upgrade pip`.


### 2) Clone or download the PECANS model

I recommend you clone the model. That way you'll be able to get any updates with a simple `git pull`
command, and if you make any changes to the code, you'll always be able to undo them. You will need
[git](https://git-scm.com/downloads) installed on your computer.  If you'd rather just download a 
compressed folder with the model, see the last paragraph.

To clone: open a terminal window and `cd` into your virtual environment folder (if you made one) or 
into a directory that you want the PECANS directory to exist in. Go to <https://github.com/firsttempora/PECANS> 
and click the "Clone or download" button. Choose "HTTPS" and copy the link. Back in your terminal, type
`git clone <link>`, pasting the link you just copied in place of `<link>` and press `enter`. Git will
create a PECANS folder and download the model into it.

To download a compressed folder, go to the [releases](https://github.com/firsttempora/PECANS/releases)
page of the GitHub repo and download the most recent release. Decompress it with your favorite extractor.


### 3) Install the required Python packages

PECANS has a number of packages it requires to run. All can be installed through the `pip3` utility.
In the top PECANS folder, there is a `requirements.txt` document. With your virtual environment activated
(if you made one), execute the command `pip3 install -r requirements.txt` from this folder. This tells
`pip3` to install each of the dependencies listed in that file.


### 4) Build a mechanism

The chemical mechanisms in PECANS are automatically coded into Cython from a mechanism file. This helps
speed up computation by allowing the bulk of the program to be compiled and optimized. The process is
automated by the `build_pecans.py` script in the second PECANS folder under the main one. With your
virtual environment activated (if you made one), execute this script with `./build_pecans.py` from the
folder that it resides in. This will present you with a list of available mechanisms. Simply choose one
and it will be build.  You can also specify a mechanism on the command line, e.q. `./build_peacns.py nox`
will build the "nox" mechanism.
