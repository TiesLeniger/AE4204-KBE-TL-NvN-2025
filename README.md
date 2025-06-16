# Glidesign

Knowledge based engineering project repository of team 28 (Niels van Nieuwland and Ties Leniger).
This project is a glider design application written in Python with ParaPy.

Installation
------------
NOTE: You must have Matlab (version >2024a) installed on your computer
  1. Clone and open the repository
  2. Open a terminal and set the working directory to the repository
  3. `python -m venv .venv`
  4. `.venv/Scripts/activate`  (windows)
  5. `pip install -r requirements.txt`

Inputs
------
GliDesign reads top-level parameters from `glider.xlsx` in the input folder. Other inputs can be changed in the UI.

Running GliDesign
-----------------
  1. Make sure your virtual environment is activated (`.venv/Scripts/activate`)
  2. `python src/glidesign/main.py`
