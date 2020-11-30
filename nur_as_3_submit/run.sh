#!/bin/bash

echo "Run NUR assignment 3 by Matthijs van Groeningen"

echo "Downloading satgal files..."
if [ ! -e satgals_m11.txt ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/files/satgals_m11.txt
  wget https://home.strw.leidenuniv.nl/~daalen/files/satgals_m12.txt
  wget https://home.strw.leidenuniv.nl/~daalen/files/satgals_m13.txt
  wget https://home.strw.leidenuniv.nl/~daalen/files/satgals_m14.txt
  wget https://home.strw.leidenuniv.nl/~daalen/files/satgals_m15.txt
fi


#Scripts for question 1
echo "Run the script for 1"
python3 ex1.py

#Scripts for question 2
echo "Run the script for 2"
python3 ex2.py


echo "Generating the pdf"

pdflatex A3.tex
bibtex A3.aux
pdflatex A3.tex
pdflatex A3.tex


