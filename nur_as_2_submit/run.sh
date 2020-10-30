#!/bin/bash

echo "Run NUR assignment 2 by Matthijs van Groeningen"

#Scripts for question 1
echo "Run the script for 1a"
python3 ex1a.py

echo "Run the script for 1b (takes ~15 seconds)"
python3 ex1b.py

echo "Run the script for 1c"
python3 ex1c.py

echo "Run the script for 1d"
python3 ex1d.py

echo "Run the script for 1e"
python3 ex1e.py

echo "Run the script for 1f"
python3 ex1f.py

#Scripts for question 2
echo "Run the script for 2a"
python3 ex2a.py

echo "Run the script for 2b (takes ~40 seconds)"
python3 ex2b.py

echo "Run the script for 2c"
python3 ex2c.py

echo "Run the script for 2d"
python3 ex2d.py

echo "Generating the pdf"

pdflatex -shell-escape -interaction=nonstopmode -file-line-error A2.tex | grep -i ".*:[0-9]*:.*\|warning" 
bibtex A2
pdflatex -shell-escape -interaction=nonstopmode -file-line-error A2.tex | grep -i ".*:[0-9]*:.*\|warning"
pdflatex -shell-escape -interaction=nonstopmode -file-line-error A2.tex | grep -i ".*:[0-9]*:.*\|warning"

