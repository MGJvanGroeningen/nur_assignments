#!/bin/bash

echo "Run NUR assignment 1 by Matthijs van Groeningen"

echo "Download cooling tables..."
if [ ! -d CoolingTables ]; then
  wget https://www.strw.leidenuniv.nl/WSS08/coolingtables_highres.tar.gz
fi

echo "Extract the cooling tables..."
if [ ! -d CoolingTables ]; then
  tar -xvf coolingtables_highres.tar.gz CoolingTables
fi

echo "Download wgs.dat and wss.dat..."
if [ ! -e wgs.dat ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/wgs.dat
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/wss.dat
fi

echo "Creating the plotting directory if it does not exist"
if [ ! -d "plots1b" ]; then
  echo "Directory does not exist create it!"
  mkdir plots1b
fi

#Scripts for question 1
echo "Run the script for 1a"
python3 ex1a.py

echo "Run the script for 1b (creates 101 plots, might take a few seconds)"
python3 ex1b.py

#Scripts for question 2
echo "Run the script for 2a"
python3 ex2a.py

echo "Run the script for 2b"
python3 ex2b.py

#Scripts for question 3
echo "Run the script for 3a"
python3 ex3a.py

echo "Run the script for 3b"
python3 ex3b.py

echo "Run the script for 3c"
python3 ex3c.py

# code that makes a movie of the movie frames
echo "Make the movie for question 1"
ffmpeg -framerate 10 -pattern_type glob -i "plots1b/snap*.png" -s:v 640x480 -c:v libx264 -profile:v high -level 4.0 -crf 10 -tune animation -preset slow -pix_fmt yuv420p -r 25 -threads 0 -f mp4 coolratesmovie.mp4

echo "Generating the pdf"

pdflatex A1.tex
bibtex A1.aux
pdflatex A1.tex
pdflatex A1.tex


