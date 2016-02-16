#!/bin/sh
#make executable chmod +x runHaplomat

python BuildAllAllelesFrom_hla_nom_g.py
python BuildP.py
python BuildH1.py
python BuildH1g.py  
python BuildAllAllelesExpanded.py
python AddGToH2.py
python PrintAllelesMissingIngCode.py
#add printed alleles by hand
python CheckH2ForGG.py

mv H1.txt ../data
mv H1g.txt ../data
mv allAllelesExpanded.txt ../data
mv H2WithAddedG.txt ../data/H2.txt
mv P.txt ../data
cp code2dna.txt ../data 


