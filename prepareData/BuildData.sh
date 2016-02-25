#!/bin/sh
#make executable chmod +x runHaplomat

python3 TransferAlphaToCode2dna.py
python3 BuildAllAllelesFrom_hla_nom_g.py
python3 BuildAllAllelesExpanded.py
python3 BuildP.py
python3 BuildH1.py
python3 BuildH1g.py  
python3 CreateH2.py
python3 AddGToH2.py
python3 PrintAllelesMissingIngCode.py




