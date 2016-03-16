#
# Hapl-O-mat: A program for HLA haplotype frequency estimation
# 
# Copyright (C) 2016, DKMS gGmbH 
# 
# This file is part of Hapl-O-mat
# 
# Hapl-O-mat is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# Hapl-O-mat is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
 
# You should have received a copy of the GNU General Public License
# along with Hapl-O-mat; see the file COPYING.  If not, see
# <http://www.gnu.org/licenses/>.
# 


python3 TransferAlphaToCode2dna.py
python3 BuildAllAllelesFrom_hla_nom_g.py
python3 BuildAllAllelesExpanded.py
python3 BuildP.py
python3 BuildH1.py
python3 BuildH1g.py  
python3 CreateH2.py
python3 AddGToH2.py
python3 PrintAllelesMissingIngCode.py




