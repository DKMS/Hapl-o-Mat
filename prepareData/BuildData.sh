#
# Hapl-o-Mat: A program for HLA haplotype frequency estimation
# 
# Copyright (C) 2016, DKMS gGmbH 
# 
# This file is part of Hapl-o-Mat
# 
# Hapl-o-Mat is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# Hapl-o-Mat is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
 
# You should have received a copy of the GNU General Public License
# along with Hapl-o-Mat; see the file COPYING.  If not, see
# <http://www.gnu.org/licenses/>.
# 

python TransferAlphaToMultipleAlleleCodes.py
python BuildAllAllelesFrom_hla_nom_g.py
python BuildAllAllelesExpanded.py
python BuildP.py
python BuildLargeG.py
python BuildSmallg.py
python BuildAmbiguity.py
python AddGToAmbiguity.py
python AddAllelesMissingIngCode.py
mv AllAllelesExpanded.txt Ambiguity.txt LargeG.txt MultipleAlleleCodes.txt P.txt Smallg.txt ../data



