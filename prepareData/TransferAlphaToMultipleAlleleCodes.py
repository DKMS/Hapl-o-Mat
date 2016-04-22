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


#Adapt format of alpha.v3.txt to MultipleAlleleCodes.txt.

with open('MultipleAlleleCodes.txt', 'w') as outFile:
    with open('alpha.v3.txt') as file:
        for line in file:
            line = line.replace('*', '')
            line = line.rstrip('\r\n')
            line = line.lstrip('\t')
            line += '\n'
            outFile.write(line)
        
