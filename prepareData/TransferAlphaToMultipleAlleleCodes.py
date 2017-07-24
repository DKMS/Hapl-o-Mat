#
# Hapl-o-Mat: A software for haplotype inference
# 
# Copyright (C) 2016, DKMS gGmbH 
# 
# Dr. Juergen Sauter
# Kressbach 1
# 72072 Tuebingen, Germany
#
# T +49 7071 943-2060
# F +49 7071 943-2090
# sauter(at)dkms.de
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
# 
# You should have received a copy of the GNU General Public License
# along with Hapl-o-Mat; see the file COPYING.  If not, see
# <http://www.gnu.org/licenses/>.
# 


#Adapt format of alpha.v3.txt to MultipleAlleleCodes.txt.

def transferAlphaToMultipleAlleleCodes():

    print('Transfer file alpha.v3.txt to file MultipleAlleleCodes.dat')
    with open('MultipleAlleleCodes.txt', 'w') as outFile:
        with open('alpha.v3.txt') as file:
            file.readline()
            file.readline()
            for line in file:
                if not line in ('\n', '\r\n'):
                    line = line.replace('*', '')
                    line = line.rstrip('\r\n')
                    line = line.lstrip('\t')
                    line += '\n'
                    outFile.write(line)



if __name__ == "__main__":

    transferAlphaToMultipleAlleleCodes()
        
