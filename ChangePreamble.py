from os import walk
from shutil import move

#search relevant files
fileEndings = ['cc', 'h', 'py']
folders = ['src', 'include', 'prepareData', 'systemTest']

filesToChange = []
for folder in folders:
    for root, dirs, files in walk(folder):
        for file in files:
            fileSplitted = file.split('.')
            if len(fileSplitted) > 1:
                fileEnding = fileSplitted[1]

                if fileEnding in fileEndings:
                    filesToChange.append(root + '/' + file)

#new preamble
preambleC = '''/*
 * Hapl-o-Mat: A software for haplotype inference
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * Christian Schäfer
 * Kressbach 1
 * 72072 Tübingen, Germany
 *
 * T +49 7071 943-2063
 * F +49 7071 943-2090
 * cschaefer(at)dkms.de
 *
 * This file is part of Hapl-o-Mat
 *
 * Hapl-o-Mat is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Hapl-o-Mat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hapl-o-Mat; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */
'''

preamblePython = '''#
# Hapl-o-Mat: A software for haplotype inference
# 
# Copyright (C) 2016, DKMS gGmbH 
# 
# Christian Schäfer
# Kressbach 1
# 72072 Tübingen, Germany
#
# T +49 7071 943-2063
# F +49 7071 943-2090
# cschaefer(at)dkms.de
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
'''

#replace old preamble with new one
indicatorLine = ' <http://www.gnu.org/licenses/>.\n'
for filename in filesToChange:
    content = ''
    contentFound = False
    with open(filename, 'r+') as file:
        for line in file:
            if contentFound:
                content += line
            else:
                if indicatorLine in line:
                    file.readline()
                    contentFound = True


    with open(filename + '.tmp', 'w') as file:
        if filename.endswith('.cc') or filename.endswith('.h'):
            file.write(preambleC)
        else:
            file.write(preamblePython)
        file.write(content)

    move(filename + '.tmp', filename)

    print('Changed '+ filename + '.')
