#Hapl-o-Mat v 1.0

##General information
Hapl-o-Mat is software for haplotype inference via an
expectation-maximization algorithm. It supports processing and
resolving various forms of HLA genotype data.

##The latest version
The latest version can be found on the Github server under
[https://github.com/DKMS/Hapl-o-Mat](https://github.com/DKMS/Hapl-o-Mat)

##Requirements
Installing Hapl-o-Mat requires a C++ compiler with C++11 support. To
produce the required data on HLA nomenclature python is required.

##Installing under Linux
Get Hapl-o-Mat via 
```
git clone https://github.com/DKMS/Hapl-o-Mat
```
Using GNU compiler collection (GCC), Hapl-o-Mat is compiled via a
Makefile:
```
make
make clean
```
If you use another compiler than GCC, adapt the compiler flags in
the Makefile, e.g. -march=native.

##Installing under Windows
A detailed explanation on how to install Hapl-o-Mat under Windows using
the Eclipse IDE for C/C++ Developers can be found in
detailedGettingStartedWindows.pdf.

##Usage
For information on how to use Hapl-o-Mat and some tutorials using Linux
follow the guide detailedGettingStartedLinux.pdf or gettingStarted if you
are a seasoned Linux user. If you use Windows, follow the guide 
detailedGettingStartedWindows.pdf

##Citation
If you use Hapl-o-Mat for your research, please cite

(Until our journal article appears)
C. Schaefer, A.H. Schmidt, J. Sauter: Hapl-O-mat: A Versatile
Software for Haplotype Frequency Estimation. HLA (2016), 87, 236–320

##Contributors
If you want to participate in actively developing Hapl-o-Mat please
join via [Github](https://github.com/DKMS/Hapl-o-Mat)

##Author
Christian Schaefer                                                                                                                    
DKMS gGmbH                                                                                                                            
Kressbach 1                                                                                                                           
72072 Tuebingen, Germany                                                                                                              

T +49 7071 943-2063                                                                                                                   
F +49 7071 943-2090                                                                                                                   
cschaefer(at)dkms.de                                                                                                                  

##License
Copyright (C) 2016, DKMS gGmbH 

DKMS gGmbH                                                                                                                            
Kressbach 1                                                                                                                           
72072 Tuebingen, Germany                                                                                                              

T +49 7071 943-2063                                                                                                                   
F +49 7071 943-2090                                                                                                                   
cschaefer(at)dkms.de   
  
Hapl-o-Mat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3, or (at your option)
any later version.
 
Hapl-o-Mat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with Hapl-o-Mat; see the file COPYING.  If not, see
<http://www.gnu.org/licenses/>.
