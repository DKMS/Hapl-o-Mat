/*
 * Hapl-o-Mat: A software for haplotype inference
 *
 * Copyright (C) 2016, DKMS gGmbH 
 *
 * Dr. Jürgen Sauter
 * Kressbach 1
 * 72072 Tuebingen, Germany
 *
 * T +49 7071 943-2060
 * F +49 7071 943-2090
 * sauter(at)dkms.de
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

#ifndef Exceptions_header
#define Exceptions_header

#include <exception>
#include <string>

class MultipleAlleleCodeException : public std::exception{

 public:
  explicit MultipleAlleleCodeException(const std::string in_multipleAlleleCode)
    : exception(),
    errorMessage("Multiple allele code " + in_multipleAlleleCode + " not found.")
    {}

  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};

class MissingAlleleException : public std::exception{

 public:
  explicit MissingAlleleException(const std::string in_missingAllele,
				  const std::string in_fileName)
    : exception(),
    errorMessage("Allele " + in_missingAllele + " does not exist in " + in_fileName + ". Please consider the information in manageInput/checkInputDeprecatedAlleles/DeprecatedMultiAlleleCodes.txt.")
    {}

  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};

class AlleleResolutionException : public std::exception{

 public:
  explicit AlleleResolutionException(const std::string in_allele)
    : exception(),
    errorMessage("Resolution of " + in_allele + " not known.")
    {}

  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};

class ResolutionException : public std::exception{


 public:
  explicit ResolutionException(const std::string in_wantedResolution, const std::string in_locus)
    : exception(),
    errorMessage("Resolution " + in_wantedResolution + " specified for locus " + in_locus + " not known.")
    {}

  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};

class ParameterAssignmentException : public std::exception{

 public:
  explicit ParameterAssignmentException(const std::string in_line)
    : exception(),
    errorMessage("Wrong format for assignment of " + in_line + ".")
    {}

  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};

class ParameterNotFoundException : public std::exception{

 public:
  explicit ParameterNotFoundException(const std::string in_parameterName)
    : exception(),
    errorMessage("Parameter " + in_parameterName + " not found.")
    {}

  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};


class InputFormatException : public std::exception{

 public:
  explicit InputFormatException()
    : exception(),
    errorMessage("Wrong input format (MAC, GLSC, GLS, READ).")
    {}

  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};

class FileException : public std::exception{

 public:
  explicit FileException(const std::string in_filename)
    : exception(),
    errorMessage("Could not open " + in_filename + ".")
    {}

  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};


class NotMatchingLociException_MAC : public std::exception{

 public:
  explicit NotMatchingLociException_MAC(const std::string in_locus)
    : exception(),
    errorMessage("Specified locus " + in_locus + " not found.")
    {}
  
  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};

class NotMatchingLociException_GLSC : public std::exception{

 public:
  explicit NotMatchingLociException_GLSC(const std::string in_locus,
					const std::string in_id)
    : exception(),
    errorMessage("Specified locus " + in_locus + " not found in id " + in_id + ".")
    {}

  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};

class MissingGenotypeException : public std::exception{
  
  virtual const char* what() const throw()
  {
    return "Encountered missing genotype.";
  }
};

class MissingGlidException : public std::exception{

 public:
  explicit MissingGlidException(const size_t in_glid)
    : exception(),
    errorMessage("GLS-id " + std::to_string(in_glid) + " not found in glid-file.")
    {}

  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};

class SplittingGenotypeException : public std::exception{

 public:
  virtual const char* what() const throw()
  {
    return "Too broad splitting of genotype.";
  }
};

class InputLineException : public std::exception{

 public:
  virtual const char* what() const throw()
  {
    return "Wrong line length.";
  }
};

class UnsafeSIZE_T_Cast : public std::exception{

 public:
    explicit UnsafeSIZE_T_Cast(const std::string convertee)
    : exception(),
    errorMessage("Could not convert unsinged integer to size_t " +
                convertee )
    {}
  
  virtual const char* what() const throw()
  {
    return errorMessage.c_str();
  }

 private:
  std::string errorMessage;
};


#endif
