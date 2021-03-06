/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel
Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. PhyloBayes is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a
copy of the GNU General Public License along with PhyloBayes. If not, see
<http://www.gnu.org/licenses/>.

**********************/
#ifndef SOURCES_STATESPACE_H_
#define SOURCES_STATESPACE_H_
#include <string>

#include "BiologicalSequences.h"
// pure interface
//
class StateSpace {
 public:
  virtual ~StateSpace() {}

  virtual std::string GetState(int state) = 0;
  virtual int GetNstate() = 0;

  virtual int GetState(std::string from) = 0;
};

// simple state space: assumes that states are referred to using a one-letter
// code
//
class SimpleStateSpace : public StateSpace {
 public:
  int GetState(std::string from);

  int GetNstate() { return Nstate; }

  std::string GetState(int state);

  char GetCharState(int state) { return AlphabetSet[state]; }

 protected:
  int Nstate;
  char* Alphabet;
  int NAlphabetSet;
  char* AlphabetSet;
};

class DNAStateSpace : public SimpleStateSpace {
 public:
  DNAStateSpace();
  ~DNAStateSpace();
};

class RNAStateSpace : public SimpleStateSpace {
 public:
  RNAStateSpace();
  ~RNAStateSpace();
};

class ProteinStateSpace : public SimpleStateSpace {
 public:
  ProteinStateSpace();
  ~ProteinStateSpace();
};

#endif  // SOURCES_STATESPACE_H_
