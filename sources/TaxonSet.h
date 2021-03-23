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

#ifndef SOURCES_TAXONSET_H_
#define SOURCES_TAXONSET_H_

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

class Tree;
class Link;

class TaxonSet {
 public:
  TaxonSet(const std::string* names, int ntaxa);
  explicit TaxonSet(const Tree* tree, const Link* subgroup = 0);
  ~TaxonSet();

  int GetNtaxa() const;
  std::string GetTaxon(int index) const;
  int GetTaxonIndex(std::string intaxon) const;
  int GetTaxonIndexWithIncompleteName(std::string intaxon) const;
  int GetTaxonIndexWithCompleteName(std::string intaxon) const;
  void ToStream(std::ostream& os);

 private:
  void RecursiveGetSubSet(const Link* from, int& index);

  int Ntaxa;
  mutable std::map<std::string, int> taxmap;
  std::string* taxlist;
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

inline int TaxonSet::GetNtaxa() const { return Ntaxa; }
inline std::string TaxonSet::GetTaxon(int index) const {
  return taxlist[index];
}
inline int TaxonSet::GetTaxonIndex(std::string intaxon) const {
  return taxmap[intaxon] - 1;
}

#endif  // SOURCES_TAXONSET_H_
