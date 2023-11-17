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

#ifndef SOURCES_TREE_H_
#define SOURCES_TREE_H_

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TaxonSet.h"

class Node {
 private:
  int index;
  std::string name;

 public:
  Node() : index(0), name("") {}
  explicit Node(std::string s) : index(0), name(s) {}
  explicit Node(const Node* from) : index(from->index), name(from->name) {}

  virtual ~Node() {}

  virtual std::string GetName() const { return name; }
  virtual void SetName(std::string inname) { name = inname; }
  int GetIndex() const { return index; }
  void SetIndex(int i) { index = i; }
};

class Branch {
 private:
  int index;
  std::string name;

 public:
  Branch() : index(0), name("") {}
  explicit Branch(std::string s) : index(0), name(s) {}
  explicit Branch(const Branch* from) : index(from->index), name(from->name) {}

  virtual ~Branch() {}

  virtual std::string GetName() const { return name; }
  virtual void SetName(std::string inname) { name = inname; }
  int GetIndex() const { return index; }
  void SetIndex(int i) { index = i; }
};

class Link {
 private:
  Link* next;
  Link* out;
  Branch* branch;
  Node* node;
  int index;

 public:
  double* tbl;

  Link() {
    tbl = 0;
    next = out = this;
    branch = 0;
    node = 0;
  }

  explicit Link(const Link* from) {
    tbl = 0;
    next = out = this;
    node = 0;
    branch = 0;
  }

  Link* Next() const { return next; }
  Link* Out() const { return out; }
  Branch* GetBranch() const { return branch; }
  Node* GetNode() const { return node; }

  void SetBranch(Branch* inbranch) {
    /*
    if (branch == inbranch)	{
            cout << "error in Link::SetBranch: branch and new branch are the
    same\n"; exit(1);
    }
    delete branch;
    */
    branch = inbranch;
  }
  void SetNode(Node* innode) {
    /*
    if (node == innode)	{
            cout << "error in Link::SetNode: node and new node are the same\n";
            exit(1);
    }
    delete node;
    */
    node = innode;
  }

  std::string GetNodeLeafSet() const {
    std::string s;
    if (isLeaf()) {
      s = GetNode()->GetName() + '|';
    } else {
      for (const Link* link = Next(); link != this; link = link->Next()) {
        s = s + link->Out()->GetNodeLeafSet();
      }
    }
    return s;
  }

  void SetIndex(int i) { index = i; }

  int GetIndex() const { return index; }

  void SetOut(Link* inout) { out = inout; }

  void SetNext(Link* innext) { next = innext; }

  void AppendTo(Link* link) {
    if (link) {
      link->next = this;
    }
  }

  void Insert(Link* link) {  // insert link after this
    link->next = next;
    next = link;
  }

  void InsertOut(Link* link) {  // insert link as out
    link->out = this;
    out = link;
  }

  // Move the next subtree after this link to be the previous of this link.
  // Should be call with 0;
  void Knit() {
    Link* previous = 0;
    for (previous = this; previous->next != this; previous = previous->Next()) {
    }
    previous->next = this->next;
    this->next = next->Next();
    previous->next->next = this;

    /*if(!link){
            link = next;
            next=next->next;
            link->next = this;
    }
    if(next == link->next){
            next = link;
    }
    else{
            next->Knit(link);
    }*/
  }

  bool isLeaf() const { return (next == this); }

  bool isUnary() const { return (next->Next() == this && !isLeaf()); }

  bool isRoot() const { return (out == this); }

  // degree : number of branches connecting to the node associated to this link
  int GetDegree() const {
    int d = 1;
    const Link* link = next;
    while (link != this) {
      d++;
      link = link->next;
    }
    return d;
  }

  const Link* GetUp(int& d) const {
    d = 1;
    const Link* link = link->Out();
    while (link->GetDegree() == 2) {
      link = link->Next()->Out();
      d++;
    }
    return link;
  }
};

class NewickTree {
 public:
  virtual ~NewickTree() {}
  virtual Link* GetRoot() = 0;
  virtual const Link* GetRoot() const = 0;

  void ToStream(std::ostream& os) const;
  void ToStream(std::ostream& os, const Link* from) const;
  double ToStreamSimplified(std::ostream& os, const Link* from) const;

  std::string GetLeftMost(const Link* from) const {
    if (from->isLeaf()) {
      return GetNodeName(from);
    }
    return GetLeftMost(from->Next()->Out());
  }

  std::string GetRightMost(const Link* from) const {
    if (from->isLeaf()) {
      return GetNodeName(from);
    }
    const Link* link = from->Next();
    while (link->Next() != from) {
      link = link->Next();
    }
    return GetRightMost(link->Out());
  }

  static void Simplify() { simplify = true; }

  virtual std::string GetNodeName(const Link* link) const = 0;
  virtual std::string GetBranchName(const Link* link) const = 0;

 protected:
  static bool simplify;
};

class Tree : public NewickTree {
 public:
  Tree();
  // default constructor: set member pointers to 0

  // void tree;
  explicit Tree(const TaxonSet* intaxset);

  void MakeRandomTree();

  explicit Tree(const Tree* from);
  // clones the entire Link structure
  // but does NOT clone the Nodes and Branches
  // calls RecursiveClone

  explicit Tree(std::string filename);
  // create a tree by reading into a file (netwick format expected)
  // calls ReadFromStream

  explicit Tree(std::istream& is);
  // create a tree by reading into a stream (netwick format expected)
  // calls ReadFromStream

  void ReadFromStream(std::istream& is);
  // reading a tree from a stream:
  // recursively invokes the two following functions

  virtual ~Tree();
  // calls RecursiveDelete
  // but does NOT delete the Nodes and Branches

  // Delete the leaf pointing by the next link and set everithing right.
  void DeleteNextLeaf(Link* previous);

  Link* Detach(Link* down, Link* up);

  void Attach(Link* down, Link* up, Link* todown, Link* toup);

  void NNIturn(Link* from);

  void RootAtRandom();
  void RootAt(Link* newroot);
  Link* ChooseLinkAtRandom();

  int CountInternalNodes(const Link* from);
  Link* ChooseInternalNode(Link* from, Link*& chosenup, int& n);
  int CountNodes(const Link* from);
  Link* ChooseNode(Link* from, Link*& chosenup, int& n);

  int DrawSubTree(Link*& down, Link*& up);
  void GrepNode(Link* from, Link*& down, Link*& up, int choose);

  bool RecursiveCheckDegree(const Link* from, int test = 3);
  bool CheckRootDegree(int test = 3);

  // Delete the unary Node wich from is paart of and set everithing right.
  void DeleteUnaryNode(Link* from);

  Link* GetRoot() { return root; }
  const Link* GetRoot() const { return root; }
  const TaxonSet* GetTaxonSet() const { return taxset; }

  void MakeTaxonSet() { taxset = new TaxonSet(this); }

  void RegisterWith(const TaxonSet* taxset, int id = 0);
  // Registers all leaves of the tree with an external TaxonSet
  // the taxon set defines a map between taxon names and indices (between 0 and
  // P-1) the tree is recursively traversed each leaf's name is looked for in
  // the map of the taxon set if not found : an error is emitted otherwise, the
  // leaf's index is set equal to the index of the corresponding taxon

  bool RegisterWith(const TaxonSet* taxset, Link* from, int& tot);
  // recursive function called by RegisterWith

  std::string GetBranchName(const Link* link) const {
    return link->GetBranch()->GetName();
  }

  std::string GetNodeName(const Link* link) const {
    return link->GetNode()->GetName();
    /*
    if (! link->isLeaf())	{
            return link->GetNode()->GetName();
    }
    std::string s = link->GetNode()->GetName();
    unsigned int l = s.length();
    unsigned int i = 0;
    while ((i < l) && (s[i] != '_')) i++;
    if (i == l)	{
            std::cerr << "error in get name\n";
            exit(1);
    }
    i++;
    return s.substr(i,l-i);
    */
  }
  // trivial accessors
  // they can be useful to override, so as to bypass Branch::GetName() and
  // Node::GetName()

  void EraseInternalNodeName();
  void EraseInternalNodeName(Link* from);

  // void Print(ostream& os,const Link* from) const ;
  // void Print(ostream& os) const;
  // printing int netwick format

  unsigned int GetSize() { return GetSize(GetRoot()); }

  int GetSize(const Link* from) const {
    if (from->isLeaf()) {
      return 1;
    } else {
      int total = 0;
      for (const Link* link = from->Next(); link != from; link = link->Next()) {
        total += GetSize(link->Out());
      }
      return total;
    }
    return 0;
  }

  int GetFullSize(const Link* from) const {
    if (from->isLeaf()) {
      return 1;
    } else {
      int total = 1;
      for (const Link* link = from->Next(); link != from; link = link->Next()) {
        total += GetFullSize(link->Out());
      }
      return total;
    }
    return 0;
  }

  Link* GetLCA(std::string tax1, std::string tax2) {
    bool found1 = false;
    bool found2 = false;
    Link* link = RecursiveGetLCA(GetRoot(), tax1, tax2, found1, found2);
    // std::cerr << tax1 << '\t' << tax2 << '\n';
    // Print(std::cerr,link);
    // std::cerr << '\n' << '\n';
    return link;
  }

  void Subdivide(Link* from, int Ninterpol);

  std::string Reduce(Link* from = 0) {
    if (!from) {
      from = GetRoot();
    }
    if (from->isLeaf()) {
      std::cerr << from->GetNode()->GetName() << '\n';
      return from->GetNode()->GetName();
    } else {
      std::string name = "None";
      for (Link* link = from->Next(); link != from; link = link->Next()) {
        std::string tmp = Reduce(link->Out());
        if (tmp == "diff") {
          name = "diff";
        } else if (name == "None") {
          name = tmp;
        } else if (name != tmp) {
          name = "diff";
        }
      }
      std::cerr << '\t' << name << '\n';
      from->GetNode()->SetName(name);
      return name;
    }
    return "";
  }

  void Print(std::ostream& os, const Link* from = 0) {
    if (!from) {
      from = GetRoot();
    }
    if (from->isLeaf()) {
      os << from->GetNode()->GetName() << ":";
    } else {
      os << '(';
      for (const Link* link = from->Next(); link != from; link = link->Next()) {
        Print(os, link->Out());
        if (link->Next() != from) {
          os << ',';
        }
      }
      if (!from->isRoot()) {
        os << "):";
      } else {
        os << ")";
      }
    }
    if (from->isRoot()) {
      os << ";\n";
    } else {
      os << from->GetBranch()->GetName();
    }
  }

  void PrintReduced(std::ostream& os, const Link* from = 0) {
    if (!from) {
      from = GetRoot();
    }
    if (from->GetNode()->GetName() != "diff") {
      os << from->GetNode()->GetName();
    } else {
      os << '(';
      for (const Link* link = from->Next(); link != from; link = link->Next()) {
        PrintReduced(os, link->Out());
        if (link->Next() != from) {
          os << ',';
        }
      }
      os << ')';
    }
    if (from->isRoot()) {
      os << ";\n";
    }
  }

  void SetIndices() {
    Nlink = 0;
    Nnode = GetSize();
    Nbranch = 1;
    linkmap.clear();
    nodemap.clear();
    branchmap.clear();
    SetIndices(GetRoot(), Nlink, Nnode, Nbranch);
  }

  int GetNlink() { return Nlink; }

  int GetNbranch() { return Nbranch; }

  int GetNnode() { return Nnode; }

  // maybe

  const Node* GetNode(int index) { return nodemap[index]; }
  const Branch* GetBranch(int index) { return branchmap[index]; }

  Link* GetLink(int index) { return linkmap[index]; }

  double GetTotalLength() { return RecursiveTotalLength(GetRoot()); }

 protected:
  std::map<int, const Node*> nodemap;
  std::map<int, const Branch*> branchmap;
  std::map<int, Link*> linkmap;

  void CheckIndices(Link* from) {
    if (!from->isRoot()) {
      if (from->GetBranch() != branchmap[from->GetBranch()->GetIndex()]) {
        std::cerr << "branch index : " << from->GetBranch()->GetIndex() << '\n';
        exit(1);
      }
    } else {
      if (branchmap[0] != 0) {
        std::cerr << "root branch index\n";
        exit(1);
      }
    }

    if (from->GetNode() != nodemap[from->GetNode()->GetIndex()]) {
      std::cerr << "node index: " << from->GetNode()->GetIndex() << '\n';
      exit(1);
    }

    if (!from->isRoot()) {
      if (from->Out() != linkmap[from->Out()->GetIndex()]) {
        std::cerr << "link index : " << from->Out()->GetIndex() << '\n';
      }
    }
    if (from != linkmap[from->GetIndex()]) {
      std::cerr << "link index : " << from->GetIndex() << '\n';
    }

    for (const Link* link = from->Next(); link != from; link = link->Next()) {
      CheckIndices(link->Out());
    }
  }

  void SetIndices(Link* from, int& linkindex, int& nodeindex,
                  int& branchindex) {
    if (!from->isRoot()) {
      from->GetBranch()->SetIndex(branchindex);
      branchmap[branchindex] = from->GetBranch();
      branchindex++;
    }

    if (!from->isLeaf()) {
      from->GetNode()->SetIndex(nodeindex);
      nodemap[nodeindex] = from->GetNode();
      nodeindex++;
    } else {
      nodemap[from->GetNode()->GetIndex()] = from->GetNode();
    }

    if (!from->isRoot()) {
      from->Out()->SetIndex(linkindex);
      linkmap[linkindex] = from->Out();
      linkindex++;
    }
    from->SetIndex(linkindex);
    linkmap[linkindex] = from;
    linkindex++;

    for (const Link* link = from->Next(); link != from; link = link->Next()) {
      SetIndices(link->Out(), linkindex, nodeindex, branchindex);
    }
  }

  // returns 0 if not found
  // returns link if found (then found1 and found2 must
  Link* RecursiveGetLCA(Link* from, std::string tax1, std::string tax2,
                        bool& found1, bool& found2) {
    Link* ret = 0;
    if (from->isLeaf()) {
      found1 |= (from->GetNode()->GetName() == tax1);
      found2 |= (from->GetNode()->GetName() == tax2);
      if (!ret) {
        if (found1 && found2) {
          ret = from;
        }
      }
    } else {
      for (Link* link = from->Next(); link != from; link = link->Next()) {
        bool tmp1 = false;
        bool tmp2 = false;
        Link* ret2 = RecursiveGetLCA(link->Out(), tax1, tax2, tmp1, tmp2);
        found1 |= tmp1;
        found2 |= tmp2;
        if (ret2) {
          if (ret) {
            std::cerr << "error : found node twice\n";
            std::cerr << tax1 << '\t' << tax2 << '\n';
            ToStream(std::cerr, ret2->Out());
            std::cerr << '\n';
            ToStream(std::cerr, ret->Out());
            std::cerr << '\n';
            exit(1);
          }
          ret = ret2;
        }
      }
      if (!ret) {
        if (found1 && found2) {
          ret = from;
        }
      }
    }
    return ret;
  }

  Link* ParseGroup(std::string input, Link* from);
  // a group is an expression of one of the two following forms:
  //
  // 	(Body)Node_name
  // 	(Body)Node_name:Branch_name
  //
  // where Body is a list, and Node_name and Branch_name are 2 std::strings
  // Node_name may be an empty std::string
  //
  // Node_name cannot contain the ':' character, but Branch_name can
  // thus, if the group reads "(BODY)A:B:C"
  // then Node_name = "A" and Branch_name = "B:C"

  Link* ParseList(std::string input, Node* node);
  // a list is an expression of the form X1,X2,...Xn
  // where Xi is a group

  void RecursiveClone(const Link* from, Link* to);
  // Used by Tree(const Tree* from)
  // only clone the Links, and their mutual relations
  // does not copy the Node or Branch objects

  // deletes the link structure
  // does not delete the Node or Branch objects
  void RecursiveDelete(Link* from);

  void SetRoot(Link* link) { root = link; }


double GetLength(const Branch* branch){
  if (!branch) {
    return 0;
  }
  return atof(branch->GetName().c_str());
}

double RecursiveTotalLength(const Link* from)	{
	double total = 0;
	if (! from->isRoot())	{
		total += GetLength(from->GetBranch());
	}
	else	{
		if (GetLength(from->GetBranch()))	{
			std::cerr << "error: non null branch length for root\n";
			std::cerr << GetLength(from->GetBranch()) << '\n';
			exit(1);
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		total += RecursiveTotalLength(link->Out());
	}
	return total;
}

  // data fields
  // just 2 pointers, to the root and to a list of taxa
  Link* root;
  const TaxonSet* taxset;
  int Nlink;
  int Nnode;
  int Nbranch;
};



#endif  // SOURCES_TREE_H_
