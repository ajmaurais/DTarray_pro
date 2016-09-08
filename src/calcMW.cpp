//
//  calcMW.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/5/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

/* #################### Peptide #################### */

Peptide::Peptide() {}

bool Peptide::operator == (string comp) const
{
	return strComp(comp, ID) == 0;
}

void Peptide::operator = (const Peptide& p)
{
	ID = p.ID;
	sequence = p.sequence;
}

string Peptide::getID() const
{
	return ID;
}

/* #################### AminoAcid #################### */

AminoAcid::AminoAcid() {}

AminoAcid::AminoAcid(string line)
{
	size_t endOfLine = line.find(";");
	line = line.substr(0, endOfLine);
	line = trim(line);
	
	vector<string> elems;
	split(line, '\t', elems);
	
	symbol = elems[0];
	avgMass = stod(elems[1]);
	monoMass = stod(elems[2]);
}

AminoAcid::AminoAcid(string str, double num1, double num2)
{
	symbol = str;
	avgMass = num1;
	monoMass = num2;
}

AminoAcid::AminoAcid(string str, double num1)
{
	symbol = str;
	avgMass = num1;
	monoMass = num1;
}

bool AminoAcid::operator < (const AminoAcid& compp) const
{
	return strComp(compp.symbol, symbol) < 0;
}

bool AminoAcid::operator > (const AminoAcid& compp) const
{
	return strComp(compp.symbol, symbol) > 0;
}

bool AminoAcid::operator == (const AminoAcid& compp) const
{
	return strComp(compp.symbol, symbol) == 0;
}

void AminoAcid::operator += (const AminoAcid& mod)
{
	avgMass += mod.avgMass;
	monoMass += mod.avgMass;
}

/* #################### mwNode #################### */

mwNode::mwNode()
{
	new struct AminoAcid;
	left = nullptr;
	right = nullptr;
}

/* #################### AATree #################### */

AATree::AATree()
{
	root = nullptr;
}

AATree::~AATree()
{
	destroyTree();
}

void AATree::destroyTree(mwNode *leaf)
{
	if(leaf != nullptr)
	{
		destroyTree(leaf->left);
		destroyTree(leaf->right);
		delete leaf;
	}
}

void AATree::destroyTree()
{
	destroyTree(root);
}

void AATree::insert(const AminoAcid& a, mwNode *leaf)
{
	if(a < leaf->aa)
	{
		if(leaf->left != nullptr)
			insert(a, leaf->left);
		else
		{
			leaf->left = new mwNode;
			leaf->left->aa = a;
			leaf->left->left=nullptr;
			leaf->left->right=nullptr;
		}
	}
	else if(a > leaf->aa || a == leaf->aa)
	{
		if(leaf->right != nullptr)
			insert(a, leaf->right);
		else
		{
			leaf->right = new mwNode;
			leaf->right->aa = a;
			leaf->right->left=nullptr;
			leaf->right->right=nullptr;
		}
	}
}

void AATree::insert(const AminoAcid& a)
{
	if(root!=nullptr)
		insert(a, root);
	else
	{
		root = new mwNode;
		root->aa = a;
		root->left = nullptr;
		root->right = nullptr;
	}
}

mwNode* AATree::search(const AminoAcid& a, mwNode *leaf) const
{
	if(leaf!=nullptr)
	{
		if(a == leaf->aa)
			return leaf;
		if(a < leaf->aa)
			return search(a, leaf->left);
		else
			return search(a, leaf->right);
	}
	else return nullptr;
}

mwNode* AATree::search(const AminoAcid& a) const
{
	return search(a, root);
}

double AATree::getMW(string a, int avgMono) const
{
	AminoAcid temp (a, 0, 0);
	mwNode* nTemp = search(temp);
	
	if(nTemp == nullptr)
		return -1;

	if(avgMono == 0)
		return nTemp->aa.avgMass;
	else if(avgMono == 1)
		return nTemp->aa.monoMass;
	else return -1;
}

double AATree::getMW(char a, int avgMono) const
{
	return getMW(string(1, a), avgMono);
}

bool AATree::readInAADB(string aaDB)
{
	ifstream inF (aaDB.c_str());
	
	if(!inF)
		return false;
	
	string line;
	
	do{
		getLineTrim(inF, line);
		if(isCommentLine(line) || line.empty())
			continue;
		else insert(AminoAcid(line));
	}while(!inF.eof());
	
	return true;
}

bool AATree::readInAAs(string staticMods, string aaDB)
{
	ifstream inF (staticMods.c_str());
	
	if(!inF || !readInAADB(aaDB))
		return false;
	
	string line;
	int i = 0;
	
	do{
		getLineTrim(inF, line);
		if(isCommentLine(line) || line.empty())
			continue;
		
		if(line == "<staticModifications>")
		{
			do{
				getLineTrim(inF, line);
				if(isCommentLine(line) || line.empty())
					continue;
				if (line != "</staticModifications>")
				{
					AminoAcid temp (line);
					if(!addStaticMod(temp))
						cout << "Could not find " << temp.symbol << "in aaDB! Using default mass." << endl;
				}
				i++;
				if(i > MAX_PARAM_ITTERATIONS)
					return false;
			} while(line != "</staticModifications>");
		}
	} while(!inF.eof() && line != "</staticModifications>");
	
	return true;
}

bool AATree::addStaticMod(const AminoAcid& mod)
{
	mwNode* nTemp = search(mod);
	
	if(nTemp == nullptr)
		return false;
	
	nTemp->aa += mod;
	
	return true;
}

/* #################### MWDB #################### */

MWDB::MWDB()
{
	new HashTable;
	new AATree;
}

MWDB::~MWDB()
{
	//peptideLibrary.destroyTable();
	//aminoAcidsDB.destroyTree();
}

double MWDB::calcMW(string sequence, int avgMono) const
{
	double mass = 0;
	int len = int(sequence.length());
	
	mass += aminoAcidsDB.getMW("N_term", avgMono);
	
	for(int j = 0; j < len; j++)
		mass += aminoAcidsDB.getMW(sequence[j], avgMono);
	
	mass += aminoAcidsDB.getMW("C_term", avgMono);
	
	return mass;
}

string MWDB::getSequence(string id) const
{
	Item* temp = peptideLibrary.getItem(id);
	if(temp == nullptr)
		return SEQ_NOT_FOUND;
	else return temp->val.sequence;
}

bool MWDB::readInPeptides(string fname)
{
	ifstream inF (fname.c_str());
	
	if(!inF)
		return false;
	
	string line;
	Peptide temp;
	
	do{
		getLineTrim(inF, line);
		if(isCommentLine(line) || line.empty()) //skip line if is comment line
			continue;
		else
		{
			if(line[0] == '>')
			{
				temp.ID = getID(line);
				getLineTrim(inF, line);
				temp.sequence = line;
				peptideLibrary.insert(temp);
			}
		}
	} while(!inF.eof());

	return true;
}

bool MWDB::readIn(string wd, const FilterFileParams& params)
{
	return readInPeptides(wd + params.peptideDBfname) && aminoAcidsDB.readInAAs(wd + params.staticModsFname, params.aaDBfanme);
}
