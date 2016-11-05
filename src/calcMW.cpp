//
//  calcMW.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/5/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "calcMW.hpp"

namespace mwDB{
	
	/* #################### Peptide #################### */
	
	Peptide::Peptide() {}
	
	bool Peptide::operator == (string comp) const
	{
		return comp == ID;
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
	
	string Peptide::getSequence() const
	{
		return sequence;
	}
	
	/* #################### AminoAcid #################### */
	
	AminoAcid::AminoAcid() {}
	
	AminoAcid::AminoAcid(string line)
	{
		size_t endOfLine = line.find(";");
		line = line.substr(0, endOfLine);
		line = util::trim(line);
		
		vector<string> elems;
		util::split(line, '\t', elems);
		
		symbol = elems[0];
		avgMass = util::toDouble(elems[1]);
		monoMass = util::toDouble(elems[2]);
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
		return symbol < compp.symbol;
	}
	
	bool AminoAcid::operator > (const AminoAcid& compp) const
	{
		return symbol > compp.symbol;
	}
	
	bool AminoAcid::operator == (const AminoAcid& compp) const
	{
		return symbol == compp.symbol;
	}
	
	void AminoAcid::operator += (const AminoAcid& mod)
	{
		avgMass += mod.avgMass;
		monoMass += mod.avgMass;
	}
	
	/* #################### MWDB #################### */
	
	MWDB::MWDB()
	{
		aminoAcidsDB = new binTree::BinTree <AminoAcid>;
	}
	
	MWDB::~MWDB()
	{
		delete aminoAcidsDB;
	}
	
	double MWDB::calcMW(string sequence, int avgMono) const
	{
		double mass = 0;
		int len = int(sequence.length());
		double temp;
		
		mass += getMW("N_term", avgMono);
		
		for(int j = 0; j < len; j++)
		{
			temp = getMW(sequence[j], avgMono);
			if(temp == -1)
				return -1;
			else mass += temp;
		}
		
		mass += getMW("C_term", avgMono);
		
		return mass;
	}
	
	double MWDB::getMW(string a, int avgMono) const
	{
		AminoAcid temp (a, 0, 0);
		binTree::Node<AminoAcid>* nTemp = aminoAcidsDB->search(temp);
		
		if(nTemp == nullptr)
			return -1;
		
		if(avgMono == 0)
			return nTemp->val.avgMass;
		else if(avgMono == 1)
			return nTemp->val.monoMass;
		else return -1;
	}
	
	double MWDB::getMW(char a, int avgMono) const
	{
		return getMW(string(1, a), avgMono);
	}
	
	bool MWDB::readInAADB(string aaDB)
	{
		ifstream inF (aaDB.c_str());
		
		if(!inF)
			return false;
		
		string line;
		
		do{
			util::getLineTrim(inF, line);
			if(util::isCommentLine(line) || line.empty())
				continue;
			else aminoAcidsDB->insert(AminoAcid(line));
		}while(!inF.eof());
		
		return true;
	}
	
	bool MWDB::readIn(string staticMods, string aaDB)
	{
		ifstream inF (staticMods.c_str());
		
		if(!inF || !readInAADB(aaDB))
			return false;
		
		string line;
		int i = 0;
		
		do{
			util::getLineTrim(inF, line);
			if(util::isCommentLine(line) || line.empty())
				continue;
			
			if(line == "<staticModifications>")
			{
				do{
					util::getLineTrim(inF, line);
					if(util::isCommentLine(line) || line.empty())
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
	
	bool MWDB::addStaticMod(const AminoAcid& mod)
	{
		binTree::Node<AminoAcid>* nTemp = aminoAcidsDB->search(mod);
		
		if(nTemp == nullptr)
			return false;
		
		nTemp->val += mod;
		
		return true;
	}
	
	SeqDB::SeqDB()
	{
		seqLibrary = new hashTable::HashTable<mwDB::Peptide>;
	}
	
	SeqDB::~SeqDB()
	{
		delete seqLibrary;
	}
	
	string SeqDB::getSequence(string id) const
	{
		hashTable::Node<Peptide>* temp = seqLibrary->getItem(id);
		if(temp == nullptr)
			return SEQ_NOT_FOUND;
		else return temp->val.getSequence();
	}
	
	bool SeqDB::readIn(string fname)
	{
		ifstream inF (fname.c_str());
		
		if(!inF)
			return false;
		
		string line;
		Peptide temp;
		
		do{
			util::getLineTrim(inF, line);
			if(util::isCommentLine(line) || line.empty()) //skip line if is comment line
				continue;
			else
			{
				if(line[0] == '>')
				{
					temp.ID = getID(line);
					util::getLineTrim(inF, line);
					temp.sequence = line;
					seqLibrary->insert(temp, temp.getID());
				}
			}
		} while(!inF.eof());
		
		return true;
	}
	
	bool MWDB_Protein::readIn(string wd, const FilterFileParams& params)
	{
		bool val1 = seqDB->readIn(wd + params.mwDBFname);
		bool val2 = MWDB::readIn(wd + params.staticModsFname, params.aaDBfanme);
		
		return val1 && val2;
	}
}
