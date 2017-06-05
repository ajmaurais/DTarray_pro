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
		line = utils::trim(line);
		
		vector<string> elems;
		utils::split(line, '\t', elems);
		
		symbol = elems[0];
		avgMass = utils::toDouble(elems[1]);
		monoMass = utils::toDouble(elems[2]);
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
	
	void AminoAcid::operator += (const AminoAcid& mod)
	{
		avgMass += mod.avgMass;
		monoMass += mod.avgMass;
	}
	
	/* #################### MWDB #################### */
	
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
		AminoAcid* aaTemp = aminoAcidsDB->getItem(temp.symbol);
		
		if(aaTemp == nullptr)
			return -1;
		
		if(avgMono == 0)
			return aaTemp->avgMass;
		else if(avgMono == 1)
			return aaTemp->monoMass;
		else return -1;
	}
	
	double MWDB::getMW(char a, int avgMono) const
	{
		return getMW(string(1, a), avgMono);
	}
	
	bool MWDB::readInAADB(string aaDB)
	{
		utils::File file;
		if(!file.read(aaDB))
			return false;
		string line;
		
		do{
			line = file.getLine_skip_trim();
			AminoAcid newAA(line);
			aminoAcidsDB->insert(newAA, newAA.symbol);
		} while(!file.end());
		
		return true;
	}
	
	bool MWDB::readIn(string staticMods, string aaDB)
	{
		utils::File file;
		
		if(!file.read(staticMods) || !readInAADB(aaDB))
			return false;
		
		string line;
		int i = 0;
		
		do{
			line = file.getLine_skip_trim();
			if(line == "<staticModifications>")
			{
				do{
					line = file.getLine_skip_trim();
					if (line != "</staticModifications>")
					{
						AminoAcid temp (line);
						if(!addStaticMod(temp))
							cout << endl << "Could not find " << temp.symbol << "in aaDB!";
					}
					i++;
					if(i > MAX_PARAM_ITTERATIONS)
						return false;
				} while(line != "</staticModifications>");
			}
		} while(!file.end() && line != "</staticModifications>");
		
		return true;
	}
	
	bool MWDB::addStaticMod(const AminoAcid& mod)
	{
		AminoAcid* aaTemp = aminoAcidsDB->getItem(mod.symbol);
		
		if(aaTemp == nullptr)
			return false;
		
		*aaTemp += mod;
		
		return true;
	}
	
	string SeqDB::getSequence(string id) const
	{
		Peptide* const temp = seqLibrary->getItem(id);
		if(temp == nullptr)
			return SEQ_NOT_FOUND;
		else return temp->getSequence();
	}
	
	bool SeqDB::readIn(string fname)
	{
		utils::File data(fname);
		
		if(!data.read(fname))
			return false;
		
		string line;
		Peptide temp;
		
		do{
			line = data.getLine_skip_trim();
			if(line[0] == '>' && !utils::strContains("Reverse", line))
			{
				temp.ID = getID(line);
				line = data.getLine_skip_trim();
				if(line[0] == '>')
					return false;
				temp.sequence = line;
				seqLibrary->insert(temp, temp.getID());
			}
		} while(!data.end());
		
		return true;
	}
	
	bool MWDB_Protein::readIn(string wd, const params::Params& params)
	{
		bool val1 = seqDB->readIn(params.mwDBFname);
		bool val2 = MWDB::readIn(wd + params.staticModsFname, params.aaDBfanme);
		
		return val1 && val2;
	}
}