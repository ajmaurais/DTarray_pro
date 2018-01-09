//
//  molecularFormula.cpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 1/5/18.
//  Copyright © 2018 Aaron Maurais. All rights reserved.
//

#include "molecularFormula.hpp"

void molFormula::Residue::calcMasses()
{
	masses = molFormula::Species(0, 0);
	for(AtomCountMapType::const_iterator it = atomCountMap.begin(); it != atomCountMap.end(); ++it)
	{
		//check that query exists in atomMassMap
		if(atomMassMap->find(it->first) == atomMassMap->end())
			throw std::runtime_error(it->first + " not found in mass map!");
		
		//add mass * count to *this masses
		masses += ((*atomMassMap)[it->first] * it->second);
	}
}

void molFormula::Residue::initalize(molFormula::AtomMassMapType* _atomMassMap,
									const molFormula::HeaderType& _header,
									const std::vector<std::string>& _elems)
{
	//check that atom count vector is same len as header
	assert(_header.size() == _elems.size());
	
	atomMassMap = _atomMassMap;
	
	//initalize atom count map
	size_t len = _header.size();
	for(size_t i = 0; i < len; i++)
		atomCountMap[_header[i]] = utils::toInt(_elems[i]);
	
	//calculate residue mass
	calcMasses();
}

double molFormula::Residue::getMass(char avg_mono) const{
	switch(avg_mono){
		case 'a' :
			return getAvg();
			break;
		case 'm' :
			return getMono();
			break;
		default :
			throw std::runtime_error("only a or m are valid arguments!");
	}
}

bool molFormula::Residues::readAtomCountTable(std::string _atomCountTableLoc)
{
	atomCountTableLoc = _atomCountTableLoc;
	return readAtomCountTable();
}

bool molFormula::Residues::readAtomCountTable()
{
	if(atomCountTableLoc.empty())
		throw std::runtime_error("atomCountTableLoc must be specified");
	
	utils::File file;
	if(!file.read(atomCountTableLoc))
		return false;
	
	std::string line;
	std::vector<std::string> elems;
	while(!file.end())
	{
		line = file.getLine_skip_trim();
		line = line.substr(0, line.find(";"));
		utils::split(line, IN_DELIM, elems);
		utils::trimAll(elems);
		
		if(elems[0] == "H" || elems[0] == "R"){
			if(elems[0] == "H"){
				//remove first and second elements
				elems.erase(elems.begin(), elems.begin() + 2);
				atomCountHeader = elems;
			}
			else if(elems[0] == "R")
			{
				//get residue symbol
				std::string residue = elems[1];
				
				//remove first and second elements
				elems.erase(elems.begin(), elems.begin() + 2);
				
				if(atomCountHeader.size() != elems.size())
					throw std::runtime_error("bad atom count table");
				
				residueMap[residue] = molFormula::Residue(&atomMassMap, atomCountHeader, elems);
			}
		}
	}
	return true;
}

void molFormula::Residue::combineAtomCountMap(AtomCountMapType& _add) const
{
	for(AtomCountMapType::const_iterator it = atomCountMap.begin(); it != atomCountMap.end(); ++it)
	{
		if(_add.find(it->first) == _add.end())
			_add[it->first] = it->second;
		else _add[it->first] += it->second;
	}
}

bool molFormula::Residues::readAtomMassTable(std::string _massTableLoc)
{
	massTableLoc = _massTableLoc;
	return readAtomCountTable();
}

bool molFormula::Residues::readAtomMassTable()
{
	if(massTableLoc.empty())
		throw std::runtime_error("atomCountTableLoc must be specified");
	
	utils::File file;
	if(!file.read(massTableLoc))
		return false;
	
	std::string line;
	std::vector<std::string> elems;
	while(!file.end())
	{
		line = file.getLine_skip_trim();
		utils::split(line, IN_DELIM, elems);
		utils::trimAll(elems);
		
		if(elems[0] == "H" || elems[0] == "A"){
			if(elems.size() != 4)
				throw std::runtime_error("Bad atom mass ");
		
			if(elems[0] == "H"){
				continue;
			}
			else if(elems[0] == "A"){
				atomMassMap[elems[1]] = molFormula::Species(utils::toDouble(elems[2]),
															utils::toDouble(elems[3]));
			}
		}
	}
	return true;
}

bool molFormula::Residues::initalize(std::string _atomCountTableLoc, std::string _massTableLoc)
{
	atomCountTableLoc = _atomCountTableLoc;
	massTableLoc = _massTableLoc;
	
	return initalize();
}

bool molFormula::Residues::initalize()
{
	bool gootAtomMassTable = readAtomMassTable();
	bool goodAtomCountTable = readAtomCountTable();
	
	return gootAtomMassTable && goodAtomCountTable;
}

double molFormula::Residues::getMass(std::string _seq, char avg_mono, bool _nterm, bool _cterm) const
{
	double mass = 0;
	size_t len = _seq.length();
	ResidueMapType::const_iterator it;
	
	if(_nterm)
	{
		it = residueMap.find(N_TERM_STR);
		if(it == residueMap.end())
			throw std::runtime_error(N_TERM_STR + " not found in residueMap");
		mass += it->second.getMass(avg_mono);
	}
	
	for(size_t i = 0; i < len; i++)
	{
		std::string aaTemp = std::string(1, _seq[i]);
		it = residueMap.find(aaTemp);
		if(it == residueMap.end())
			throw std::runtime_error(aaTemp + " is not a valid amino acid!");
		mass += it->second.getMass(avg_mono);
	}
	
	if(_cterm)
	{
		it = residueMap.find(C_TERM_STR);
		if(it == residueMap.end())
			throw std::runtime_error(C_TERM_STR + " not found in residueMap");
		mass += it->second.getMass(avg_mono);
	}
	
	return mass;
}

std::string molFormula::Residues::getFormula(std::string _seq, bool _nterm, bool _cterm) const
{
	std::string formula = "";
	ResidueMapType::const_iterator it;
	
	HeaderType atoms (FORMULA_RESIDUE_ORDER, FORMULA_RESIDUE_ORDER + FORMULA_RESIDUE_ORDER_LEN);
	for(std::string::iterator it = _seq.begin(); it != _seq.end(); ++it)
	{
		std::string aaTemp = std::string(1, *it);
		bool found = false;
		for(size_t i = 0; i < FORMULA_RESIDUE_ORDER_LEN; i++)
		{
			if(FORMULA_RESIDUE_ORDER[i] == aaTemp)
			{
				found = true;
				break;
			}
		}
		if(!found)
			atoms.push_back(aaTemp);
	}
	
	AtomCountMapType atomCounts;
	ResidueMapType::const_iterator resMapIt;
	if(_nterm)
	{
		resMapIt = residueMap.find(N_TERM_STR);
		if(it == residueMap.end())
			throw std::runtime_error(N_TERM_STR + " not found in residueMap");
		resMapIt->second.combineAtomCountMap(atomCounts);
	}
	for(std::string::iterator it = _seq.begin(); it != _seq.end(); ++it)
	{
		std::string aaTemp = std::string(1, *it);
		resMapIt = residueMap.find(aaTemp);
		if(resMapIt == residueMap.end())
			throw std::runtime_error(aaTemp + " not foud in residue map!");
		resMapIt->second.combineAtomCountMap(atomCounts);
	}
	if(_cterm)
	{
		resMapIt = residueMap.find(C_TERM_STR);
		if(it == residueMap.end())
			throw std::runtime_error(C_TERM_STR + " not found in residueMap");
		resMapIt->second.combineAtomCountMap(atomCounts);
	}
	
	for(HeaderType::iterator it = atoms.begin(); it != atoms.end(); ++it)
	{
		if(atomCounts[*it] == 0)
			continue;
		else if(atomCounts[*it] == 1)
			formula += *it;
		else formula += *it + utils::toString(atomCounts[*it]);
	}
	
	return formula;
}








