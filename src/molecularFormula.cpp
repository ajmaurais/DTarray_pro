//
//  molecularFormula.cpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 1/5/18.
//  Copyright © 2018 Aaron Maurais. All rights reserved.
//

#include <molecularFormula.hpp>

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

//removes atoms from atomCountMap which have count of 0
void molFormula::Residue::removeZeros()
{
	AtomCountMapType::iterator it = atomCountMap.begin();
	while(it != atomCountMap.end())
	{
		if(it->second == 0){
			atomCountMap.erase(it ++);
		} else {
			++it;
		}
	}
}

void molFormula::Residue::initalize(molFormula::AtomMassMapType* _atomMassMap,
									const molFormula::HeaderType& _header,
									const std::vector<std::string>& _elems)
{
	//check that atom count std::vector is same len as header
	assert(_header.size() == _elems.size());
	
	atomMassMap = _atomMassMap;
	
	//initalize atom count map
	size_t len = _header.size();
	for(size_t i = 0; i < len; i++)
		atomCountMap[_header[i]] = utils::toInt(_elems[i]);
	
	calcMasses(); //calculate residue mass
	removeZeros(); //remove 0 atom counts
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
		
		if(elems[0] == "H" || elems[0] == "A")
		{
			if(elems.size() != 4)
				throw std::runtime_error("Bad atom mass ");
		
			if(elems[0] == "H"){
				continue;
			}
			else if(elems[0] == "A"){
				atomMassMap[elems[1]] = molFormula::Species(utils::toDouble(elems[2]),
															utils::toDouble(elems[3]));
			}
		}//end if
	}//end while
	return true;
}//end fxn

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

double molFormula::Residues::calcMass(std::string _seq, char avg_mono, bool _nterm, bool _cterm) const
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
			return -1;
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

std::string molFormula::Residues::calcFormula(std::string _seq, bool unicode,
											  bool _nterm, bool _cterm) const
{
	std::string formula = "";
	ResidueMapType::const_iterator it;
	
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
			return BAD_AMINO_ACID;
		resMapIt->second.combineAtomCountMap(atomCounts);
	}
	if(_cterm)
	{
		resMapIt = residueMap.find(C_TERM_STR);
		if(it == residueMap.end())
			throw std::runtime_error(C_TERM_STR + " not found in residueMap");
		resMapIt->second.combineAtomCountMap(atomCounts);
	}
	
	return getFormulaFromMap(atomCounts, unicode);
}

/*std::string molFormula::symbolToUnicode(std::string _symbol)
{
	if(!(utils::strContains(")", _symbol) && _symbol[0] == '('))
		return _symbol;
	std::string ret;
	
	size_t endParen = _symbol.find(")");
	
	std::string temp = _symbol.substr(1, endParen - 1);
	
	ret = utils::toSuperscript(utils::toInt(temp));
	ret += _symbol.substr(endParen + 1);
	
	return ret;
}*/

std::string molFormula::getFormulaFromMap(const molFormula::AtomCountMapType& atomCountMap, bool unicode)
{
	std::string formula;
	
	//make atom count map which can keep track of already printed atoms
	typedef std::pair<int, bool> PairType;
	typedef std::map<std::string, PairType> AtomCountGraphType;
	AtomCountGraphType atomCountGraph;
	for(AtomCountMapType::const_iterator it = atomCountMap.begin(); it != atomCountMap.end(); ++it)
		atomCountGraph[it->first] =  PairType(it->second, false);
	
	//first print atoms in FORMULA_RESIDUE_ORDER
	for(size_t i = 0; i < FORMULA_RESIDUE_ORDER_LEN; i++)
	{
		if(atomCountGraph[FORMULA_RESIDUE_ORDER[i]].first == 0)
		{
			atomCountGraph[FORMULA_RESIDUE_ORDER[i]].second = true;
			continue;
		}
		else if(atomCountGraph[FORMULA_RESIDUE_ORDER[i]].first == 1) {
			formula += FORMULA_RESIDUE_ORDER[i];
		}
		else {
			if(unicode){
				//formula += molFormula::symbolToUnicode(FORMULA_RESIDUE_ORDER[i]);
				formula += FORMULA_RESIDUE_ORDER[i];
				formula += utils::toSubscript(atomCountGraph[FORMULA_RESIDUE_ORDER[i]].first);
			}
			else {
				formula += FORMULA_RESIDUE_ORDER[i];
				formula += utils::toString(atomCountGraph[FORMULA_RESIDUE_ORDER[i]].first);
			}
		}
		atomCountGraph[FORMULA_RESIDUE_ORDER[i]].second = true;
	}
	
	//next itterate through atomCountGraph and print atoms not added to formula
	for(AtomCountGraphType::iterator it = atomCountGraph.begin();
		it != atomCountGraph.end(); ++it)
	{
		if(it->second.second) //check if atom has already been added to formula
			continue;
		
		formula += it->first;
		if(it->second.first > 1) {
			if(unicode)
				formula += utils::toSubscript(it->second.first);
			else formula += utils::toString(it->second.first);
		}
		it->second.second = true;
	}
	
	return formula;
}








