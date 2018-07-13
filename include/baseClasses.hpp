//
//  parentClasses.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 11/4/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#pragma once

#include <params.hpp>
#include <map>
//#include <hashTable.hpp>

size_t const DATA_SIZE = 500;
std::string const BLANK_STR = "null";
std::string const BLANK_VAL = "-1";

class ProteinTemplate;
template <class _Tp> class ProteinDataTemplate;
template <class _Tp> class DBTemplate;
class SampleData;
class SampleData_peptide;
class SampleData_protein;

class ProteinTemplate{
protected:
	std::string ID, protein, description;
	
public:
	ProteinTemplate() {}
	~ProteinTemplate() {}
	
	inline bool operator == (const ProteinTemplate& comp) const{
		return comp.ID == ID;
	}
	inline bool operator == (std::string comp) const{
		return comp == ID;
	}
	
	std::string getID(){
		return ID;
	}
};

template <class _Tp>
class ProteinDataTemplate{
public:
	ProteinDataTemplate(params::Params* const _par) {
		par = _par;
		supDataAdded = false;
	}
	void initialize(const std::vector<_Tp>&, size_t, size_t*);
	ProteinDataTemplate(){}
	~ProteinDataTemplate(){}
	
protected:
	static size_t colSize;
	static params::Params* par;
	
	std::vector<_Tp> col;
	double avgMass, monoMass;
	std::string formula;
	std::string length;
	std::string sequence;
	std::string matchDirrection;
	static size_t* colIndex;
	bool supDataAdded;
};

template <class _Tp> size_t* ProteinDataTemplate<_Tp>::colIndex = nullptr;
template <class _Tp> size_t ProteinDataTemplate<_Tp>::colSize = 0;
template <class _Tp> params::Params* ProteinDataTemplate<_Tp>::par = nullptr;

template<class _Tp>
class DBTemplate{
protected:
	static size_t colIndex;
	//hashTable::HashTable <_Tp>* data;
	typedef std::map<std::string, _Tp> DataType;
	DataType data;
	
public:
	std::vector<std::string> colNames;
	
	//constructor
	DBTemplate(){
		//data = new hashTable::HashTable <_Tp>(DATA_SIZE);
	}
	DBTemplate(const params::Params& par, size_t dataSize){
		//data = new hashTable::HashTable <_Tp>(dataSize);
		
		for(int i = 0; i < par.getNumFiles(); i++)
			colNames.push_back(par.getFileColname(i));
	}
	~DBTemplate(){
		//delete data;
	}
};

template <class _Tp> size_t DBTemplate<_Tp>::colIndex = 0;

//stores the data pertaining to a specific filter file (or MS run) for each protein
class SampleData{
public:
	std::string colname;
	int count;
	int modPeptidesSC;
	
	//constructor
	SampleData (std::string _colName)
	{
		colname = _colName;
		count = 0;
		modPeptidesSC = 0;
	}
	SampleData(){
		colname = BLANK_STR;
		count = 0;
		modPeptidesSC = 0;
	}
	~SampleData() {}
	
	inline bool isNull() const{
		return count == 0;
	}
};

class SampleData_peptide : public SampleData {
public:
	std::string parentFile, scan, obsMH, xCorr;
	
	SampleData_peptide(std::string colName) : SampleData(colName)
	{
		parentFile = BLANK_STR;
		scan = BLANK_VAL;
		obsMH = BLANK_VAL;
		xCorr = BLANK_VAL;
	}
	SampleData_peptide() : SampleData() {
		parentFile = BLANK_STR;
		scan = BLANK_VAL;
		obsMH = BLANK_VAL;
		xCorr = BLANK_VAL;
	}
	~SampleData_peptide() {}
};

class SampleData_protein : public SampleData {
public:
	std::string coverage, sequenceCount;
	unsigned int modPeptides, uniquePeptides;
	
	SampleData_protein(std::string colName) : SampleData(colName){
		
		coverage = "0";
		sequenceCount = "0";
		modPeptides = 0;
		uniquePeptides = 0;
		
	}
	SampleData_protein() : SampleData(){
		coverage = "0";
		sequenceCount = "0";
		modPeptides = 0;
		uniquePeptides = 0;
	}
	~SampleData_protein() {}
};

/* parentClasses_hpp */
