//
//  parentClasses.hpp
//  DTarray_pro
// -----------------------------------------------------------------------------
// Copyright 2018 Aaron maurais
// -----------------------------------------------------------------------------
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
// -----------------------------------------------------------------------------
//

#pragma once

#include <params.hpp>
#include <map>

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
	std::string _ID, _protein, _description;
	
public:
	ProteinTemplate() {}
	~ProteinTemplate() {}
	
	inline bool operator == (const ProteinTemplate& comp) const{
		return comp._ID == _ID;
	}
	inline bool operator == (std::string comp) const{
		return comp == _ID;
	}
	
	std::string getID(){
		return _ID;
	}
};

template <class _Tp>
class ProteinDataTemplate{
public:
	ProteinDataTemplate(params::Params* const par) {
		_par = par;
		_supDataAdded = false;
	}
	void initialize(const std::vector<_Tp>&, size_t, size_t*);
	ProteinDataTemplate(){}
	~ProteinDataTemplate(){}
	
protected:
	static size_t _colSize;
	static params::Params* _par;
	
	std::vector<_Tp> _col;
	double _avgMass, _monoMass;
	std::string _formula;
	std::string _length;
	std::string _sequence;
	std::string _matchDirrection;
	static size_t* _colIndex;
	bool _supDataAdded;
};

template <class _Tp> size_t* ProteinDataTemplate<_Tp>::_colIndex = nullptr;
template <class _Tp> size_t ProteinDataTemplate<_Tp>::_colSize = 0;
template <class _Tp> params::Params* ProteinDataTemplate<_Tp>::_par = nullptr;

template<class _Tp>
class DBTemplate{
protected:
	static size_t _colIndex;
	typedef std::map<std::string, _Tp> DataType;
	DataType data;
	
public:
	std::vector<std::string> _colNames;
	
	//constructor
	DBTemplate(){}
	DBTemplate(const params::Params& par)
	{
		for(int i = 0; i < par.getNumFiles(); i++)
			_colNames.push_back(par.getFileColname(i));
	}
	~DBTemplate(){}
};

template <class _Tp> size_t DBTemplate<_Tp>::_colIndex = 0;

//stores the data pertaining to a specific filter file (or MS run) for each protein
class SampleData{
public:
	std::string _colname;
	int _count;
	int _modPeptidesSC;
	
	//constructor
	SampleData (std::string _colName)
	{
		_colname = _colName;
		_count = 0;
		_modPeptidesSC = 0;
	}
	SampleData(){
		_colname = BLANK_STR;
		_count = 0;
		_modPeptidesSC = 0;
	}
	~SampleData() {}
	
	inline bool isNull() const{
		return _count == 0;
	}
};

class SampleData_peptide : public SampleData {
public:
	std::string _parentFile, _scan, _obsMH, _xCorr;
	
	SampleData_peptide(std::string colName) : SampleData(colName)
	{
		_parentFile = BLANK_STR;
		_scan = BLANK_VAL;
		_obsMH = BLANK_VAL;
		_xCorr = BLANK_VAL;
	}
	SampleData_peptide() : SampleData() {
		_parentFile = BLANK_STR;
		_scan = BLANK_VAL;
		_obsMH = BLANK_VAL;
		_xCorr = BLANK_VAL;
	}
	~SampleData_peptide() {}
};

class SampleData_protein : public SampleData {
public:
	std::string _coverage, sequenceCount;
	int _modPeptides, _uniquePeptides;
	
	SampleData_protein(std::string colName) : SampleData(colName){
		
		_coverage = "0";
		sequenceCount = "0";
		_modPeptides = 0;
		_uniquePeptides = 0;
		
	}
	SampleData_protein() : SampleData(){
		_coverage = "0";
		sequenceCount = "0";
		_modPeptides = 0;
		_uniquePeptides = 0;
	}
	~SampleData_protein() {}
};

/* parentClasses_hpp */
