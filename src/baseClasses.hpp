//
//  parentClasses.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 11/4/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef parentClasses_hpp
#define parentClasses_hpp

#include "FilterFile.hpp"
#include "../lib/hashTable.hpp"

size_t const DATA_SIZE = 500;

using namespace std;

class ProteinTemplate;
template <class _Tp> class ProteinDataTemplate;
template <class _Tp> class DBTemplate;

class ProteinTemplate{
protected:
	string ID, protein, description;
	
public:
	ProteinTemplate() {}
	~ProteinTemplate() {}
	
	inline bool operator == (const ProteinTemplate& comp) const{
		return comp.ID == ID;
	}
	inline bool operator == (string comp) const{
		return comp == ID;
	}
	inline bool operator > (const ProteinTemplate& comp) const{
		return comp.ID > ID;
	}
	inline bool operator < (const ProteinTemplate& comp) const{
		return comp.ID < ID;
	}
	
	string get_ID(){
		return ID;
	}
};

template <class _Tp>
class ProteinDataTemplate{
public:
	ProteinDataTemplate(FilterFileParams* const _par) {
		par = _par;
		supDataAdded = false;
	}
	void initialize(const vector<_Tp>&, size_t, size_t*);
	ProteinDataTemplate(){}
	~ProteinDataTemplate(){}
	
protected:
	static size_t colSize;
	static FilterFileParams* par;
	static size_t supInfoNum;
	
	vector<_Tp> col;
	double avgMass, monoMass;
	string sequence;
	static size_t* colIndex;
	bool supDataAdded;
};

template <class _Tp> size_t* ProteinDataTemplate<_Tp>::colIndex = nullptr;
template <class _Tp> size_t ProteinDataTemplate<_Tp>::colSize = 0;
template <class _Tp> FilterFileParams* ProteinDataTemplate<_Tp>::par = nullptr;
template <class _Tp> size_t ProteinDataTemplate<_Tp>::supInfoNum = 0;

template<class _Tp>
class DBTemplate{
protected:
	size_t colIndex;
	hashTable::HashTable <_Tp>* data;
	
public:
	vector<string> colNames;
	
	//constructor
	DBTemplate() : colIndex(0) {
		data = new hashTable::HashTable <_Tp>(DATA_SIZE);
	}
	DBTemplate(const FilterFileParams& par, size_t dataSize) : colIndex(0){
		data = new hashTable::HashTable <_Tp>(dataSize);
		
		for (int i = 0; i < par.numFiles; i++)
			colNames.push_back(par.getFileColname(i));
	}
	~DBTemplate(){
		delete data;
	}
};

#endif /* parentClasses_hpp */
