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
	
	string getID(){
		return ID;
	}
};

template <class _Tp>
class ProteinDataTemplate{
public:
	ProteinDataTemplate(filterFile::FilterFileParams* const _par) {
		par = _par;
		supDataAdded = false;
	}
	void initialize(const vector<_Tp>&, size_t, size_t*);
	ProteinDataTemplate(){}
	~ProteinDataTemplate(){}
	
protected:
	static size_t colSize;
	static filterFile::FilterFileParams* par;
	
	vector<_Tp> col;
	double avgMass, monoMass;
	string length;
	string sequence;
	static size_t* colIndex;
	bool supDataAdded;
};

template <class _Tp> size_t* ProteinDataTemplate<_Tp>::colIndex = nullptr;
template <class _Tp> size_t ProteinDataTemplate<_Tp>::colSize = 0;
template <class _Tp> filterFile::FilterFileParams* ProteinDataTemplate<_Tp>::par = nullptr;

template<class _Tp>
class DBTemplate{
protected:
	static size_t colIndex;
	hashTable::HashTable <_Tp>* data;
	
public:
	vector<string> colNames;
	
	//constructor
	DBTemplate(){
		data = new hashTable::HashTable <_Tp>(DATA_SIZE);
	}
	DBTemplate(const filterFile::FilterFileParams& par, size_t dataSize){
		data = new hashTable::HashTable <_Tp>(dataSize);
		
		for(int i = 0; i < par.numFiles; i++)
			colNames.push_back(par.getFileColname(i));
	}
	~DBTemplate(){
		delete data;
	}
};

template <class _Tp> size_t DBTemplate<_Tp>::colIndex = 0;

#endif /* parentClasses_hpp */
