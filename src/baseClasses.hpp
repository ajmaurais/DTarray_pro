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
#include "hashTable.hpp"

size_t const MAX_NUM_FILES = 50;
size_t const DATA_SIZE = 500;

using namespace std;

class ProteinTemplate;
class ProteinDataTemplate;
template <class T> class DBTemplate;

class ProteinTemplate{
protected:
	string ID, protein, description;
	
public:
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

class ProteinDataTemplate{
public:
	ProteinDataTemplate(FilterFileParams* const _par) {
		par = _par;
		supDataAdded = false;
	}
	ProteinDataTemplate(){}
	~ProteinDataTemplate(){
		//free(
	}
	
protected:
	static size_t colSize;
	static FilterFileParams* par;
	
	double avgMass, monoMass;
	string sequence;
	static size_t* colIndex;
	bool supDataAdded;
};

size_t* ProteinDataTemplate::colIndex = nullptr;
size_t ProteinDataTemplate::colSize = 0;
FilterFileParams* ProteinDataTemplate::par = nullptr;

template<class T>
class DBTemplate{
protected:
	size_t colIndex;
	hashTable::HashTable <T>* data;
	
public:
	string colNames[MAX_NUM_FILES];
	
	//constructor
	DBTemplate(){
		colIndex = 0;
		data = new hashTable::HashTable <T>(DATA_SIZE);
	}
	DBTemplate(const FilterFileParams& par){
		colIndex = 0;
		data = new hashTable::HashTable <T>(DATA_SIZE);
		
		for (int i = 0; i < par.numFiles; i++)
			colNames[i] = par.getFileColname(i);
	}
	~DBTemplate(){
		delete data;
	}
};

#endif /* parentClasses_hpp */
