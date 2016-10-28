//
//  DTarray_AJM.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 3/25/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef DTarray_AJM_hpp
#define DTarray_AJM_hpp

//#define nullptr NULL //to compile on pleiades this line must be included

#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <string>
#include <algorithm>
#include "utils.hpp"
#include "FilterFile.hpp"
#include "BinTree.cpp"
#include "hashTable.cpp"
#include "calcMW.hpp"

using namespace std;

/******************************/
/* globally scoped constants */
/*****************************/

/* dtafilter.cpp */
bool const INCLUDE_FULL_DESCRIPTION = true;
string const DEFAULT_COL_NAMES [] = {"Full_description", "ID", "Protein", "Description", "Mass(Da)", "subcellular_location"};
int const DEFAULT_COL_NAMES_LENGTH = 5;
string const COLUMN_HEADER_LINE_ELEMENTS[] = {"Unique", "FileName", "XCorr", "DeltCN", "Conf%", "M+H+",
	"CalcM+H+", "TotalIntensity", "SpR", "ZScore", "IonProportion", "Redundancy", "Sequence"};
int const COLUMN_HEADER_LINE_ELEMENTS_LENGTH = 13;
string const PARAM_ERROR_MESSAGE = " is an invalid arguement for: ";
int const MAX_PARAM_ITTERATIONS = 100;
string const SEQ_NOT_FOUND = "SEQUENCE_NOT_FOUND_IN_DB";

//editable params for DB output format
string const DEFAULT_COL_NAMES_DB [] = {"Full_description", "ID", "Protein", "Description", "Mass(Da)", "Long_sample_name",
	"Spectral_counts", "Sample", "Replicate"};
int const DEFAULT_COL_NAMES_DB_LENGTH = 9;
string const DEFAULT_COL_NAMES_DB_LOC [] = {"Protein","ID", "Protein", "Description", "Mass(Da)", "subcellular_location", "Long_sample_name",
	"Spectral_counts", "Sample", "Replicate"};
int const DEFAULT_COL_NAMES_DB_LOC_LENGTH = 10;
string const UNIQUE_PEPTIDE_HEADERS[] = {"SC", "Unique_pep_SC"};
string const MWCALC_HEADERS [] = {"avg_mass", "monoisotopic_mass", "sequence"};
int const MWCALC_HEADERS_LENGTH = 2;

/* subcelluarLoc.cpp */
string const LOC_NOT_FOUND = "NOT_FOUND_IN_DB";


/**********************/
/* class definitions */
/*********************/

class Protein;
class Proteins;
class ProteinTemplate;
template <class T> class DBTemplate;


/* ############ parent classes ############*/
class ProteinTemplate{
protected:
	string ID, protein, description, loc;
	
public:
	inline bool operator == (const ProteinTemplate&) const;
	inline bool operator == (string) const;
	inline bool operator > (const ProteinTemplate&) const;
	inline bool operator < (const ProteinTemplate&) const;
	
	string get_ID(){
		return ID;
	}
};

class ProteinDataTemplate{
public:
	void consolidate(const ProteinDataTemplate&, int);
	void calcMW(const mwDB::MWDB&, string);
	
protected:
	vector <FilterFile> col;
	
	string MW;
	double avgMass, monoMass;
	string sequence;
};

template<class T>
class DBTemplate{
protected:
	vector <T>* data;
	int colIndex;
	
	//modifers
	long insert(const T&);
	
public:
	vector<string> colNames;
	
	//constructor
	DBTemplate(){
		colIndex = 0;
		data = new vector<T>;
	}
	DBTemplate(const FilterFileParams& par){
		colIndex = 0;
		data = new vector<T>;
		
		for (int i = 0; i < par.numFiles; i++)
			colNames.push_back(par.getFileColname(i));
	}
	~DBTemplate(){
		delete data;
	}
	
	//properities
	virtual bool writeOut(string, const FilterFileParams&) const =0;
	virtual bool writeOutDB(string, const FilterFileParams&) const =0;
	
	//modifiers
	void calcMW(const mwDB::MWDB&);
};

/* #################### subcelluarLoc.cpp #################### */
class DBProtein: public ProteinTemplate{
	friend class Protein;
	friend class Proteins;
	
	//modifers
	void clear();
	
public:
	//constructor
	DBProtein(string);
	DBProtein();
	
	//modifers
	void operator = (const DBProtein&);
};

//stores data for each protein found in filter file
class Protein : public ProteinTemplate , public ProteinDataTemplate {
friend class Proteins;
private:
	
	//modifier
	void initialize(const vector<string>&);
	bool getProteinData(string, int);
	inline void getProteinAndDescr(string);
	DBProtein toDBprotein() const;
	void addLoc(string);
	
	string fullDescription, matchDirrection;
};

//stores data for all proteins found in DTA filter files
class Proteins : public DBTemplate<Protein>{
	hashTable::HashTable<DBProtein>* locDB;
	
	//modifers
	bool readIn(string, const FilterFileParam&, bool);
	//long insert(const Protein&);
	
public:
	//constructor
	Proteins(const FilterFileParams& pars) : DBTemplate<Protein>(pars){
		locDB = new hashTable::HashTable<DBProtein>;
	}
	Proteins() : DBTemplate<Protein>(){
		locDB = new hashTable::HashTable<DBProtein>;
	}
	~Proteins(){
		delete locDB;
	}
	
	bool readInLocDB(string);
	string locSearch(const DBProtein&) const;
	
	//properities
	bool writeOut(string, const FilterFileParams&) const;
	bool writeOutDB(string, const FilterFileParams&) const;
	
	//modifiers
	bool readIn(string, const FilterFileParams&);
	void addSubcelluarLoc();
	void addSeq(const mwDB::SeqDB&);
};

/*************/
/* functions */
/*************/

/* dtafilter.cpp */
bool isColumnHeaderLine(const vector<string>&);
string parseSample(string, string, bool);
int parsePeptideSC(string);
//void parsePeptide(string);
string parseReplicate(string);
string getID(string);

#endif /* DTarray_AJM_hpp */

