//
//  DTarray_AJM.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 3/25/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef DTarray_AJM_hpp
#define DTarray_AJM_hpp

//deal with older c++ compilers
#if (__cplusplus == 199711L || __cplusplus == 1)
	#define nullptr NULL
#endif

#define OUT_DELIM '\t'
#define IN_DELIM '\t'

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
size_t const DEFAULT_COL_NAMES_LENGTH = 5;
string const COLUMN_HEADER_LINE_ELEMENTS[] = {"Unique", "FileName", "XCorr", "DeltCN", "Conf%", "M+H+",
	"CalcM+H+", "TotalIntensity", "SpR", "ZScore", "IonProportion", "Redundancy", "Sequence"};
size_t const COLUMN_HEADER_LINE_ELEMENTS_LENGTH = 13;
string const PARAM_ERROR_MESSAGE = " is an invalid arguement for: ";
size_t const MAX_PARAM_ITTERATIONS = 100;
string const SEQ_NOT_FOUND = "SEQUENCE_NOT_FOUND_IN_DB";
size_t const MAX_NUM_FILES = 50;

//editable params for DB output format
string const DEFAULT_COL_NAMES_DB [] = {"Full_description", "ID", "Protein", "Description", "Mass(Da)", "Long_sample_name",
	"Spectral_counts", "Sample", "Replicate"};
size_t const DEFAULT_COL_NAMES_DB_LENGTH = 9;
string const DEFAULT_COL_NAMES_DB_LOC [] = {"Protein","ID", "Protein", "Description", "Mass(Da)", "subcellular_location", "Long_sample_name",
	"Spectral_counts", "Sample", "Replicate"};
size_t const DEFAULT_COL_NAMES_DB_LOC_LENGTH = 10;
string const SUP_INFO_HEADERS[] = {"SC", "Unique_pep_SC", "coverage"};
string const MWCALC_HEADERS [] = {"avg_mass", "monoisotopic_mass", "sequence"};
size_t const MWCALC_HEADERS_LENGTH = 2;
string const DEFALUT_PEPTIDE_COLNAMES [] = {"protein_ID", "parent_protein", "protein_description", "sequence", "charge", "unique", "calcMH"};
size_t const DEFALUT_PEPTIDE_COLNAMES_LEN = 7;
string const DEFALUT_PEPTIDE_DB_COLNAMES [] = {"protein_ID", "parent_protein", "protein_description", "sequence", "charge", "unique", "calcMH", "obsMH", "scan", "parent_file", "Long_sample_name", "Spectral_counts", "Sample", "Replicate"};
size_t const DEFALUT_PEPTIDE_DB_COLNAMES_LEN = 12;

/* subcelluarLoc.cpp */
string const LOC_NOT_FOUND = "NOT_FOUND_IN_DB";

/**********************/
/* class definitions */
/*********************/

class Protein;
class Proteins;
class ProteinTemplate;
class Peptide;
class ProteinDataTemplate;
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
	ProteinDataTemplate(FilterFileParams* const _par) {
		par = _par;
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
		data = new hashTable::HashTable <T>;
	}
	DBTemplate(const FilterFileParams& par){
		colIndex = 0;
		data = new hashTable::HashTable <T>;
		
		for (int i = 0; i < par.numFiles; i++)
			colNames[i] = par.getFileColname(i);
	}
	~DBTemplate(){
		delete data;
	}
	
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

class Peptide : public ProteinDataTemplate {
	friend class Proteins;
	friend class Peptides;
public:
	Peptide (FilterFileParams * par) : ProteinDataTemplate(par) {}
	Peptide () : ProteinDataTemplate() {}
	~Peptide() {
		//delete col;
	}
	
	//modifers
	inline void clear();
	void calcMW(mwDB::MWDB* const, string);
	
	//properities
	bool operator == (const Peptide& comp) const {
		return comp.key == key;
	}
	string makeKey() const {
		return proteinID + "_" + sequence + "_" + charge;
	}
	inline bool operator > (const Peptide& comp) const{
		return key > comp.key;
	}
	inline bool operator < (const Peptide& comp) const{
		return key < comp.key;
	}
	void consolidate(const Peptide&);
	void calcMW(mwDB::MWDB* const);
	void write(ofstream&) const;
	
private:
	string key;
	string proteinID, calcMH, scan, protein, description, charge;
	bool unique;
	FilterFileData_peptide* col;
	
	void initialize(FilterFileData_peptide* const, size_t, size_t*);
	void parsePeptide(const string&);
	void addSupData(mwDB::MWDB_Protein* const);
};

class Peptides : public DBTemplate<Peptide> {
	friend class Proteins;
private:
	mwDB::MWDB* mwdb;
public:
	Peptides(const FilterFileParams& pars) : DBTemplate<Peptide>(pars) {}
	Peptides() : DBTemplate<Peptide>(){}
	~Peptides(){}
	
	//properities
	bool writeOut(string, const FilterFileParams&) const;
	bool writeOutDB(string, const FilterFileParams&);
};

//stores data for each protein found in filter file
class Protein : public ProteinTemplate , public ProteinDataTemplate {
friend class Proteins;
private:
	string MW;
	string fullDescription, matchDirrection;
	int sequenceCount;
	FilterFileData_protein* col;
	
	//modifier
	bool getProteinData(string, size_t);
	inline void getProteinAndDescr(string);
	DBProtein toDBprotein() const;
	void initialize(FilterFileData_protein* const, size_t, size_t*);
	
	void addSupData(hashTable::HashTable<DBProtein>* const, mwDB::MWDB_Protein* const, mwDB::SeqDB* const);
	inline void clear();
	
public:
	Protein(FilterFileParams* const pars) : ProteinDataTemplate(pars){}
	Protein() : ProteinDataTemplate() {}
	~Protein(){
		//delete col;
	}
	
	void consolidate(const Protein&);
	void calcMW(mwDB::MWDB_Protein* const);
	void addSeq(mwDB::SeqDB* const);
	void addLoc(hashTable::HashTable<DBProtein>* const);
	void write(ofstream&) const;
};

//stores data for all proteins found in DTA filter files
class Proteins : public DBTemplate<Protein>{
	hashTable::HashTable<DBProtein>* locDB;
	mwDB::MWDB_Protein* mwdb;
	mwDB::SeqDB* seqDB;
	
	//modifers
	bool readIn(string, FilterFileParams* const, const FilterFileParam&, FilterFileData_protein* const, FilterFileData_peptide* const, Peptides* const);
public:
	//constructor
	Proteins(const FilterFileParams& pars) : DBTemplate<Protein>(pars){
		locDB = nullptr;
		seqDB = nullptr;
		mwdb = nullptr;
	}
	Proteins() : DBTemplate<Protein>(){
		locDB = nullptr;
		seqDB = nullptr;
		mwdb = nullptr;
	}
	~Proteins(){
		delete locDB;
	}
	
	bool readInLocDB(string);
	bool readInMWdb(string, const FilterFileParams&);
	bool readInSeqDB(string);
	
	//string locSearch(const DBProtein&) const;
	
	//properities
	bool writeOut(string, const FilterFileParams&) const;
	bool writeOutDB(string, const FilterFileParams&);
	
	//modifiers
	bool readIn(string, FilterFileParams&, Peptides* const);
	void addSubcelluarLoc();
};

/*************/
/* functions */
/*************/

/* dtafilter.cpp */
bool isColumnHeaderLine(const vector<string>&);
string parseSample(string, string, bool);
int parsePeptideSC(string);
string parseReplicate(string);
string getID(string);

#endif /* DTarray_AJM_hpp */

