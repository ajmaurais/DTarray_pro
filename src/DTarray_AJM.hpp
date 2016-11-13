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
#include "baseClasses.hpp"
#include "dbase.hpp"
#include "utils.hpp"
#include "FilterFile.hpp"
#include "hashTable.cpp"
#include "calcMW.hpp"
#include <stdexcept>

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
size_t const PROTEINS_DATA_SIZE = 500;
size_t const PEPTIDES_DATA_SIZE = 2500;

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

/**********************/
/* class definitions */
/*********************/

class Protein;
class Proteins;
class Peptide;

/* #################### subcelluarLoc.cpp #################### */

class Peptide : public ProteinDataTemplate<FilterFileData_peptide> {
	friend class Proteins;
	friend class Peptides;
public:
	Peptide (FilterFileParams * par, mwDB::MWDB* const _mwdb) : ProteinDataTemplate <FilterFileData_peptide> (par) {
		mwdb = _mwdb;
	}
	Peptide () : ProteinDataTemplate <FilterFileData_peptide> () {}
	~Peptide() {}
	
	//modifers
	inline void clear();
	void calcMW();
	inline void operator = (const Peptide&);
	
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
	void write(ofstream&);
	
private:
	string key, calcSequence;
	string proteinID, calcMH, scan, protein, description, charge;
	bool unique;
	
	static mwDB::MWDB* mwdb;

	void parsePeptide(const string&);
	void addSupData(mwDB::MWDB_Protein* const);
};

mwDB::MWDB* Peptide::mwdb = nullptr;

class Peptides : public DBTemplate<Peptide> {
	friend class Proteins;
private:
	mwDB::MWDB* mwdb;
public:
	Peptides(const FilterFileParams& pars) : DBTemplate<Peptide>(pars, PEPTIDES_DATA_SIZE) {}
	Peptides() : DBTemplate<Peptide>(){}
	~Peptides(){}
	
	//properities
	bool writeOut(string, const FilterFileParams&);
	bool writeOutDB(string, const FilterFileParams&);
};

//stores data for each protein found in filter file
class Protein : public ProteinTemplate , public ProteinDataTemplate<FilterFileData_protein> {
	friend class Proteins;
	friend class hashTable::HashTable<Protein>;
	friend class hashTable::LinkedList<Protein>;
private:
	string MW, loc, fxn;
	string fullDescription, matchDirrection;
	int sequenceCount;
	
	//pointers to Proteins data
	static Dbase* locDB;
	static mwDB::MWDB_Protein* mwdb;
	static mwDB::SeqDB* seqDB;
	static Dbase* fxnDB;
	
	//modifier
	bool getProteinData(string, size_t);
	inline void getProteinAndDescr(string);
	DBProtein toDBprotein() const;
	inline void clear();
	
public:
	Protein(FilterFileParams* const pars, Dbase* const _locDB, Dbase* const _fxnDB, mwDB::MWDB_Protein* const _mwdb, mwDB::SeqDB* const _seqDB) : ProteinDataTemplate<FilterFileData_protein>(pars) {
		locDB = _locDB;
		mwdb = _mwdb;
		seqDB = _seqDB;
		fxnDB = _fxnDB;
	}
	Protein() : ProteinDataTemplate<FilterFileData_protein>() {}
	~Protein(){}
	
	inline void operator = (const Protein&);
	
	void consolidate(const Protein&);
	void calcMW();
	void addSeq();
	void addLoc();
	void addFxn();
	void write(ofstream&);
};

Dbase* Protein::locDB = nullptr;
Dbase* Protein::fxnDB = nullptr;
mwDB::MWDB_Protein* Protein::mwdb = nullptr;
mwDB::SeqDB* Protein::seqDB = nullptr;

//stores data for all proteins found in DTA filter files
class Proteins : public DBTemplate<Protein>{
	Dbase* locDB;
	Dbase* fxnDB;
	mwDB::MWDB_Protein* mwdb;
	mwDB::SeqDB* seqDB;
	
	//modifers
	bool readIn(string, FilterFileParams* const, const FilterFileParam&,
				const vector <FilterFileData_protein>&, const vector<FilterFileData_peptide>&,
				Peptides* const);
public:
	//constructor
	Proteins(const FilterFileParams& pars) : DBTemplate<Protein>(pars, PROTEINS_DATA_SIZE){
		locDB = nullptr;
		seqDB = nullptr;
		mwdb = nullptr;
		fxnDB = nullptr;
	}
	Proteins() : DBTemplate<Protein>(){
		locDB = nullptr;
		seqDB = nullptr;
		mwdb = nullptr;
		fxnDB = nullptr;
	}
	~Proteins(){
		delete locDB;
		delete mwdb;
		delete seqDB;
		delete fxnDB;
	}
	
	bool readInLocDB(string);
	bool readInMWdb(string, const FilterFileParams&);
	bool readInSeqDB(string);
	bool readInFxnDB(string);
	
	//properities
	bool writeOut(string, const FilterFileParams&);
	bool writeOutDB(string, const FilterFileParams&);
	
	//modifiers
	bool readIn(string, FilterFileParams&, Peptides* const);
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
inline string parseSequence(string);

#endif /* DTarray_AJM_hpp */

