//
//  dtafilter.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 3/25/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef dtafilter_hpp
#define dtafilter_hpp

#define OUT_DELIM '\t'
#define IN_DELIM '\t'

#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <string>
#include <algorithm>
#include <stdexcept>
#include "baseClasses.hpp"
#include "dbase.hpp"
#include "../lib/utils.hpp"
#include "FilterFile.hpp"
#include "calcMW.hpp"
#include "saintOutput.hpp"

using namespace std;

/******************************/
/* globally scoped constants */
/*****************************/

bool const INCLUDE_FULL_DESCRIPTION = true;
string const DEFAULT_COL_NAMES [] = {"Full_description", "ID", "Protein", "Description", "pI", "length(aa)", "Mass(Da)"};
size_t const DEFAULT_COL_NAMES_LENGTH = 7;
size_t const COLUMN_HEADER_LINE_ELEMENTS_LENGTH = 13;
size_t const MAX_PARAM_ITTERATIONS = 100;
string const SEQ_NOT_FOUND = "SEQUENCE_NOT_FOUND_IN_DB";
size_t const PROTEINS_DATA_SIZE = 500;
size_t const PEPTIDES_DATA_SIZE = 2500;

//editable params for DB output format
string const DEFAULT_COL_NAMES_DB [] = {"Full_description", "ID", "Protein", "Description", "pI", "length(aa)", "Mass(Da)",
	"Long_sample_name", "Spectral_counts"};
size_t const DEFAULT_COL_NAMES_DB_LENGTH = 9;
string const PARSE_SAMPLE_NAME_HEADERS [] = {"Sample", "Replicate"};
size_t const PARSE_SAMPLE_NAME_HEADERS_LEN = 2;
string const SUP_INFO_HEADERS[] = {"SC", "Unique_pep_SC", "coverage", "sequence_count"};
string const MWCALC_HEADERS [] = {"avg_mass", "monoisotopic_mass", "sequence"};
size_t const MWCALC_HEADERS_LENGTH = 2;
string const DEFALUT_PEPTIDE_COLNAMES [] = {"protein_ID", "parent_protein", "protein_description", "sequence", "length(aa)", "unique", "calcMH"};
size_t const DEFALUT_PEPTIDE_COLNAMES_LEN = 7;
string const DEFALUT_PEPTIDE_DB_COLNAMES [] = {"protein_ID", "parent_protein", "protein_description", "sequence", "length(aa)", "unique",
	"calcMH", "Long_sample_name", "Spectral_counts", "Sample", "Replicate"};
size_t const DEFALUT_PEPTIDE_DB_COLNAMES_LEN = 9;

/**********************/
/* class definitions */
/*********************/

class Protein;
class Proteins;
class Peptide;

class Peptide : public ProteinDataTemplate<filterFile::FilterFileData_peptide> {
	friend class Proteins;
	friend class Peptides;
public:
	Peptide (filterFile::FilterFileParams * par, mwDB::MWDB* const _mwdb) : ProteinDataTemplate <filterFile::FilterFileData_peptide> (par) {
		mwdb = _mwdb;
	}
	Peptide () : ProteinDataTemplate <filterFile::FilterFileData_peptide> () {}
	~Peptide() {}
	
	//modifers
	inline void clear();
	void calcMW();
	inline void operator = (const Peptide&);
	
	//properities
	inline bool operator == (const Peptide& comp) const {
		return comp.key == key;
	}
	inline bool operator > (const Peptide& comp) const{
		return key > comp.key;
	}
	inline bool operator < (const Peptide& comp) const{
		return key < comp.key;
	}
	inline string makeKey() const;
	void consolidate(const Peptide&);
	void write(ofstream&, int);
	
private:
	string key, calcSequence;
	string proteinID, calcMH, fileName, protein, description, charge;
	bool unique;
	
	static mwDB::MWDB* mwdb;

	void parsePeptide(const string&);
};

mwDB::MWDB* Peptide::mwdb = nullptr;

class Peptides : public DBTemplate<Peptide> {
	friend class Proteins;
private:
	mwDB::MWDB* mwdb;
public:
	Peptides(const filterFile::FilterFileParams& pars) : DBTemplate<Peptide>(pars, PEPTIDES_DATA_SIZE) {}
	Peptides() : DBTemplate<Peptide>(){}
	~Peptides(){}
	
	//properities
	bool writeOut(string, const filterFile::FilterFileParams&);
	bool writeOutDB(string, const filterFile::FilterFileParams&);
};

//stores data for each protein found in filter file
class Protein : public ProteinTemplate , public ProteinDataTemplate<filterFile::FilterFileData_protein> {
	friend class Proteins;
	friend class hashTable::HashTable<Protein>;
	friend class hashTable::LinkedList<Protein>;
private:
	string MW, loc, fxn;
	string fullDescription, pI;
	
	//pointers to Proteins data
	static Dbase* locDB;
	static mwDB::MWDB_Protein* mwdb;
	static mwDB::SeqDB* seqDB;
	static Dbase* fxnDB;
	static saint::BaitFile* baitFile;
	
	//modifier
	void getProteinData(string);
	inline void getProtein(string);
	DBProtein toDBprotein() const;
	inline void clear();
	void calcMW();
	void addSeq();
	void addLoc();
	void addFxn();
	
	void writeCount(ofstream&) const;
	void writeUnique(ofstream&) const;
	void writeCoverage(ofstream&) const;
	void writeSequenceCount(ofstream&) const;
	
public:
	Protein(filterFile::FilterFileParams* const pars,
			Dbase* const _locDB,
			Dbase* const _fxnDB,
			mwDB::MWDB_Protein* const _mwdb,
			mwDB::SeqDB* const _seqDB,
			saint::BaitFile* const _baitFile)
		: ProteinDataTemplate<filterFile::FilterFileData_protein>(pars) {
		locDB = _locDB;
		mwdb = _mwdb;
		seqDB = _seqDB;
		fxnDB = _fxnDB;
		baitFile = _baitFile;
	}
	Protein() : ProteinDataTemplate<filterFile::FilterFileData_protein>() {}
	~Protein(){}
	
	inline void operator = (const Protein&);
	
	void consolidate(const Protein&);
	void write(ofstream&, int);
	void writeProtein(ofstream&);
	void writePrey(ofstream&) const;
	void writeInteractions(ofstream&) const;
};

Dbase* Protein::locDB = nullptr;
Dbase* Protein::fxnDB = nullptr;
mwDB::MWDB_Protein* Protein::mwdb = nullptr;
mwDB::SeqDB* Protein::seqDB = nullptr;
saint::BaitFile* Protein::baitFile = nullptr;

//stores data for all proteins found in DTA filter files
class Proteins : public DBTemplate<Protein>{
	friend class saint::BaitFile;
	Dbase* locDB;
	Dbase* fxnDB;
	mwDB::MWDB_Protein* mwdb;
	mwDB::SeqDB* seqDB;
	saint::BaitFile* baitFile;
	
	//modifers
	bool readIn(string, filterFile::FilterFileParams* const,
				const vector <filterFile::FilterFileData_protein>&,
				const vector<filterFile::FilterFileData_peptide>&,
				Peptides* const);
public:
	//constructor
	Proteins(const filterFile::FilterFileParams& pars) : DBTemplate<Protein>(pars, PROTEINS_DATA_SIZE){
		locDB = nullptr;
		seqDB = nullptr;
		mwdb = nullptr;
		fxnDB = nullptr;
		baitFile = nullptr;
	}
	Proteins() : DBTemplate<Protein>(){
		locDB = nullptr;
		seqDB = nullptr;
		mwdb = nullptr;
		fxnDB = nullptr;
		baitFile = nullptr;
	}
	~Proteins(){
		delete locDB;
		delete mwdb;
		delete seqDB;
		delete fxnDB;
	}
	
	bool readInLocDB(string);
	bool readInMWdb(string, const filterFile::FilterFileParams&);
	bool readInSeqDB(string);
	bool readInFxnDB(string);
	bool readBaitFile(string);
	
	//properities
	bool writeOut(string, const filterFile::FilterFileParams&);
	bool writeOutDB(string, const filterFile::FilterFileParams&);
	bool writeSaint(string, int) const;
	
	//modifiers
	bool readIn(string, filterFile::FilterFileParams&, Peptides* const);
};

/*************/
/* functions */
/*************/

inline bool isColumnHeaderLine(const string&);
string parseSample(string, string, bool, bool);
int parsePeptideSC(string);
string parseReplicate(string);
string getID(string);
inline string parseSequence(string);

#endif /* dtafilter_hpp */

