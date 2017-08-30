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
#include "params.hpp"
#include "calcMW.hpp"
#include "saintOutput.hpp"
#include "locReport.hpp"

using namespace std;

/******************************/
/* globally scoped constants */
/*****************************/

bool const INCLUDE_FULL_DESCRIPTION = true;
string const DEFAULT_COL_NAMES [] = {"Full_description", "ID", "Protein", "Description", "pI", "Length(aa)", "Mass(Da)"};
size_t const DEFAULT_COL_NAMES_LENGTH = 7;
string const SEQ_NOT_FOUND = "SEQUENCE_NOT_FOUND_IN_DB";
size_t const PROTEINS_DATA_SIZE = 500;
size_t const PEPTIDES_DATA_SIZE = 2500;
string const DEFAULT_COL_NAMES_DB [] = {"Full_description", "ID", "Protein", "Description", "pI", "Length(aa)", "Mass(Da)",
	"Long_sample_name", "Spectral_counts"};
size_t const DEFAULT_COL_NAMES_DB_LENGTH = 9;
string const PARSE_SAMPLE_NAME_HEADERS [] = {"Sample", "Replicate"};
size_t const PARSE_SAMPLE_NAME_HEADERS_LEN = 2;
string const SUP_INFO_HEADERS[] = {"SC", "Unique_pep_SC", "Coverage", "Sequence_count", "Num_mod_pep", "SC_mod_pep"};
size_t const SUP_INFO_HEADERS_LEN = 6;
string const PEP_SUP_INFO_HEADERS[] = {"SC", "Mod_pep_SC"};
size_t const PEP_SUP_INFO_HEADERS_LEN = 2;
string const MWCALC_HEADERS [] = {"Avg_mass", "Monoisotopic_mass", "Sequence"};
size_t const MWCALC_HEADERS_LENGTH = 2;
string const DEFALUT_PEPTIDE_COLNAMES [] = {"Protein_ID", "Parent_protein", "Protein_description", "Sequence", "Length(aa)", "Unique", "CalcMH"};
string const DEFALUT_PEPTIDE_DB_COLNAMES [] = {"Protein_ID", "Parent_protein", "Protein_description", "Sequence", "Length(aa)", "Unique",
	"CalcMH", "Long_sample_name", "Spectral_counts", "Sample", "Replicate"};
size_t const DEFALUT_PEPTIDE_DB_COLNAMES_LEN = 9;
string const REVERSE_MATCH = "Reverse_";
const char* DIFFMODS = "*";

char const DB_DELIM = ';';
string const LOC_REPORT_HEADERS [] = {"Count", "Sum_SC", "Sum_seq_count"};

/**********************/
/* class definitions */
/*********************/

class Protein;
class Proteins;
class Peptide;
class Peptides;

class Peptide : public ProteinDataTemplate<SampleData_peptide> {
	friend class Proteins;
	friend class Peptides;
public:
	Peptide (params::Params* const par, mwDB::MWDB* const _mwdb) : ProteinDataTemplate <SampleData_peptide>(par)
	{
		mwdb = _mwdb;
	}
	Peptide () : ProteinDataTemplate <SampleData_peptide> () {}
	~Peptide() {}
	
	//modifers
	inline void clear();
	void calcMW();
	inline void operator = (const Peptide&);
	
	//properities
	inline bool operator == (const Peptide& comp) const{
		return comp.key == key;
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
	inline void parseSequence(const string&);
};

mwDB::MWDB* Peptide::mwdb = nullptr;

class Peptides : public DBTemplate<Peptide> {
	friend class Proteins;
private:
	mwDB::MWDB* mwdb;
public:
	Peptides(const params::Params& pars) : DBTemplate<Peptide>(pars, PEPTIDES_DATA_SIZE) {}
	Peptides() : DBTemplate<Peptide>(){}
	~Peptides(){}
	
	//properities
	bool writeOut(string, const params::Params&);
	bool writeOutDB(string, const params::Params&);
};

//stores data for each protein found in filter file
class Protein : public ProteinTemplate , public ProteinDataTemplate<SampleData_protein> {
	friend class Proteins;
private:
	string MW, loc, fxn;
	string fullDescription, pI;
	
	//pointers to Proteins data
	static Dbase* locDB;
	static mwDB::MWDB_Protein* mwdb;
	static mwDB::SeqDB* seqDB;
	static Dbase* fxnDB;
	static saint::BaitFile* baitFile;
	static locReport::LocDB* locTable;
	
	//modifier
	void getProteinData(string);
	inline void getProtein(string);
	inline void clear();
	void calcMW();
	void addSeq();
	void addLoc();
	void addFxn();
	void addSupData();
	void addLocToTable();
	
	void writeCount(ofstream&) const;
	void writeUnique(ofstream&) const;
	void writeCoverage(ofstream&) const;
	void writeSequenceCount(ofstream&) const;
	void writeModStat(ofstream&) const;
	
	void writeProtein(ofstream&); //write fxn 0
	void writePrey(ofstream&) const; //write fxn 1
	void writeInteractions(ofstream&) const; //write fxn 2
	
public:
	Protein(params::Params* const pars,
			Dbase* const _locDB,
			Dbase* const _fxnDB,
			mwDB::MWDB_Protein* const _mwdb,
			mwDB::SeqDB* const _seqDB,
			saint::BaitFile* const _baitFile,
			locReport::LocDB* const _locTable)
		: ProteinDataTemplate<SampleData_protein>(pars) {
		locDB = _locDB;
		mwdb = _mwdb;
		seqDB = _seqDB;
		fxnDB = _fxnDB;
		baitFile = _baitFile;
		locTable = _locTable;
	}
	Protein() : ProteinDataTemplate<SampleData_protein>() {}
	~Protein(){}
	
	inline void operator = (const Protein&);
	
	void consolidate(const Protein&);
	void write(ofstream&, int);
	void apply(int);
	locReport::LocDat toLocDat(const string&) const;
};

Dbase* Protein::locDB = nullptr;
Dbase* Protein::fxnDB = nullptr;
mwDB::MWDB_Protein* Protein::mwdb = nullptr;
mwDB::SeqDB* Protein::seqDB = nullptr;
saint::BaitFile* Protein::baitFile = nullptr;
locReport::LocDB* Protein::locTable = nullptr;

//stores data for all proteins found in DTA filter files
class Proteins : public DBTemplate<Protein>{
	friend class saint::BaitFile;
	Dbase* locDB;
	Dbase* fxnDB;
	mwDB::MWDB_Protein* mwdb;
	mwDB::SeqDB* seqDB;
	saint::BaitFile* baitFile;
	locReport::LocDB* locTable;
	
	//modifers
	bool readIn(params::Params* const,
				const vector <SampleData_protein>&,
				const vector<SampleData_peptide>&,
				Peptides* const);
public:
	//constructor
	Proteins(const params::Params& pars) : DBTemplate<Protein>(pars, PROTEINS_DATA_SIZE){
		locDB = nullptr;
		seqDB = nullptr;
		mwdb = nullptr;
		fxnDB = nullptr;
		baitFile = nullptr;
		locTable = nullptr;
	}
	Proteins() : DBTemplate<Protein>(){
		locDB = nullptr;
		seqDB = nullptr;
		mwdb = nullptr;
		fxnDB = nullptr;
		baitFile = nullptr;
		locTable = nullptr;
	}
	~Proteins(){
		delete locDB;
		delete mwdb;
		delete seqDB;
		delete fxnDB;
		delete baitFile;
		delete locTable;
	}
	
	bool readInLocDB(string);
	bool readInMWdb(string, const params::Params&);
	bool readInSeqDB(string);
	bool readInFxnDB(string);
	bool readBaitFile(string);
	void buildLocTable();
	
	//properities
	bool writeOut(string, const params::Params&);
	bool writeOutDB(string, const params::Params&);
	bool writeSaint(string, int) const;
	bool writeWideLocTable(string, const params::Params&) const;
	bool writeLongLocTable(string, const params::Params&) const;
	
	//modifiers
	bool readIn(params::Params* const, Peptides* const);
};

/*************/
/* functions */
/*************/

string parseSample(string, string, bool, bool);
int parsePeptideSC(string);
int parseModPeptide(string);
string parseReplicate(string);
string getID(string);
inline string parseSequence(string);

#endif /* dtafilter_hpp */

