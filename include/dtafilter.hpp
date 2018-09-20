//
//  dtafilter.hpp
//  DTarray_AJM
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

#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <string>
#include <algorithm>
#include <stdexcept>

#include <dtarray_pro.hpp>
#include <baseClasses.hpp>
#include <dbase.hpp>
#include <utils.hpp>
#include <params.hpp>
#include <calcMW.hpp>
#include <saintOutput.hpp>
#include <locReport.hpp>
#include <molecularFormula.hpp>

/******************************/
/* globally scoped constants */
/*****************************/

bool const INCLUDE_FULL_DESCRIPTION = true;
std::string const DEFAULT_COL_NAMES [] = {"Full_description", "ID", "Protein", "Description", "pI",
	"Length(aa)", "Mass(Da)"};
size_t const DEFAULT_COL_NAMES_LENGTH = 7;
size_t const PROTEINS_DATA_SIZE = 500;
size_t const PEPTIDES_DATA_SIZE = 2500;
std::string const DEFAULT_COL_NAMES_DB [] = {"Full_description", "ID", "Protein", "Description",
	"pI", "Length(aa)", "Mass(Da)",	"Long_sample_name", "Spectral_counts"};
size_t const DEFAULT_COL_NAMES_DB_LENGTH = 9;
std::string const PARSE_SAMPLE_NAME_HEADERS [] = {"Sample", "Replicate"};
size_t const PARSE_SAMPLE_NAME_HEADERS_LEN = 2;
std::string const SUP_INFO_HEADERS[] = {"SC", "Unique_pep_SC", "Coverage", "Sequence_count",
	"Num_mod_pep", "SC_mod_pep"};
size_t const SUP_INFO_HEADERS_LEN = 6;
std::string const PEP_SUP_INFO_HEADERS[] = {"SC", "Mod_pep_SC"};
size_t const PEP_SUP_INFO_HEADERS_LEN = 2;
std::string const MWCALC_HEADERS [] = {"Avg_mass", "Monoisotopic_mass", "Formula", "Sequence"};
size_t const MWCALC_HEADERS_LENGTH = 3;
std::string const DEFALUT_PEPTIDE_COLNAMES [] = {"Protein_ID", "Parent_protein", "Protein_description",
	"Flanking_sequence", "Sequence", "Length(aa)", "Unique", "CalcMH"};
std::string const DEFALUT_PEPTIDE_DB_COLNAMES [] = {"Protein_ID", "Parent_protein", "Protein_description", "Sequence", "Length(aa)", "Unique", "CalcMH", "Long_sample_name", "Spectral_counts", "Sample",
	"Replicate"};
size_t const DEFALUT_PEPTIDE_DB_COLNAMES_LEN = 9;

char const DB_DELIM = ';';
std::string const LOC_REPORT_HEADERS [] = {"Count", "Sum_SC", "Sum_seq_count"};

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
	Peptide (params::Params* const par, molFormula::Residues* const _mwdb) : ProteinDataTemplate <SampleData_peptide>(par){
		mwdb = _mwdb;
	}
	Peptide () : ProteinDataTemplate <SampleData_peptide> () {}
	~Peptide() {}
	
	//modifers
	void clear();
	void calcMW();
	void operator = (const Peptide&);
	
	//properities
	bool operator == (const Peptide& comp) const{
		return comp.key == key;
	}
	std::string makeKey() const;
	void consolidate(const Peptide&);
	void write(std::ofstream&);
	
private:
	std::string key, calcSequence;
	std::string proteinID, calcMH, fileName, protein, description, charge;
	bool unique;
	
	static molFormula::Residues* mwdb;

	void parsePeptide(const std::string&);
	void parseSequence(const std::string&);
};

class Peptides : public DBTemplate<Peptide> {
	friend class Proteins;
public:
	Peptides(const params::Params& pars) : DBTemplate<Peptide>(pars, PEPTIDES_DATA_SIZE) {}
	Peptides() : DBTemplate<Peptide>(){}
	~Peptides(){}
	
	//properities
	bool writeOut(std::string, const params::Params&);
	bool writeOutDB(std::string, const params::Params&);
};

//stores data for each protein found in filter file
class Protein : public ProteinTemplate , public ProteinDataTemplate<SampleData_protein> {
	friend class Proteins;
private:
	std::string MW, loc, fxn;
	std::string fullDescription, pI;
	
	//pointers to Proteins data
	static Dbase* locDB;
	static mwDB::MWDB_Protein* mwdb;
	static mwDB::SeqDB* seqDB;
	static Dbase* fxnDB;
	static saint::BaitFile* baitFile;
	static locReport::LocDB* locTable;
	
	//modifier
	bool getProteinData(std::string);
	bool parse_matchDir_ID_Protein(std::string);
	void clear();
	void calcMW();
	void addSeq();
	void addLoc();
	void addFxn();
	void addSupData();
	void addLocToTable();
	
	void writeCount(std::ofstream&) const;
	void writeUnique(std::ofstream&) const;
	void writeCoverage(std::ofstream&) const;
	void writeSequenceCount(std::ofstream&) const;
	void writeModStat(std::ofstream&) const;
	
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
	
	void operator = (const Protein&);
	
	void consolidate(const Protein&);
	//void apply(int);
	void writeProtein(std::ofstream&);
	void writePrey(std::ofstream&) const;
	void writeInteractions(std::ofstream&) const;
	locReport::LocDat toLocDat(const std::string&) const;
};

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
				const std::vector <SampleData_protein>&,
				const std::vector<SampleData_peptide>&,
				Peptides* const);
public:
	enum OutputFiles {preyFile, interactionFile};
	
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
	
	bool readInLocDB(std::string);
	bool readInMWdb(const params::Params&);
	bool readInSeqDB(std::string);
	bool readInFxnDB(std::string);
	bool readBaitFile(std::string);
	void buildLocTable();
	
	//properities
	bool writeOut(std::string, const params::Params&);
	bool writeOutDB(std::string, const params::Params&);
	bool writeSaint(std::string, OutputFiles) const;
	bool writeWideLocTable(std::string, const params::Params&) const;
	bool writeLongLocTable(std::string, const params::Params&) const;
	
	//modifiers
	bool readIn(params::Params* const, Peptides* const);
};

/*************/
/* functions */
/*************/

int parsePeptideSC(std::string);
int parseModPeptide(std::string);
//std::string getID(std::string);
//inline std::string parseSequence(std::string);

/* dtafilter_hpp */

