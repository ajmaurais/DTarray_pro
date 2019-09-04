//
//  dtafilter.hpp
//  DTarray_pro
// -----------------------------------------------------------------------------
// Copyright 2018 Aaron Maurais
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
#include <regex>
#include <set>

#include <dtarray_pro.hpp>
#include <baseClasses.hpp>
#include <utils.hpp>
#include <params.hpp>
#include <calcMW.hpp>
#include <saintOutput.hpp>
#include <locReport.hpp>
#include <fastaFile.hpp>
#include <tsvFile.hpp>
#include <molecularFormula.hpp>

/******************************/
/* globally scoped constants */
/*****************************/

bool const INCLUDE_FULL_DESCRIPTION = true;
std::string const DEFAULT_COL_NAMES [] = {"full_description", "id", "protein", "description", "pi",
	"length_aa", "mass_da"};
size_t const DEFAULT_COL_NAMES_LENGTH = 7;
std::string const DEFAULT_COL_NAMES_DB [] = {"full_description", "id", "protein", "description",
	"pi", "length_aa", "mass_da", "long_sample_name", "spectral_counts"};
size_t const DEFAULT_COL_NAMES_DB_LENGTH = 9;
std::string const PARSE_SAMPLE_NAME_HEADERS [] = {"sample", "replicate"};
size_t const PARSE_SAMPLE_NAME_HEADERS_LEN = 2;
std::string const SUP_INFO_HEADERS[] = {"sc", "nsaf", "unique_pep_sc", "coverage", "sequence_count",
	"num_mod_pep", "sc_mod_pep"};
size_t const SUP_INFO_HEADERS_LEN = 6;
std::string const PEP_SUP_INFO_HEADERS[] = {"sc", "mod_pep_sc"};
size_t const PEP_SUP_INFO_HEADERS_LEN = 2;
std::string const MWCALC_HEADERS [] = {"avg_mass", "mono_mass", "formula", "sequence"};
size_t const MWCALC_HEADERS_LENGTH = 3;
std::string const DEFALUT_PEPTIDE_COLNAMES [] = {"protein_id", "parent_protein", "protein_description",
	"full_sequence", "sequence", "length_aa", "unique", "calc_mh"};
std::string const DEFALUT_PEPTIDE_DB_COLNAMES [] = {"protein_id", "parent_protein", "protein_description",
	"full_sequence", "sequence", "length_aa", "unique", "calc_mh", "long_sample_name", "spectral_counts",
	"sample", "replicate"};
size_t const DEFALUT_PEPTIDE_DB_COLNAMES_LEN = 10;

std::string const PEP_SEQ_HEADERS [] = {"begin", "end", "mods"};
size_t const PEP_SEQ_HEADERS_LEN = 3;

char const DB_DELIM = ';';
std::string const LOC_REPORT_HEADERS [] = {"count", "sum_sc", "sum_seq_count"};

std::string const DAT_NOT_FOUND = "NOT_FOUND_IN_DB";

/**********************/
/* class definitions */
/*********************/

class Protein;
class Proteins;
class Peptide;
class Peptides;

class Peptide : public base::ProteinDataTemplate<base::SampleData_peptide> {
	friend class Proteins;
	friend class Peptides;
public:
	Peptide () : ProteinDataTemplate <base::SampleData_peptide>(){}
	~Peptide() {}
	
	//modifers
	void clear();
	void calcMW();
	void calcMod();
	void operator = (const Peptide&);
	
	//properties
	bool operator == (const Peptide& comp) const{
		return comp._key == _key;
	}
	std::string makeKey() const;
	std::string getMods(std::string delim = "|") const;
	void consolidate(const Peptide&);
	void write(std::ostream&);
	
private:
	//!Key used to compare peptides when grouping by charge or mod state
	std::string _key;
	//!Sequence with pre cleavage N and C terminal residues removed
	std::string _calcSequence;
	//!calcSequence with/without modification depending on compare method
	std::string _compareSequence;
	//!calcSequence without modification
	std::string _baseSequence;
	std::string _proteinID, _calcMH, _fileName, _protein, _description, _charge;
	bool unique;
	//!Contains string representations of sites of modifications in peptdie
	std::set<std::string> _mods;
	size_t _begin, _end;
	
	static utils::Residues* _mwdb;
	static utils::FastaFile* _seqDB;
	
	void _parsePeptide(const std::string&);
	void _parseSequence(const std::string&);
};

class Peptides : public base::DBTemplate<Peptide> {
	friend class Proteins;
private:
	utils::Residues _mwdb;
	utils::FastaFile* _seqDB;
	
public:
	Peptides(const params::Params& pars) : DBTemplate<Peptide>(pars){
		_seqDB = new utils::FastaFile(pars.getSeqDBfname());
	}
	Peptides() : DBTemplate<Peptide>(){
		_seqDB = new utils::FastaFile();
	}
	~Peptides(){}
	
	//modifers
	bool readInMWdb(const params::Params&);
	void setMWdb(utils::Residues mwdb){
		_mwdb = mwdb;
	}
	void setSeqDB(utils::FastaFile* const seqdb){
		_seqDB = seqdb;
	}
	
	//properties
	bool writeWide(std::string, const params::Params&);
	bool writeLong(std::string, const params::Params&);
};

//stores data for each protein found in filter file
class Protein : public base::ProteinTemplate, public base::ProteinDataTemplate<base::SampleData_protein> {
	friend class Proteins;
private:
	std::string MW, loc, fxn;
	std::string fullDescription, pI;
	
	//pointers to Proteins data
	static base::StringMap* _locDB;
	static base::StringMap* _fxnDB;
	static utils::Residues* _mwdb;
	static utils::FastaFile* _seqDB;
	static saint::BaitFile* _baitFile;
	static locReport::LocDB* _locTable;
	
	//modifier
	bool getProteinData(std::string);
	bool parse_matchDir_ID_Protein(std::string);
	void clear();
	void calcMW();
	void addSeq();
	void addLoc();
	void addFxn();
	void addNSAF();
	void addSupData();
	void addLocToTable();
	
	void writeCount(std::ostream&) const;
	void writeNSAF(std::ostream&) const;
	void writeUnique(std::ostream&) const;
	void writeCoverage(std::ostream&) const;
	void writeSequenceCount(std::ostream&) const;
	void writeModStat(std::ostream&) const;
	
public:
	Protein() : ProteinDataTemplate<base::SampleData_protein>() {}
	~Protein(){}
	
	void operator = (const Protein&);
	
	void consolidate(const Protein&);
	void writeProtein(std::ostream&);
	void writePrey(std::ostream&) const;
	void writeInteractions(std::ofstream&) const;
	locReport::LocDat toLocDat(const std::string&) const;
};

//stores data for all proteins found in DTA filter files
class Proteins : public base::DBTemplate<Protein>{
	friend class saint::BaitFile;
private:
	base::StringMap _locDB;
	base::StringMap _fxnDB;
	utils::Residues _mwdb;
	utils::FastaFile _seqDB;
	saint::BaitFile _baitFile;
	locReport::LocDB _locTable;
	
	//modifers
	bool readIn(params::Params* const,
				const std::vector <base::SampleData_protein>&,
				const std::vector<base::SampleData_peptide>&,
				Peptides&);
	
	static bool _readDB(std::string fname,
						base::StringMap& dat,
						const std::string columns [2]);
public:
	enum OutputFiles {preyFile, interactionFile};
	
	//constructor
	Proteins(const params::Params& pars) : DBTemplate<Protein>(pars), _seqDB() {}
	Proteins() : DBTemplate<Protein>(){}
	~Proteins(){}
	
	bool readInLocDB(std::string fname, std::string col);
	bool readInMWdb(const params::Params&);
	bool readInSeqDB(std::string);
	bool readInFxnDB(std::string);
	bool readBaitFile(std::string);
	void calcNSAF();
	void calcEMPAI();
	void buildLocTable(bool summary = false);
	
	//properties
	utils::Residues& get_mwdb(){
		return _mwdb;
	}
	utils::FastaFile* get_seqdb(){
		return &_seqDB;
	}
	bool writeWide(std::string, const params::Params&);
	bool writeLong(std::string, const params::Params&);
	bool writeSaint(std::string, OutputFiles) const;
	bool writeWideLocTable(std::string, const params::Params&) const;
	bool writeLongLocTable(std::string, const params::Params&) const;
	
	//modifiers
	bool readIn(params::Params* const, Peptides&);
};

/*************/
/* functions */
/*************/

int parsePeptideSC(std::string);
int parseModPeptide(std::string);

/* dtafilter_hpp */

