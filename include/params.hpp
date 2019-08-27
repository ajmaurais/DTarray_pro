//
//  params.hpp
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

#include <cassert>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>

#include <dtarray_pro.hpp>
#include <utils.hpp>
#include <config.hpp>

namespace params{
	
	/******************************/
	/* namespace scoped constants */
	/******************************/
	
	const std::string  PARAM_ERROR_MESSAGE = " is an invalid argument for: ";
	const std::string INVALID_ARG = " is an invalid option! Exiting...\nUse DTarray -h for help.";
	const std::string VNUM_STR = "<versionNum>";
	const std::string END_VNUM_STR = "</versionNum>";
	const std::string DTAFILTER_NAME = "DTASelect-filter.txt";
	
	//default file locations
    const std::string PROG_WD_HOME = CONFIG_PROG_WD_DIR;
	const std::string PROG_WD_DB = PROG_WD_HOME + "/db";
	const std::string PROG_LOC_DB_FNAME = PROG_WD_DB + "/humanLoc.tsv";
	const std::string PROG_SEQ_DB_FNAME = PROG_WD_DB + "/humanProteome.fasta";
	const std::string PROG_FXN_DB_FNAME = PROG_WD_DB + "/humanFxn.tsv";
	const std::string PROG_ATOM_COUNT_TABLE_FNAME = PROG_WD_DB + "/defaultResidueAtoms.txt";
	const std::string PROG_ATOM_MASS_TABLE_FNAME = PROG_WD_DB + "/atomMasses.txt";
	const std::string PROG_HELP_FILE_FNAME = PROG_WD_DB + "/helpFile.man";
	const std::string PROG_USAGE_FNAME = PROG_WD_DB + "/usage.txt";
	
	//default file names
	const std::string DEFAULT_FLIST_NAME = "dtarray_pro_flist.txt";
	const std::string DEFAULT_ATOM_COUNT_TABLE_FNAME = "atomCountTable.txt";
	const std::string DTAFILTER_EXT = ".dtafilter";
	const std::string OFNAME = "DTarray_pro.tsv";
	const std::string DB_OFNAME = "DTarray_long.tsv";
	const std::string PEPTIDE_OFNAME = "peptideList.tsv";
	const std::string PEPTIDE_DB_OFNAME = "peptideList_long.tsv";
	const std::string SAINT_PREY_FILE = "prey_file.txt";
	const std::string SAINT_INTERACTION_FILE = "interaction_file.txt";
	const std::string LOC_TABLE_FNAME = "loc_report.tsv";
	const std::string LOC_TABLE_LONG_FNAME = "loc_report_long.tsv";
	
	/**********************/
	/* class definitions */
	/*********************/
	
	class Params;
	class FilterFileParam;
	
	class FilterFileParam{
	private:
		std::string path;
		std::string colname;
	public:
		//constructor
		FilterFileParam(std::string);
		~FilterFileParam() {}
		
		//properties
		std::string getPath() const{
			return path;
		}
		std::string getColname() const{
			return colname;
		}
	};
		
	/**
	 stores names and locations for DTA filter files and output parameters found in
	 params file.
	 */
	class Params{
	public:
		enum OutputFormat {none, wideFormat, longFormat, both};
		enum PeptideGroupFormat {byScan, byProtein, byCharge};
		
		std::string flistName;
		std::string versionNum;
		std::string inputFormat;
		static OutputFormat outputFormat;
		std::string sampleNamePrefix;
		bool parseSampleName;
		bool includeUnique;
		bool getSubCelluarLoc;
		bool rewriteFlist;
		std::string locDBfname;
		bool calcMW;
		std::string atomMassTableFname, atomCountTableFname;
		bool unicode;
		std::string ofname, dbOfname, dbPeptideOfFname, peptideOfFname;
		bool printSeq;
		bool getSeq;
		bool includeCoverage;
		bool includeSequenceCount;
		bool getNSAF;
		bool getEMPAI;
		bool includePeptides;
		bool includeProteins;
		static OutputFormat peptideOutput;
		bool getFxn;
		std::string fxnDBfname;
		bool useDefaultSeqDB;
		bool includeNullPeptides;
		int supInfoOutput;
		unsigned int supInfoNum, peptideSupInfoNum;
		PeptideGroupFormat peptideGroupMethod;
		bool includeSaint;
		std::string saintBaitFile, saintPreyFname, saintInteractionFname;
		bool includeReverse;
		bool wdSpecified;
		int modGroupMethod;
		bool includeModStat;
		std::string locTableFname, locTableLongFname;
		static OutputFormat locOutput;
		int locSupInfoNum;
		bool locSummary;
		
		Params ()
		{
			//change default parameters here
			wd = "";
			flistName = DEFAULT_FLIST_NAME;
			versionNum = "";
			inputFormat = "std";
			numFiles = 0;
			sampleNamePrefix = "";
			parseSampleName = false;
			includeUnique = false;
			getSubCelluarLoc = false;
			rewriteFlist = false;
			locDBfname = PROG_LOC_DB_FNAME;
			calcMW = false;
			atomMassTableFname = PROG_ATOM_MASS_TABLE_FNAME;
			atomCountTableFname = PROG_ATOM_COUNT_TABLE_FNAME;
			unicode = false;
			ofname = OFNAME;
			dbOfname = DB_OFNAME;
			dbPeptideOfFname = PEPTIDE_DB_OFNAME;
			peptideOfFname = PEPTIDE_OFNAME;
			printSeq = false;
			getSeq = false;
			seqDBfname = PROG_SEQ_DB_FNAME;
			includeCoverage = false;
			includeSequenceCount = false;
			includePeptides = false;
			includeProteins = true;
			getFxn = false;
			getEMPAI = false;
			getNSAF = false;
			fxnDBfname = PROG_FXN_DB_FNAME;
			useDefaultSeqDB = true;
			includeNullPeptides = false;
			supInfoOutput = 0;
			peptideGroupMethod = byProtein;
			supInfoNum = 0;
			peptideSupInfoNum = 0;
			saintBaitFile = "";
			saintPreyFname = SAINT_PREY_FILE;
			saintInteractionFname = SAINT_INTERACTION_FILE;
			includeSaint = false;
			includeReverse = true;
			wdSpecified = false;
			modGroupMethod = 0;
			includeModStat = false;
			locTableFname = LOC_TABLE_FNAME;
			locTableLongFname = LOC_TABLE_LONG_FNAME;
			locSupInfoNum = 0;
			locCol = "loc";
			locSummary = false;
			excludeStr = "";
			filter = false;
			matchRegex = true;
			toLower = true;
		}
		
		//modifiers
		bool readFlist(std::string, std::string);
		bool getOpts(int, const char* const []);
		
		//properties
		bool writeFlist();
		bool optionsCompatable() const;
		bool writeAtomCountTable(std::string) const;
		void printGitVersion() const;
		std::string getFilePath(size_t index) const {
			return file[index].getPath();
		}
		std::string getFileColname(size_t index) const {
			return file[index].getColname();
		}
		
		size_t getNumFiles() const{
			return numFiles;
		}
		std::string getwd() const{
			return wd;
		}
		std::string getExcludeStr() const{
			return excludeStr;
		}
		std::string getAddStr() const{
			return addStr;
		}
		bool getFilter() const{
			return filter;
		}
		bool getMatchRegex() const{
			return matchRegex;
		}
		bool getToLower() const{
			return toLower;
		}
		std::string getSeqDBfname() const{
			return seqDBfname;
		}
		std::string getLocCol() const;
		static std::string groupFormatString(PeptideGroupFormat);
	private:
		std::vector<FilterFileParam> file;
		
		std::string parseVersionNum(std::string) const;
		OutputFormat intToOutputFormat(int) const;
		PeptideGroupFormat intToGroupFormat(int) const;
		
		void displayHelp() const{
			utils::systemCommand("man " + PROG_HELP_FILE_FNAME);
		}
		void usage(std::ostream& out = std::cerr) const;
		bool writeStdFlist(std::ofstream&) const;
		bool writeSubdirFlist(std::ofstream&) const;
		void purgeDir(std::string) const;
		
		std::string wd;
		std::string seqDBfname;
		std::string excludeStr, addStr;
		bool filter, matchRegex, toLower;
		std::string locCol;
		
		size_t numFiles;
	};
}

/* params_hpp */
