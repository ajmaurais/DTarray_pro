//
//  params.hpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef params_hpp
#define params_hpp

#include <cassert>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <ctime>

#include "../lib/utils.hpp"

using namespace std;

namespace params{
	
	/******************************/
	/* namespace scoped constants */
	/******************************/
	
	string const PARAM_ERROR_MESSAGE = " is an invalid arguement for: ";
	const string INVALID_ARG = " is an invalid option! Exiting...\nUse DTarray -h for help.";
	enum OutputFormat {none, wideFormat, longFormat, both};
	enum PeptideGroupFormat {byScan, byProtein, byCharge};
	const string VNUM_STR = "<versionNum>";
	const string END_VNUM_STR = "</versionNum>";
	const string DTAFILTER_NAME = "DTASelect-filter.txt";
	
	//default file locations
	const string PROG_WD_HOME = string(getenv("HOME")) + "/scripts/DTarray_pro";
	const string PROG_WD_DB = PROG_WD_HOME + "/db";
	const string LOC_DB_FNAME = PROG_WD_DB + "/humanLoc.tsv";
	const string SEQ_DB_FNAME = PROG_WD_DB + "/humanProteome.fasta";
	const string MW_DB_FNAME = SEQ_DB_FNAME;
	const string AA_DB_FNAME = PROG_WD_DB + "/aaMasses.txt";
	const string FXN_DB_FNAME = PROG_WD_DB + "/humanFxn.tsv";
	const string STATIC_MOD_FNAME = PROG_WD_DB + "/staticModifications.txt";
	const string HELP_FILE_FNAME = PROG_WD_DB + "/helpFile.man";
	const string USAGE_FNAME = PROG_WD_DB + "/usage.txt";
	
	//default file names
	const string DEFAULT_FLIST_NAME = "dtarray_pro_flist.txt";
	const string DEFAULT_SMOD_NAME = "staticModifications.txt";
	const string DTAFILTER_EXT = ".dtafilter";
	const string OFNAME = "DTarray_pro.tsv";
	const string DB_OFNAME = "DTarray_long.tsv";
	const string PEPTIDE_OFNAME = "peptideList.tsv";
	const string PEPTIDE_DB_OFNAME = "peptideList_long.tsv";
	const string SAINT_PREY_FILE = "prey_file.txt";
	const string SAINT_INTERACTION_FILE = "interaction_file.txt";
	const string LOC_TABLE_FNAME = "loc_summary.tsv";
	const string LOC_TABLE_LONG_FNAME = "loc_summary_long.tsv";
	
	/**********************/
	/* class definitions */
	/*********************/
	
	class Params;
	class FilterFileParam;
	
	class FilterFileParam{
	private:
		string path;
		string colname;
	public:
		//constructor
		FilterFileParam(string);
		~FilterFileParam() {}
		
		//properities
		string getPath() const{
			return path;
		}
		string getColname() const{
			return colname;
		}
	};
		
	//stores names and locations for DTA filter files and output paramaters found in
	//params file.
	class Params{
	private:
		friend class Proteins;
		vector<FilterFileParam> file;
		
		string parseVersionNum(string) const;
		
		void displayHelp() const;
		void usage() const;
		bool writeStdFlist(ofstream&) const;
		bool writeSubdirFlist(ofstream&) const;
		void purgeDir(string) const;
		
		string wd;
		int numFiles;
		
	public:
		string flistName;
		string versionNum;
		string inputFormat;
		static OutputFormat outputFormat;
		string sampleNamePrefix;
		bool parseSampleName;
		bool includeUnique;
		bool getSubCelluarLoc;
		bool rewriteFlist;
		bool rewriteSmod;
		string locDBfname;
		bool calcMW;
		string aaDBfanme, mwDBFname, staticModsFname;
		string ofname, dbOfname, dbPeptideOfFname, peptideOfFname;
		bool getSeq;
		string seqDBfname;
		bool includeCoverage;
		bool includeSequenceCount;
		bool includePeptides;
		bool includeProteins;
		static OutputFormat peptideOutput;
		bool getFxn;
		string fxnDBfname;
		bool useDefaultSeqDB;
		bool includeNullPeptides;
		int supInfoOutput;
		unsigned int supInfoNum, peptideSupInfoNum;
		PeptideGroupFormat peptideGroupMethod;
		bool includeSaint;
		string saintBaitFile, saintPreyFname, saintInteractionFname;
		bool includeReverse;
		bool wdSpecified;
		int modGroupMethod;
		bool includeModStat;
		string locTableFname, locTableLongFname;
		static OutputFormat locOutput;
		int locSupInfoNum;
		
		Params ()
		{
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
			rewriteSmod = false;
			locDBfname = LOC_DB_FNAME;
			calcMW = false;
			aaDBfanme = AA_DB_FNAME;
			mwDBFname = MW_DB_FNAME;
			staticModsFname = DEFAULT_SMOD_NAME;
			ofname = OFNAME;
			dbOfname = DB_OFNAME;
			dbPeptideOfFname = PEPTIDE_DB_OFNAME;
			peptideOfFname = PEPTIDE_OFNAME;
			getSeq = false;
			seqDBfname = SEQ_DB_FNAME;
			includeCoverage = false;
			includeSequenceCount = false;
			includePeptides = false;
			includeProteins = true;
			getFxn = false;
			fxnDBfname = FXN_DB_FNAME;
			useDefaultSeqDB = true;
			includeNullPeptides = false;
			supInfoOutput = 0;
			peptideGroupMethod = byProtein;
			supInfoNum = 0;
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
		}
		
		//modifiers
		bool readFlist(string, string);
		bool getOpts(int,const char* const []);
		
		//properties
		bool writeFlist();
		bool optionsCompatable() const;
		string getFilePath(size_t index) const {
			return file[index].getPath();
		}
		string getFileColname(size_t index) const {
			return file[index].getColname();
		}
		bool writeSmod(string) const;
		
		int getNumFiles() const{
			return numFiles;
		}
		string getwd() const{
			return wd;
		}
	};
	
	OutputFormat Params::outputFormat = wideFormat;
	OutputFormat Params::peptideOutput = none;
	OutputFormat Params::locOutput = none;
	
	OutputFormat intToOutputFormat(int);
	PeptideGroupFormat intToGroupFormat(int);
	string groupFormatString(PeptideGroupFormat);
}

#endif /* params_hpp */
