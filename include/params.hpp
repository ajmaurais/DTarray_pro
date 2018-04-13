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

#include <utils.hpp>
#include <dtafilter.hpp>

namespace params{
	
	/******************************/
	/* namespace scoped constants */
	/******************************/
	
	const std::string  PARAM_ERROR_MESSAGE = " is an invalid arguement for: ";
	const std::string INVALID_ARG = " is an invalid option! Exiting...\nUse DTarray -h for help.";
	enum OutputFormat {none, wideFormat, longFormat, both};
	enum PeptideGroupFormat {byScan, byProtein, byCharge};
	const std::string VNUM_STR = "<versionNum>";
	const std::string END_VNUM_STR = "</versionNum>";
	const std::string DTAFILTER_NAME = "DTASelect-filter.txt";
	
	//default file locations
	const std::string PROG_WD_HOME = std::string(getenv("HOME")) + "/local/DTarray_pro";
	const std::string PROG_WD_DB = PROG_WD_HOME + "/db";
	const std::string PROG_LOC_DB_FNAME = PROG_WD_DB + "/humanLoc.tsv";
	const std::string PROG_SEQ_DB_FNAME = "/humanProteome.fasta";
	const std::string PROG_MW_DB_FNAME = PROG_SEQ_DB_FNAME;
	//const std::string PROG_AA_DB_FNAME = PROG_WD_DB + "/aaMasses.txt";
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
	const std::string LOC_TABLE_FNAME = "loc_summary.tsv";
	const std::string LOC_TABLE_LONG_FNAME = "loc_summary_long.tsv";
	
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
		
		//properities
		std::string getPath() const{
			return path;
		}
		std::string getColname() const{
			return colname;
		}
	};
		
	//stores names and locations for DTA filter files and output paramaters found in
	//params file.
	class Params{
	private:
		//friend class Proteins;
		std::vector<FilterFileParam> file;
		
		std::string parseVersionNum(std::string) const;
		
		void displayHelp() const{
			utils::systemCommand("man " + PROG_HELP_FILE_FNAME);
		}
		void usage() const;
		bool writeStdFlist(std::ofstream&) const;
		bool writeSubdirFlist(std::ofstream&) const;
		void purgeDir(std::string) const;
		
		std::string wd;
		std::string mwDBFname;
		std::string seqDBfname;
		bool mwDBFnameSpecified, seqDBFnameSpecified;
		
		int numFiles;
		
	public:
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
		bool getSeq;
		bool includeCoverage;
		bool includeSequenceCount;
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
		
		Params ()
		{
			//change default paramaters here
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
			mwDBFname = PROG_MW_DB_FNAME;
			mwDBFnameSpecified = false;
			seqDBFnameSpecified = false;
			atomMassTableFname = PROG_ATOM_MASS_TABLE_FNAME;
			atomCountTableFname = PROG_ATOM_COUNT_TABLE_FNAME;
			unicode = false;
			ofname = OFNAME;
			dbOfname = DB_OFNAME;
			dbPeptideOfFname = PEPTIDE_DB_OFNAME;
			peptideOfFname = PEPTIDE_OFNAME;
			getSeq = false;
			seqDBfname = PROG_SEQ_DB_FNAME;
			includeCoverage = false;
			includeSequenceCount = false;
			includePeptides = false;
			includeProteins = true;
			getFxn = false;
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
		}
		
		//modifiers
		bool readFlist(std::string, std::string);
		bool getOpts(int, const char* const []);
		
		//properties
		bool writeFlist();
		bool optionsCompatable() const;
		bool writeAtomCountTable(std::string) const;
		std::string getFilePath(size_t index) const {
			return file[index].getPath();
		}
		std::string getFileColname(size_t index) const {
			return file[index].getColname();
		}
		
		int getNumFiles() const{
			return numFiles;
		}
		std::string getwd() const{
			return wd;
		}
		std::string getmwDBFname() const{
			/*return (mwDBFnameSpecified ? PROG_WD_DB + mwDBFname :
					getwd() + mwDBFname);*/
			
			std::string ret = (mwDBFnameSpecified ? getwd() + mwDBFname :
							   PROG_WD_DB + mwDBFname);
			return ret;
		}
		std::string getSeqDBfname() const{
			/*return (seqDBFnameSpecified ? PROG_WD_DB + seqDBfname :
					getwd() + seqDBfname);*/
			
			std::string ret = (seqDBFnameSpecified ? getwd() + seqDBfname :
							   PROG_WD_DB + seqDBfname);
			return ret;
		}
	};
	
	OutputFormat Params::outputFormat = wideFormat;
	OutputFormat Params::peptideOutput = none;
	OutputFormat Params::locOutput = none;
	
	OutputFormat intToOutputFormat(int);
	PeptideGroupFormat intToGroupFormat(int);
	std::string groupFormatString(PeptideGroupFormat);
}

#endif /* params_hpp */
