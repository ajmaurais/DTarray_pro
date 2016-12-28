//
//  FilterFile.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef FilterFile_hpp
#define FilterFile_hpp

#include <cassert>
#include <vector>
#include <stdexcept>
#include "../lib/utils.hpp"

using namespace std;

namespace filterFile{
	
	string const BLANK_STR = "null";
	string const BLANK_VAL = "-1";
	string const PARAM_ERROR_MESSAGE = " is an invalid arguement for: ";
	enum OutputFormat {none, wideFormat, longFormat, both};
	enum PeptideGroupFormat {NA, byScan, byProtein, byCharge};
	
	/**********************/
	/* class definitions */
	/*********************/
	
	class FilterFileParams;
	class FilterFileParam;
	class FilterFileData;
	
	//stores the data pertaining to a specific filter file (or MS run) for each protein
	class FilterFileData{
	public:
		string colname;
		unsigned int count;
		string uniquePeptides;
		
		//constructor
		FilterFileData (string);
		FilterFileData(){
			colname = BLANK_STR;
			count = 0;
			uniquePeptides = "0";
		}
		~FilterFileData() {}
		
		inline bool isNull() const{
			return count == 0;
		}
	};
	
	class FilterFileData_peptide : public FilterFileData {
	public:
		string parentFile, scan, obsMH;
		
		FilterFileData_peptide(string colName) : FilterFileData(colName)
		{
			parentFile = BLANK_STR;
			scan = BLANK_VAL;
			obsMH = BLANK_VAL;
		}
		FilterFileData_peptide() : FilterFileData() {
			parentFile = BLANK_STR;
			scan = BLANK_VAL;
			obsMH = BLANK_VAL;
		}
		~FilterFileData_peptide() {}
	};
	
	class FilterFileData_protein : public FilterFileData {
	public:
		string coverage;
		
		FilterFileData_protein(string colName) : FilterFileData(colName){
			coverage = "0";
		}
		FilterFileData_protein() : FilterFileData(){
			coverage = "0";
		}
		~FilterFileData_protein() {}
	};
	
	struct FilterFileParam{
		string path;
		string colname;
	};
	
	struct Param {
		string param;
		string value;
		
		//constructor
		Param (string line){
			size_t posStart = line.find("=");
			
			param = line.substr(0, posStart);
			value = line.substr(posStart + 1);
		}
	};
	
	//stores names and locations for DTA filter files and output paramaters found in
	//params file.
	class FilterFileParams{
		friend class Proteins;
		vector<FilterFileParam> file;
	public:
		int numFiles;
		static OutputFormat outputFormat;
		string sampleNamePrefix;
		bool includeUnique;
		bool getSubCelluarLoc;
		string locDBfname;
		bool calcMW;
		string aaDBfanme, mwDBFname, staticModsFname;
		string ofname, dbOfname, dbPeptideOfFname, peptideOfFname;
		bool getSeq;
		string seqDBfname;
		bool includeCoverage;
		bool includePeptides;
		bool includeProteins;
		static OutputFormat peptideOutput;
		bool getFxn;
		string fxnDBfname;
		bool useDefaultSeqDB;
		bool includeNullPeptides;
		int supInfoOutput;
		PeptideGroupFormat peptideGroupMethod;
		
		FilterFileParams ()
		{
			numFiles = 0;
			sampleNamePrefix = "";
			includeUnique = false;
			getSubCelluarLoc = false;
			locDBfname = "";
			calcMW = false;
			aaDBfanme = "";
			mwDBFname = "";
			staticModsFname = "";
			ofname = "";
			dbOfname = "";
			dbPeptideOfFname = "";
			peptideOfFname = "";
			getSeq = false;
			seqDBfname = "";
			includeCoverage = false;
			includePeptides = false;
			includeProteins = false;
			getFxn = false;
			fxnDBfname = "";
			useDefaultSeqDB=false;
			includeNullPeptides = false;
			supInfoOutput = 0;
			peptideGroupMethod = byProtein;
		}
		
		//modifiers
		bool readDTParams(string, string);
		bool readFlist(string, string);
		
		//properties
		string getFilePath(size_t index) const
		{
			return file[index].path;
		}
		string getFileColname(size_t index) const
		{
			return file[index].colname;
		}
		bool optionsCompatable() const;
	};
	
	OutputFormat FilterFileParams::outputFormat = none;
	OutputFormat FilterFileParams::peptideOutput = none;
	
	OutputFormat intToOutputFormat(int);
	PeptideGroupFormat intToGroupFormat(int);
	string groupFormatString(PeptideGroupFormat);
}

#endif /* FilterFile_hpp */
