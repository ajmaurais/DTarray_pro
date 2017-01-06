//
//  FilterFile.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "FilterFile.hpp"

namespace filterFile{
	
	OutputFormat intToOutputFormat(int val)
	{
		switch(val){
			case 0 : return none;
				break;
			case 1 : return wideFormat;
				break;
			case 2 : return longFormat;
				break;
			case 3 : return both;
				break;
			default: throw runtime_error("Invalid type!");
		}
	}
	
	PeptideGroupFormat intToGroupFormat(int val)
	{
		switch(val){
			case 0 : return byScan;
				break;
			case 1 : return byProtein;
				break;
			case 2 : return byCharge;
				break;
			default: throw runtime_error("Invalid type!");
		}
	}
	
	string groupFormatString(PeptideGroupFormat format)
	{
		switch(format){
			case byScan : return "by scan";
				break;
			case byProtein : return "by protein";
				break;
			case byCharge : return "by protein and charge";
				break;
		}
	}

	//read in files to combine and output paramaters from params file
	bool FilterFileParams::readDTParams(string fname, string path)
	{
		utils::File file(path + fname);
		if(!file.read(path + fname))
			return false;
		string line;
		
		do{
			line = file.getLine_skip_trim();
			if(line == "<params>")
				do{
					line = file.getLine_skip_trim();
					if(utils::strContains('=', line)) //find lines containing params by = symbol
					{
						Param param (line);
						if(param.param == "versionNum")
						{
							versionNum = param.value;
							assert(MIN_BIN_VERSION_NUM <= utils::toDouble(versionNum));
							continue;
						}
						if(param.param == "sampleNamePrefix")
						{
							sampleNamePrefix = param.value;
							continue;
						}
						if(param.param == "outputFormat")
						{
							if(!(param.value == "0" || param.value == "1" || param.value == "2" || param.value == "3"))
							{
								cout << param.value << PARAM_ERROR_MESSAGE << "outputFormat" << endl;
								return false;
							}
							outputFormat = intToOutputFormat(utils::toInt(param.value));
							includeProteins = (outputFormat != none);
							continue;
						}
						if(param.param == "locDBfname")
						{
							locDBfname = param.value;
							continue;
						}
						if(param.param == "includeUnique")
						{
							assert(param.value == "0" || param.value == "1");
							includeUnique = utils::toInt(param.value);
							continue;
						}
						if(param.param == "getSubCelluarLoc")
						{
							assert(param.value == "0" || param.value == "1");
							getSubCelluarLoc = utils::toInt(param.value);
							continue;
						}
						if(param.param == "calcMW")
						{
							assert(param.value == "0" || param.value == "1");
							calcMW = utils::toInt(param.value);
							continue;
						}
						if(param.param == "mwDBFname")
						{
							mwDBFname = param.value;
							continue;
						}
						if(param.param == "aaDBfanme")
						{
							aaDBfanme = param.value;
							continue;
						}
						if(param.param == "staticModsFname")
						{
							staticModsFname = param.value;
							continue;
						}
						if(param.param == "ofname")
						{
							ofname = param.value;
							continue;
						}
						if(param.param == "dbOfname")
						{
							dbOfname = param.value;
							continue;
						}
						if(param.param == "peptideOfFname")
						{
							peptideOfFname = param.value;
							continue;
						}
						if(param.param == "dbPeptideOfFname")
						{
							dbPeptideOfFname = param.value;
							continue;
						}
						if(param.param == "seqDBfname")
						{
							seqDBfname = param.value;
							continue;
						}
						if(param.param == "getSeq")
						{
							assert(param.value == "0" || param.value == "1");
							getSeq = utils::toInt(param.value);
							continue;
						}
						if(param.param == "includeCoverage")
						{
							assert(param.value == "0" || param.value == "1");
							includeCoverage = utils::toInt(param.value);
							continue;
						}
						if(param.param == "includePeptides")
						{
							if(!(param.value == "0" || param.value == "1" || param.value == "2" || param.value == "3"))
							{
								cout << param.value << PARAM_ERROR_MESSAGE << "includePeptides" << endl;
								return false;
							}
							peptideOutput = intToOutputFormat(utils::toInt(param.value));
							includePeptides = (peptideOutput != 0);
							continue;
						}
						if(param.param == "getFxn")
						{
							assert(param.value == "0" || param.value == "1");
							getFxn = utils::toInt(param.value);
							continue;
						}
						if(param.param == "fxnDBfname")
						{
							fxnDBfname = param.value;
							continue;
						}
						if(param.param == "useDefaultSeqDB")
						{
							assert(param.value == "0" || param.value == "1");
							useDefaultSeqDB = utils::toInt(param.value);
							continue;
						}
						if(param.param == "includeNullPeptides")
						{
							assert(param.value == "0" || param.value == "1");
							includeNullPeptides = utils::toInt(param.value);
							continue;
						}
						if(param.param == "supInfoOutput")
						{
							if(!(param.value == "0" || param.value == "1"))
							{
								cout << param.value << PARAM_ERROR_MESSAGE << "supInfoOutput" << endl;
								return false;
							}
							supInfoOutput = utils::toInt(param.value);
							continue;
						}
						if(param.param == "groupPeptides")
						{
							if(!(param.value == "0" || param.value == "1" || param.value == "2"))
							{
								cout << param.value << PARAM_ERROR_MESSAGE << "groupPeptides" << endl;
								return false;
							}
							peptideGroupMethod = intToGroupFormat(utils::toInt(param.value));
							continue;
						}
						else return false;
					}
					else if(line != "</params>")
						return false;
				} while(line != "</params>");
			
		} while(!file.end() && line != "</paramsFile>");
		
		if(calcMW && !getSeq)
			seqDBfname = mwDBFname;
		
		return true;
	}
	
	bool FilterFileParams::readFlist(string fname, string path)
	{
		utils::File data;
		if(!data.read(path + fname))
			return false;
		
		int i = 0;
		numFiles = 0;
		string line;
		string versionNumStr = "<versionNum>";
		
		do{
			line = data.getLine_skip_trim();
			if(line.substr(0, versionNumStr.length()) == versionNumStr)
			{
				assert(MIN_BIN_VERSION_NUM <= utils::toDouble(parseVersionNum(line)));
				continue;
			}
			if(line == "<flist>")
			{
				do{
					line = data.getLine_skip_trim();
					if(line == "</flist>")
						continue;
					vector<string>elems;
					utils::split(line, '\t', elems);
					if(elems.size() == 2)
					{
						FilterFileParam blank;
						file.push_back(blank);
						file[i].colname = elems[0];
						file[i].path = elems[1];
						numFiles++;
						i++;
					}
					else if (elems.size() != 0 || data.end())
						return false;
				} while(line != "</flist>");
				continue;
			}
		} while(!data.end());
		
		return true;
	}
	
	string FilterFileParams::parseVersionNum(string line) const
	{
		string versionNumStr = "<versionNum>";
		string endVersionNumStr = "</versionNum>";
		
		size_t before = line.find(versionNumStr);
		size_t end = line.find(endVersionNumStr);
		if(end == string::npos || before != 0)
			return "-1";
		
		size_t len = line.length() - (versionNumStr.length() + endVersionNumStr.length());
		string mid = line.substr(before + versionNumStr.length(), len);
		
		return mid;
	}
	
	bool FilterFileParams::optionsCompatable() const
	{
		assert(MIN_BIN_VERSION_NUM <= utils::toDouble(versionNum));
		if(BIN_VERSION_NUM != versionNum)
		{
			cout << "DTarray_pro cpp binary file was compiled from version: " << BIN_VERSION_NUM << endl;
			cout << "DTarray_pro.sh is version: " << versionNum << endl;
			cout << "Version incompatability may cause unpredictable behavior.\n\n";
		}
		if(peptideGroupMethod == byScan && ((peptideOutput == wideFormat) || (peptideOutput == both)))
		{
			cout << "groupPeptides and peptideOutput options are incompatable!" << endl
				<< "Use DTarray -h for more info." << endl << endl;
			return false;
		}
		return true;
	}
	
	FilterFileData::FilterFileData(string _colName)
	{
		colname = _colName;
		count = 0;
		uniquePeptides = "0";
	}
}
