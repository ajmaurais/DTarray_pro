//
//  FilterFile.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "FilterFile.hpp"

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
						outputFormat = utils::toInt(param.value);
						includeProteins = (outputFormat == 1 || outputFormat == 2 || outputFormat == 3);
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
						peptideOutput = utils::toInt(param.value);
						includePeptides = (peptideOutput == 1 || peptideOutput == 2 || peptideOutput == 3);
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
	
	while(!data.end())
	{
		line = data.getLine_skip_trim();
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
		else if (elems.size() != 0)
			return false;
	}
	
	return true;
}

inline string FilterFileParams::getFileColname(int index) const
{
	return file[index].colname;
}

FilterFileData::FilterFileData(string _colName)
{
	colname = _colName;
	count = "0";
	uniquePeptides = "0";
}
