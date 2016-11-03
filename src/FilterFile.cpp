//
//  FilterFile.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "FilterFile.hpp"

Param::Param(string line)
{
	size_t posStart = line.find("=");
	
	param = line.substr(0, posStart);
	value = line.substr(posStart + 1);
}

//read in files to combine and output paramaters from params file
bool FilterFileParams::readDTParams(string fname, string path)
{
	ifstream inF ((path + fname).c_str());
	
	if (!inF)
		return false;
	
	string line;
	
	do{
		util::getLineTrim(inF, line);
		if(util::isCommentLine(line) || line.empty()) //skip line if is comment line
			continue;
		
		if(line == "<params>")
			do{
				util::getLineTrim(inF, line);
				if(util::strContains('=', line)) //find lines containing params by = symbol
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
						outputFormat = util::toInt(param.value);
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
						includeUnique = util::toInt(param.value);
						continue;
					}
					if(param.param == "getSubCelluarLoc")
					{
						assert(param.value == "0" || param.value == "1");
						getSubCelluarLoc = util::toInt(param.value);
						continue;
					}
					if(param.param == "calcMW")
					{
						assert(param.value == "0" || param.value == "1");
						calcMW = util::toInt(param.value);
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
					if(param.param == "includeSeq")
					{
						assert(param.value == "0" || param.value == "1");
						includeSeq = util::toInt(param.value);
						continue;
					}
					if(param.param == "includeCoverage")
					{
						assert(param.value == "0" || param.value == "1");
						includeCoverage = util::toInt(param.value);
						continue;
					}
					if(param.param == "includePeptides")
					{
						if(!(param.value == "0" || param.value == "1" || param.value == "2" || param.value == "3"))
						{
							cout << param.value << PARAM_ERROR_MESSAGE << "includePeptides" << endl;
							return false;
						}
						peptideOutput = util::toInt(param.value);
						includePeptides = (peptideOutput == 1 || peptideOutput == 2 || peptideOutput == 3);
						continue;
					}
					else return false;
				}
				else if(util::isCommentLine(line) || line.empty())
					continue;
				else if(line != "</params>")
					return false;
			} while(line != "</params>");
		
	} while(!inF.eof() && line != "</paramsFile>");
	return true;
}

bool FilterFileParams::readFlist(string fname, string path)
{
	ifstream inF ((path + fname).c_str());
	
	if (!inF)
		return false;
	
	int i = 0;
	numFiles = 0;
	string line;
	
	while(!inF.eof())
	{
		getline(inF, line);
		line = util::trim(line);
		if(util::isCommentLine(line) || line.empty()) //skip line if is comment line
			continue;
		else //else line contains data for filter file
		{
			vector<string>elems;
			util::split(line, '\t', elems);
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
	}
	
	if(numFiles > MAX_NUM_FILES)
		throw runtime_error("\n\nMaxium number of filter files is " + util::toString(MAX_NUM_FILES) + "!\n"
							"You can change the max number of files this program can read by changeing \n"+
							"the value of MAX_NUM_FILES in DTarray_AJM.hpp and recompiling the program\n" +
							"with DTarray --recompile.\n");
	
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
