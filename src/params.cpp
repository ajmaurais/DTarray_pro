//
//  params.cpp
//  DTarray_Pro
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

#include <params.hpp>
#include <gitVersion.hpp>

namespace params{
	
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
			default: throw std::runtime_error("Invalid type!");
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
			default: throw std::runtime_error("Invalid type!");
		}
	}
	
	std::string groupFormatString(PeptideGroupFormat format)
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
	
	OutputFormat Params::outputFormat = wideFormat;
	OutputFormat Params::peptideOutput = none;
	OutputFormat Params::locOutput = none;
	
	bool Params::getOpts(int argc, const char* const argv [])
	{
		//get wd
		wd = utils::pwd();
		assert(utils::dirExists(wd));
		
		//get options
		for(int i = 1; i < argc; i++)
		{
			if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help"))
			{
				displayHelp();
				return false;
			}
			if(!strcmp(argv[i], "-i") || !strcmp(argv[i], "--in"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				rewriteFlist = true;
				inputFormat = std::string(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "-o") || !strcmp(argv[i], "--out"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1") || !strcmp(argv[i], "2") || !strcmp(argv[i], "3")))
				{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "outputFormat" << std::endl;
					return false;
				}
				outputFormat = intToOutputFormat(std::stoi(argv[i]));
				includeProteins = (outputFormat != none);
				continue;
			}
			if(!strcmp(argv[i], "-d") || !strcmp(argv[i], "--dir"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				wdSpecified = true;
				wd = utils::absPath(argv[i]);
				if(!utils::dirExists(wd))
				{
					std::cerr << "Specified direectory does not exist." << std::endl;
					return false;
				}
				continue;
			}
			if(!strcmp(argv[i], "-flist"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				flistName = std::string(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "-rw"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!strcmp(argv[i], "flist"))
					rewriteFlist = true;
				else{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "rewrite" << std::endl;
					return false;
				}
				continue;
			}
			if(!strcmp(argv[i], "-loc"))
			{
				getSubCelluarLoc = true;
				continue;
			}
			if(!strcmp(argv[i], "-fxn"))
			{
				getFxn = true;
				continue;
			}
			if(!strcmp(argv[i], "-mw"))
			{
				if(!utils::isArg(argv[i+1]))
				{
					mwDBFname = seqDBfname;
				}
				else if(utils::isArg(argv[i+1])){
					mwDBFname = std::string(argv[++i]);
					mwDBFnameSpecified = true;
				}
				else throw std::runtime_error("bad opts!");
				calcMW = true;
				continue;
			}
			if(!strcmp(argv[i], "-seq"))
			{
				if(!utils::isArg(argv[i+1])){
					seqDBfname = PROG_SEQ_DB_FNAME;
				}
				else if(utils::isArg(argv[i+1])){
					seqDBfname = argv[++i];
					seqDBFnameSpecified = true;
				}
				else throw std::runtime_error("bad opts!");
				getSeq = true;
				continue;
			}
			if(!strcmp(argv[i], "-mact") || !strcmp(argv[i], "--makeAtomCountTable"))
			{
				if(!writeAtomCountTable(wd))
					std::cerr << "Could not write atomCountTable!" << std::endl;
				return false;
			}
			if(!strcmp(argv[i], "-act") || !strcmp(argv[i], "--atomCountTable"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				atomCountTableFname = std::string(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "--unicode"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
				{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "unicode" << std::endl;
					return false;
				}
				unicode = std::stoi(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "-p") || !strcmp(argv[i], "--peptides"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1") || !strcmp(argv[i], "2") || !strcmp(argv[i], "3")))
				{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "peptideOutput" << std::endl;
					return false;
				}
				peptideOutput = intToOutputFormat(std::stoi(argv[i]));
				includePeptides = (peptideOutput != none);
				continue;
			}
			if(!strcmp(argv[i], "-lr") || !strcmp(argv[i], "--locReport"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1") || !strcmp(argv[i], "2") || !strcmp(argv[i], "3")))
				{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "peptideOutput" << std::endl;
					return false;
				}
				locOutput = intToOutputFormat(std::stoi(argv[i]));
				continue;
			}
			if(!strcmp(argv[i], "-g") || !strcmp(argv[i], "--group"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1") || !strcmp(argv[i], "2")))
				{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "peptideGroupMethod" << std::endl;
					return false;
				}
				peptideGroupMethod = intToGroupFormat(std::stoi(argv[i]));
				continue;
			}
			if(!strcmp(argv[i], "-modG"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
				{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "modGroupMethod" << std::endl;
					return false;
				}
				modGroupMethod = std::stoi(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "-modS"))
			{
				includeModStat = true;
				supInfoNum += 2;
				peptideSupInfoNum++;
				continue;
			}
			if(!strcmp(argv[i], "-u"))
			{
				includeUnique = true;
				supInfoNum++;
				continue;
			}
			if(!strcmp(argv[i], "-c") || !strcmp(argv[i], "--coverage"))
			{
				includeCoverage = true;
				supInfoNum++;
				continue;
			}
			if(!strcmp(argv[i], "-seqC"))
			{
				includeSequenceCount = true;
				supInfoNum++;
				continue;
			}
			if(!strcmp(argv[i], "-s"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
				{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "supInfoOutput" << std::endl;
					return false;
				}
				supInfoOutput = std::stoi(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "-n") || !strcmp(argv[i], "--nullp"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				includeNullPeptides = true;
				continue;
			}
			if(!strcmp(argv[i], "-saint"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				includeSaint = true;
				saintBaitFile = std::string(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "-rev"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
				{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "includeReverse" << std::endl;
					return false;
				}
				includeReverse = std::stoi(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "-f") || !strcmp(argv[i], "--prefix"))
			{
				if(!utils::isArg(argv[i+1]))
				{
					parseSampleName = true;
				}
				else if(utils::isArg(argv[i+1]))
				{
					parseSampleName = true;
					sampleNamePrefix = std::string(argv[++i]);
				}
				else throw std::runtime_error("bad opts!");
				continue;
			}
			if(!strcmp(argv[i], "-v") || !strcmp(argv[i], "--version"))
			{
				std::cerr << "DTarray_pro " << BIN_VERSION_NUM << std::endl;
				std::cerr << "Last git commit: " << GIT_DATE << std::endl;
				std::cerr << "git revision: " << GIT_COMMIT << std::endl;
				return false;
			}
			if(!strcmp(argv[i], "--purge"))
			{
				purgeDir(wd);
				return false;
			}
			if(!strcmp(argv[i], "-pswd"))
			{
				std::cerr << PROG_WD_HOME << std::endl;
				return false;
			}
			if(!strcmp(argv[i], "-oswd"))
			{
				utils::systemCommand("open " + PROG_WD_HOME);
				return false;
			}
			else{
				std::cerr << argv[i] << INVALID_ARG << std::endl;
				usage();
				return false;
			}
		}
		
		//fix options
		if(wd[wd.length() - 1] != '/')
			wd += "/";
		if(calcMW && !getSeq)
			seqDBfname = mwDBFname;
		if(locOutput != none)
		{
			getSubCelluarLoc = true;
			if(includeUnique)
				locSupInfoNum += 1;
		}
		
		return true;
	}
	
	bool Params::writeAtomCountTable(std::string _wd) const
	{
		
		if(_wd[_wd.length() - 1] != '/')
			_wd += "/";
		std::ofstream outF((_wd + DEFAULT_ATOM_COUNT_TABLE_FNAME).c_str());
		std::ifstream inF(PROG_ATOM_COUNT_TABLE_FNAME.c_str());
		if(!inF || !outF)
			return false;
		
		inF.seekg(0, inF.end);
		const long size = inF.tellg();
		inF.seekg(0, inF.beg);
		char* actBuff = (char*) calloc(size, sizeof(char));
		if(!inF.read(actBuff, size))
			return false;
		
		if(wdSpecified)
			std::cerr << std::endl << "Generating " << _wd << DEFAULT_ATOM_COUNT_TABLE_FNAME << std::endl;
		else std::cerr << std::endl <<"Generating ./" << DEFAULT_ATOM_COUNT_TABLE_FNAME << std::endl;
		
		outF << utils::COMMENT_SYMBOL << " Residue atom counts for DTarray_pro" << std::endl
		<< utils::COMMENT_SYMBOL << " File generated on: " << utils::ascTime() << std::endl
		<< std::endl << actBuff << std::endl;
		
		return true;
	}
	
	//print program usage information located in PROG_USAGE_FNAME
	void Params::usage() const
	{
		utils::File file(PROG_USAGE_FNAME);
		assert(file.read());
		
		while(!file.end())
			std::cerr << file.getLine() << std::endl;
	}
	
	//removes all DTarray_pro generated files with default flie names in
	//working dirrectory
	void Params::purgeDir(std::string _wd) const
	{
		if(utils::dirExists(_wd))
		{
			if(_wd[_wd.length() - 1] != '/')
				_wd += "/";
			
			std::string deleteFiles [] = {DEFAULT_FLIST_NAME, DEFAULT_ATOM_COUNT_TABLE_FNAME, OFNAME,
				DB_OFNAME, PEPTIDE_OFNAME, PEPTIDE_DB_OFNAME, SAINT_PREY_FILE,
				SAINT_INTERACTION_FILE, LOC_TABLE_FNAME, LOC_TABLE_LONG_FNAME};
			
			for(std::string* p = utils::begin(deleteFiles); p != utils::end(deleteFiles); ++p)
			{
				if(utils::fileExists(_wd + *p))
				{
					if(wdSpecified)
						std::cerr << "Removed " << _wd << *p << std::endl;
					else std::cerr << "Removed ./" << *p << std::endl;
					utils::systemCommand("rm -f " + _wd + *p);
				}
			}
		}
		else{
			std::cerr << "Dir does not exist!" << std::endl;
		}
	}
	
	bool Params::writeFlist()
	{
		if(wd[wd.length() - 1] != '/')
			wd += "/";
		std::ofstream outF((wd + flistName).c_str());
		if(!outF)
		{
			std::cout << "Could not write flist!" << std::endl;
			return false;
		}
		
		outF << utils::COMMENT_SYMBOL << "File list for DTarray_pro" << std::endl
			<< utils::COMMENT_SYMBOL << "File List generated on: " << utils::ascTime() << std::endl;
			outF << std::endl << VNUM_STR << BIN_VERSION_NUM << END_VNUM_STR << "\n<flist>\n\n";
		
		if(inputFormat == "std")
			return writeStdFlist(outF);
		else if(inputFormat == "subdir")
			return writeSubdirFlist(outF);
		else{
			std::cout << inputFormat << " is not a valid input format!" << std::endl;
			return false;
		}//end of else
	} //end of functon
	
	bool Params::writeStdFlist(std::ofstream& outF) const
	{
		assert(outF);
		std::vector<std::string> filterFiles;
		if(!utils::ls(wd.c_str(), filterFiles, DTAFILTER_EXT))
		{
			std::cerr << "\nDTA-filter files could not be found in the specified directory! Exiting..." << std::endl;
			return false;
		}
		
		for(std::vector<std::string>::iterator it = filterFiles.begin(); it != filterFiles.end(); ++it)
			outF << '\t' << (*it).substr(0, (*it).length() - DTAFILTER_EXT.length()) << OUT_DELIM << *it << std::endl;
		
		outF << "\n</flist>\n";
		
		return true;
	}
	
	bool Params::writeSubdirFlist(std::ofstream& outF) const
	{
		assert(outF);
		std::vector<std::string> files;
		std::vector<std::string> filterFiles;
		if(!utils::ls(wd.c_str(), files))
			return false;
		
		for(std::vector<std::string>::iterator it = files.begin(); it != files.end(); ++it)
			if(utils::dirExists(wd + *it))
				if(utils::fileExists(wd + *it + "/" + DTAFILTER_NAME))
					filterFiles.push_back(*it);
		
		for(std::vector<std::string>::iterator it = filterFiles.begin(); it != filterFiles.end(); ++it)
			outF << *it << OUT_DELIM << *it << "/" << DTAFILTER_NAME << std::endl;
		
		outF << "\n</flist>\n";
		
		return true;
	}

	FilterFileParam::FilterFileParam(std::string line)
	{
		std::vector<std::string>elems;
		utils::split(line, '\t', elems);
		assert(elems.size() == 2);
		colname = elems[0];
		path = elems[1];
	}
	
	bool Params::readFlist(std::string fname, std::string path)
	{
		utils::File data;
		if(!data.read(path + fname))
			return false;

		numFiles = 0;
		std::string line;
		
		do{
			line = data.getLine_skip_trim();
			if(utils::startsWith(line, VNUM_STR)) //check line begins with VNUM_STR
			{
				versionNum = parseVersionNum(line);
				if(!(MIN_BIN_VERSION_NUM <= std::stod(versionNum)))
				{
					std::cerr << "File list was generated under binary version: " << versionNum
						<< std::endl << "Min flist version is: " << MIN_BIN_VERSION_NUM << std::endl;
					return false;
				}
				continue;
			}
			if(line == "<flist>")
			{
				do{
					if(data.end())
						return false;
					line = data.getLine_skip_trim();
					if(line == "</flist>")
						continue;
					file.push_back(FilterFileParam(line));
					numFiles++;
				} while(line != "</flist>");
				continue;
			}
		} while(!data.end());
		
		return true;
	}
	
	std::string Params::parseVersionNum(std::string line) const
	{
		size_t before = line.find(VNUM_STR);
		size_t end = line.find(END_VNUM_STR);
		if(end == std::string::npos || before != 0)
			return "-1";
		
		size_t len = line.length() - (VNUM_STR.length() + END_VNUM_STR.length());
		std::string mid = line.substr(before + VNUM_STR.length(), len);
		
		return mid;
	}
	
	bool Params::optionsCompatable() const
	{
		bool good = true;
		assert(MIN_BIN_VERSION_NUM <= std::stod(versionNum));
		if(peptideGroupMethod == byScan && ((peptideOutput == wideFormat) || (peptideOutput == both)))
		{
			std::cerr << std::endl << "peptideGroupMethod and peptideOutput options are incompatable!" << std::endl
				<< "Use DTarray -h for more info." << std::endl << std::endl;
			good = false;
		}
		if(modGroupMethod == 1 && (peptideGroupMethod == byScan || peptideOutput == none))
		{
			std::cerr << std::endl << "modGroupMethod and peptide output options are incompatable!" << std::endl
				<< "Use DTarray -h for more info." << std::endl << std::endl;
			good = false;
		}
		if(supInfoOutput == 1 && (supInfoNum <= 0 && peptideSupInfoNum <= 0 && locOutput == none))
		{
			std::cerr << std::endl <<"Non zero supInfoOutput with zero supInfoNum." << std::endl
				<< "Use DTarray -h for more info." << std::endl << std::endl;
			good = false;
		}
		if(!getSubCelluarLoc && locOutput != none)
		{
			std::cerr << std::endl << "Zero getSubCelluarLoc with nonzero locOutput." << std::endl
			<< "Use DTarray -h for more info." << std::endl << std::endl;
			good = false;
		}
		return good;
	} //end of function
}//end of namespace
