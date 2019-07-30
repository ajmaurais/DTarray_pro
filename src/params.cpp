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

namespace params{
	
	Params::OutputFormat Params::intToOutputFormat(int val) const
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
	
	Params::PeptideGroupFormat Params::intToGroupFormat(int val) const
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
	
	std::string Params::groupFormatString(PeptideGroupFormat format)
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
	
	/**
	 Get name of location column.
	 
	 \return one of "subcelluar_loc", "go_cellular_component", "all_locations".
	 */
	std::string Params::getLocCol() const{
		//Options are: subcelluar_loc go_cellular_component all_locations
		if(locCol == "loc")
			return "subcelluar_loc";
		else if(locCol == "go")
			return "go_cellular_component";
		else if(locCol == "both")
			return "all_locations";
		else
			throw std::runtime_error("Unknown column!");
	}
	
	Params::OutputFormat Params::outputFormat = wideFormat;
	Params::OutputFormat Params::peptideOutput = none;
	Params::OutputFormat Params::locOutput = none;
	
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
			if(!strcmp(argv[i], "-lc") || !strcmp(argv[i], "--locCol"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!(!strcmp(argv[i], "loc") || !strcmp(argv[i], "go") || !strcmp(argv[i], "both")))
				{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "locCol" << std::endl;
					return false;
				}
				locCol = std::string(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "-fxn"))
			{
				getFxn = true;
				continue;
			}
			if(!strcmp(argv[i], "-mw"))
			{
				calcMW = true;
				getSeq = true;
				continue;
			}
			if(!strcmp(argv[i], "-seq"))
			{
				printSeq = true;
				getSeq = true;
				continue;
			}
			if(!strcmp(argv[i], "-fasta"))
			{
				if(!utils::isArg(argv[i++])){
					usage();
					return false;
				}
				seqDBfname = utils::absPath(argv[i]);
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
			if(!strcmp(argv[i], "-ls") || !strcmp(argv[i], "--locSummary"))
			{
				locSummary = true;
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
				else return false;
				continue;
			}
			if(!strcmp(argv[i], "-e") || !strcmp(argv[i], "--exclude"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				excludeStr = std::string(argv[i]);
				filter = true;
				continue;
			}
			if(!strcmp(argv[i], "-a") || !strcmp(argv[i], "--add"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				addStr = std::string(argv[i]);
				filter = true;
				continue;
			}
			if(!strcmp(argv[i], "-r") || !strcmp(argv[i], "--regex"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
				{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "matchRegex" << std::endl;
					return false;
				}
				matchRegex = std::stoi(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "-l") || !strcmp(argv[i], "--toLower"))
			{
				if(!utils::isArg(argv[++i]))
				{
					usage();
					return false;
				}
				if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
				{
					std::cerr << argv[i] << PARAM_ERROR_MESSAGE << "toLower" << std::endl;
					return false;
				}
				toLower = std::stoi(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "-v") || !strcmp(argv[i], "--version"))
			{
				printGitVersion();
				return false;
			}
			if(!strcmp(argv[i], "--purge"))
			{
				purgeDir(wd);
				return false;
			}
			if(!strcmp(argv[i], "-pswd"))
			{
				std::cout << PROG_WD_HOME << std::endl;
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
			std::cout << std::endl << "Generating " << _wd << DEFAULT_ATOM_COUNT_TABLE_FNAME << std::endl;
		else std::cout << std::endl <<"Generating ./" << DEFAULT_ATOM_COUNT_TABLE_FNAME << std::endl;
		
		outF << utils::COMMENT_SYMBOL << " Residue atom counts for DTarray_pro" << std::endl
		<< utils::COMMENT_SYMBOL << " File generated on: " << utils::ascTime() << std::endl
		<< std::endl << actBuff << std::endl;
		
		return true;
	}
	
	/**
	 \brief Print usage to \p out
	 
	 Prints contents of ParamsBase::_usageFile to \p out
	 \param out ostream to print to.
	 */
	void Params::usage(std::ostream& out) const
	{
		std::ifstream inF(PROG_USAGE_FNAME);
		std::string line;
		while(utils::safeGetline(inF, line))
			out << line << NEW_LINE;
	}
	
	//removes all DTarray_pro generated files with default file names in
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
			
			for(std::string* p = std::begin(deleteFiles); p != std::end(deleteFiles); ++p)
			{
				if(utils::fileExists(_wd + *p))
				{
					if(wdSpecified)
						std::cout << "Removed " << _wd << *p << std::endl;
					else std::cout << "Removed ./" << *p << std::endl;
					utils::systemCommand("rm -f " + _wd + *p);
				}
			}
		}
		else{
			std::cerr << "Dir does not exist!" << std::endl;
		}
	}
	
	/**
	 Print git version and date and time of last commit to std::out.
	 */
	void params::Params::printGitVersion() const{
		std::cout << "DTarray_pro " << BIN_VERSION_NUM << std::endl;
		std::cout << "Last git commit: " << GIT_VERSION << std::endl;
		std::cout << "git revision: " << GIT_DATE << std::endl;
	}
	
	bool Params::writeFlist()
	{
		if(wd[wd.length() - 1] != '/')
			wd += "/";
		std::ofstream outF((wd + flistName).c_str());
		if(!outF)
		{
			std::cerr << "Could not write flist!" << std::endl;
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
			std::cerr << inputFormat << " is not a valid input format!" << std::endl;
			return false;
		}//end of else
	} //end of function
	
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
		std::ifstream inF(path + fname);
		if(!inF) return false;

		numFiles = 0;
		std::string line;
		
		while(utils::safeGetline(inF, line))
		{
			//line = data.getLine_skip_trim();
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
				std::streampos pos;
				while(utils::safeGetline(inF, line, pos))
				{
					line = utils::trim(line);
					if(line.empty() || utils::isCommentLine(line)) continue;
					if(line == "</flist>"){
						inF.seekg(pos);
						break;
					}
					file.push_back(FilterFileParam(line));
					numFiles++;
				}
				continue;
			}
		}
		
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
			std::cerr << std::endl << "peptideGroupMethod and peptideOutput options are incompatible!" << std::endl
				<< "Use DTarray -h for more info." << std::endl << std::endl;
			good = false;
		}
		if(modGroupMethod == 1 && peptideGroupMethod == byScan)
		{
			std::cerr << std::endl << "modGroupMethod and peptideGroupMethod options are incompatible!" << std::endl
				<< "Use DTarray -h for more info." << std::endl << std::endl;
			good = false;
		}
		if(supInfoOutput == 1 && (supInfoNum <= 0 && peptideSupInfoNum <= 0 && locOutput == none))
		{
			std::cerr << std::endl <<"Non zero supInfoOutput with zero supInfoNum." << std::endl
				<< "Use DTarray -h for more info." << std::endl << std::endl;
			good = false;
		}
		return good;
	} //end of function
}//end of namespace
