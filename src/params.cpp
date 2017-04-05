//
//  params.cpp
//  DTarray_Pro
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "params.hpp"

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
				utils::systemCommand("man " + HELP_FILE_FNAME);
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
				inputFormat = string(argv[i]);
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
					cerr << argv[i] << PARAM_ERROR_MESSAGE << "outputFormat" << endl;
					return false;
				}
				outputFormat = intToOutputFormat(utils::toInt(argv[i]));
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
					cerr << "Specified direectory does not exist." << endl;
					return false;
				}
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
				else if(!strcmp(argv[i], "smod"))
					rewriteSmod = true;
				else{
					cerr << argv[i] << PARAM_ERROR_MESSAGE << "rewrite" << endl;
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
					useDefaultSeqDB = true;
				}
				else if(utils::isArg(argv[i+1]))
					mwDBFname = utils::absPath(argv[++i]);
				else throw runtime_error("bad opts!");
				calcMW = true;
				continue;
			}
			if(!strcmp(argv[i], "-seq"))
			{
				if(!utils::isArg(argv[i+1]))
					seqDBfname = SEQ_DB_FNAME;
				else if(utils::isArg(argv[i+1]))
					seqDBfname = utils::absPath(argv[++i]);
				else throw runtime_error("bad opts!");
				getSeq = true;
				continue;
			}
			if(!strcmp(argv[i], "-smod"))
			{
				if(!writeSmod(wd))
					cerr << "Could not write new smod file!" << endl;
				return false;
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
					cerr << argv[i] << PARAM_ERROR_MESSAGE << "peptideOutput" << endl;
					return false;
				}
				peptideOutput = intToOutputFormat(utils::toInt(argv[i]));
				includePeptides = (peptideOutput != none);
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
					cerr << argv[i] << PARAM_ERROR_MESSAGE << "peptideGroupMethod" << endl;
					return false;
				}
				peptideGroupMethod = intToGroupFormat(utils::toInt(argv[i]));
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
					cerr << argv[i] << PARAM_ERROR_MESSAGE << "modGroupMethod" << endl;
					return false;
				}
				modGroupMethod = utils::toInt(argv[i]);
				continue;
			}
			if(!strcmp(argv[i], "-modS"))
			{
				includeModStat = true;
				supInfoNum += 2;
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
					cerr << argv[i] << PARAM_ERROR_MESSAGE << "supInfoOutput" << endl;
					return false;
				}
				supInfoOutput = utils::toInt(argv[i]);
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
				saintBaitFile = utils::absPath(argv[i]);
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
					cerr << argv[i] << PARAM_ERROR_MESSAGE << "includeReverse" << endl;
					return false;
				}
				includeReverse = utils::toInt(argv[i]);
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
					sampleNamePrefix = string(argv[++i]);
				}
				else throw runtime_error("bad opts!");
				continue;
			}
			if(!strcmp(argv[i], "-v") || !strcmp(argv[i], "--version"))
			{
				cerr << "DTarray_pro " << BIN_VERSION_NUM << endl;
				cerr << "Last git commit: " << GIT_DATE << endl;
				cerr << "git revision: " << GIT_COMMIT << endl;
				return false;
			}
			if(!strcmp(argv[i], "--purge"))
			{
				purgeDir(wd);
				return false;
			}
			if(!strcmp(argv[i], "--pswd"))
			{
				cerr << PROG_WD_HOME << endl;
				return false;
			}
			if(!strcmp(argv[i], "--oswd"))
			{
				utils::systemCommand("open " + PROG_WD_HOME);
				return false;
			}
			else{
				cerr << argv[i] << INVALID_ARG << endl;
				usage();
				return false;
			}
		}
		
		//fix options
		if(wd[wd.length() - 1] != '/')
			wd += "/";
		if(calcMW && !getSeq)
			seqDBfname = mwDBFname;
		
		return true;
	}
	
	bool Params::writeSmod(string _wd) const
	{
		if(_wd[_wd.length() - 1] != '/')
			_wd += "/";
		ofstream outF((_wd + DEFAULT_SMOD_NAME).c_str());
		utils::File staticMods(STATIC_MOD_FNAME);
		if(!outF || !staticMods.read())
			return false;
		
		if(wdSpecified)
			cerr << "Generating " << _wd << DEFAULT_SMOD_NAME << endl;
		else cerr << "Generating ./" << DEFAULT_SMOD_NAME << endl;
		
		outF << utils::COMMENT_SYMBOL << " Static modifications for DTarray_pro" << endl
		<< utils::COMMENT_SYMBOL << " File generated on: " << utils::ascTime() << endl
		<< "<staticModifications>" << endl;
		
		while(!staticMods.end())
			outF << staticMods.getLine() << endl;
		
		outF << endl << "</staticModifications>" << endl;
		
		return true;
	}
	
	void Params::usage() const
	{
		utils::systemCommand("cat " + USAGE_FNAME);
		cerr << endl;
	}
	
	//removes all DTarray_pro generated files with default flie names in
	//working dirrectory
	void Params::purgeDir(string _wd) const
	{
		if(utils::dirExists(_wd))
		{
			if(_wd[_wd.length() - 1] != '/')
				_wd += "/";
			
			string deleteFiles [] = {DEFAULT_FLIST_NAME, DEFAULT_SMOD_NAME,
				OFNAME, DB_OFNAME, PEPTIDE_OFNAME, PEPTIDE_DB_OFNAME,
				SAINT_PREY_FILE, SAINT_INTERACTION_FILE};
			
			for(string* p = utils::begin(deleteFiles); p != utils::end(deleteFiles); ++p)
			{
				if(utils::fileExists(_wd + *p))
				{
					if(wdSpecified)
						cerr << "Removed " << _wd << *p << endl;
					else cerr << "Removed ./" << *p << endl;
					utils::systemCommand("rm -f " + _wd + *p);
				}
			}
		}
		else{
			cerr << "Dir does not exist!" << endl;
		}
	}
	
	bool Params::writeFlist()
	{
		if(wd[wd.length() - 1] != '/')
			wd += "/";
		ofstream outF((wd + flistName).c_str());
		if(!outF)
			return false;
		
		outF << utils::COMMENT_SYMBOL << "File list for DTarray_pro" << endl
		<< utils::COMMENT_SYMBOL << "File List generated on: " << utils::ascTime() << endl;
		outF << endl << VNUM_STR << BIN_VERSION_NUM << END_VNUM_STR << "\n<flist>\n\n";
		
		if(inputFormat == "standard")
			return writeStdFlist(outF);
		else if(inputFormat == "subdir")
			return writeSubdirFlist(outF);
		else return false;
	}
	
	bool Params::writeStdFlist(ofstream& outF) const
	{
		assert(outF);
		vector<string> filterFiles;
		if(!utils::ls(wd.c_str(), filterFiles, DTAFILTER_EXT))
		{
			cerr << "\nDTA-filter files could not be found in the specified directory! Exiting..." << endl;
			return false;
		}
		
		for(vector<string>::iterator it = filterFiles.begin(); it != filterFiles.end(); ++it)
			outF << (*it).substr(0, (*it).length() - DTAFILTER_EXT.length()) << OUT_DELIM << *it << endl;
		
		outF << "\n</flist>\n";
		
		return true;
	}
	
	bool Params::writeSubdirFlist(ofstream& outF) const
	{
		assert(outF);
		vector<string> files;
		vector<string> filterFiles;
		if(!utils::ls(wd.c_str(), files))
			return false;
		
		for(vector<string>::iterator it = files.begin(); it != files.end(); ++it)
			if(utils::dirExists(wd + *it))
				if(utils::fileExists(wd + *it + "/DTASelect-filter.txt"))
					filterFiles.push_back(*it);
		
		for(vector<string>::iterator it = filterFiles.begin(); it != filterFiles.end(); ++it)
			outF << *it << OUT_DELIM << *it << "/DTASelect-filter.txt" << endl;
		
		outF << "\n</flist>\n";
		
		return true;
	}

	FilterFileParam::FilterFileParam(string line)
	{
		vector<string>elems;
		utils::split(line, '\t', elems);
		assert(elems.size() == 2);
		colname = elems[0];
		path = elems[1];
	}
	
	bool Params::readFlist(string fname, string path)
	{
		utils::File data;
		if(!data.read(path + fname))
			return false;

		numFiles = 0;
		string line;
		
		do{
			line = data.getLine_skip_trim();
			if(line.substr(0, VNUM_STR.length()) == VNUM_STR) //check line begins with VNUM_STR
			{
				versionNum = parseVersionNum(line);
				if(!(MIN_BIN_VERSION_NUM <= utils::toDouble(versionNum)))
				{
					cerr << "File list was generated under binary version: " << versionNum
						<< endl << "Min flist version is: " << MIN_BIN_VERSION_NUM << endl;
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
	
	string Params::parseVersionNum(string line) const
	{
		size_t before = line.find(VNUM_STR);
		size_t end = line.find(END_VNUM_STR);
		if(end == string::npos || before != 0)
			return "-1";
		
		size_t len = line.length() - (VNUM_STR.length() + END_VNUM_STR.length());
		string mid = line.substr(before + VNUM_STR.length(), len);
		
		return mid;
	}
	
	bool Params::optionsCompatable() const
	{
		bool good = true;
		assert(MIN_BIN_VERSION_NUM <= utils::toDouble(versionNum));
		if(peptideGroupMethod == byScan && ((peptideOutput == wideFormat) || (peptideOutput == both)))
		{
			cerr << endl << "peptideGroupMethod and peptideOutput options are incompatable!" << endl
				<< "Use DTarray -h for more info." << endl << endl;
			good = false;
		}
		if(modGroupMethod == 1 && (peptideGroupMethod == byScan || peptideOutput == none))
		{
			cerr << endl << "modGroupMethod and peptide output options are incompatable!" << endl
				<< "Use DTarray -h for more info." << endl << endl;
			good = false;
		}
		if(supInfoOutput == 1 && supInfoNum <= 0)
		{
			cerr << endl <<"Non zero supInfoOutput with zero supInfoNum." << endl
				<< "Use DTarray -h for more info." << endl << endl;
			good = false;
		}
		return good;
	}
}
