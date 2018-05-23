//
//  calcMW.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 9/5/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include <calcMW.hpp>

namespace mwDB{
	
	std::string SeqDB::getSequence(std::string id) const
	{
		seqLibraryType::const_iterator it = seqLibrary.find(id);
		if(it == seqLibrary.end())
			return SEQ_NOT_FOUND;
		else return it->second;
	}
	
	std::string getID(std::string str)
	{
		size_t firstBar = str.find("|");
		size_t secBar = str.find("|", firstBar+1);
		return str.substr(firstBar+1, secBar-firstBar-1);
	}
	
	bool SeqDB::readIn(std::string fname)
	{
		utils::File data(fname);
		if(!data.read(fname))
			return false;
		
		std::string line;
		std::string tempID;
		while(!data.end())
		{
			line = data.getLine_skip_trim();
			if(line[0] == '>' && !utils::strContains("Reverse", line))
			{
				tempID = mwDB::getID(line);
				line = data.getLine_skip_trim();
				if(line[0] == '>')
					return false;
				seqLibrary[tempID] = line;
			}//end if
		}//end while
		return true;
	}//end fxn
	
	bool MWDB_Protein::initalize(const params::Params& params)
	{
		bool val1 = seqDB->readIn(params.getmwDBFname());
		bool val2 = molFormula::Residues::initalize(params.atomCountTableFname,
													params.atomMassTableFname);
		
		return val1 && val2;
	}
}
