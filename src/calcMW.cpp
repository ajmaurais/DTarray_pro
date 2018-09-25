//
//  calcMW.cpp
//  DTarray_AJM
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
	
	bool MWDB_Protein::initialize(const params::Params& params)
	{
		bool val1 = seqDB->readIn(params.getmwDBFname());
		bool val2 = molFormula::Residues::initialize(params.atomCountTableFname,
													params.atomMassTableFname);
		
		return val1 && val2;
	}
}
