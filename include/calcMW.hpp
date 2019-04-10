//
//  calcMW.hpp
//  DTarray_pro
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

#pragma once

#include <map>

#include <params.hpp>
#include <molecularFormula.hpp>

namespace mwDB{
	class SeqDB;
	
	std::string const SEQ_NOT_FOUND = "SEQUENCE_NOT_FOUND_IN_DB";
	
	class SeqDB{
		typedef std::map<std::string, std::string> seqLibraryType;
		seqLibraryType seqLibrary;
	public:
		SeqDB(){}
		~SeqDB(){}
		
		//modifers
		bool readIn(std::string);
		std::string getSequence(std::string) const;
	};
	
	class MWDB_Protein : public utils::Residues{
	public:
		SeqDB* seqDB;
		
		MWDB_Protein() : utils::Residues(){
			seqDB = new SeqDB;
		}
		~MWDB_Protein(){
			delete seqDB;
		}
		
		bool initialize(const params::Params&);
	};
	
	std::string getID(std::string str);
}

/* calcMW_h */
