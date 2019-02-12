//
//  saintOutput.cpp
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

#include <saintOutput.hpp>

namespace saint{
	
	bool BaitFile::read(std::string fname)
	{
		_fname = fname;
		return read();
	}
	
	bool BaitFile::read()
	{
		if(_fname.empty()){
			std::cerr << "\nBaitFile fname not specified!\n";
			return false;
		}
		
		std::ifstream inF(_fname);
		if(!inF) return false;
		
		std::string line;
		std::vector<std::string> elems;
		
		while(utils::safeGetline(inF, line))
		{
			utils::split(line, '\t', elems);
			assert(elems.size() == 3);
			
			std::string ipName = elems[0];
			std::string baitName = elems[1];
			
			dat[ipName] = baitName;
		}
		return true;
	}
	
	std::string BaitFile::getBaitName(std::string _key) const
	{	
		DatType::const_iterator it = dat.find(_key);
		if(it == dat.end())
			return "BAIT_NOT_FOUND";
		else return it->second;
	}
}
