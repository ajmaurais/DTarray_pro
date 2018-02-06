//
//  saintOutput.cpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 1/11/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#include <saintOutput.hpp>

namespace saint{
	
	bool BaitFile::read()
	{
		utils::File file(fname);
		if(!file.read())
			return false;
		
		std::string line;
		std::vector<std::string> elems;
		
		while(!file.end())
		{
			line = file.getLine();
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
