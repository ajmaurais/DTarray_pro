//
//  saintOutput.cpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 1/11/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#include "saintOutput.hpp"

namespace saint{
	
	BaitFileData::BaitFileData(string line)
	{
		vector<string> elems;
		utils::split(line, '\t', elems);
		assert(elems.size() == 3);
		ipName = elems[0];
		baitName = elems[1];
		tc = elems[2];
	}
	
	bool BaitFile::read()
	{
		utils::File file(fname);
		if(!file.read())
			return false;
		
		string line;
		
		while(!file.end())
		{
			line = file.getLine();
			BaitFileData newBaitDat(line);
			dat->insert(newBaitDat, newBaitDat.getIPname());
		}
		return true;
	}
	
	string BaitFile::getBaitName(string _key) const
	{
		BaitFileData* item = dat->getItem(_key);
		if(item == nullptr)
			return "BAIT_NOT_FOUND";
		return item->getBaitName();
	}
}
