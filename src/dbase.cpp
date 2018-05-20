//
//  subCelluarLoc.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include <dbase.hpp>

DBProtein::DBProtein()
{
	ID = "";
	protein = "";
	description = "";
	dat = "";
}

void DBProtein::clear()
{
	ID.clear();
	protein.clear();
	description.clear();
	dat.clear();
}

DBProtein::DBProtein(std::string line)
{
	if(line.empty())
		clear();
	else
	{
		std::vector<std::string>elems;
		utils::split(line, '\t', elems);
		assert(elems.size() == 4);
		
		ID = elems[0];
		protein = elems[1];
		description = elems[2];
		dat = elems[3];
	}
}

void DBProtein::operator = (const DBProtein& pget)
{
	ID = pget.ID;
	protein = pget.protein;
	description = pget.description;
	dat = pget.dat;
}

bool Dbase::readIn(std::string fname)
{
	utils::File file(fname);
	if(!file.read(fname))
		return false;
	
	std::string line;
	
	while(!file.end()){
		line = file.getLine_skip_trim();
		DBProtein newDBProtein(line);
		if (newDBProtein.getID() != "ID")  //skip line if it is header line
		{
			//db->insert(newDBProtein, newDBProtein.getID());
			db[newDBProtein.getID()] = newDBProtein;
		}
	}
	return true;
}

std::string Dbase::getDat(std::string key) const
{
	/*DBProtein* const tempNode = db->getItem(key);
	if(tempNode == nullptr)
		return DAT_NOT_FOUND;
	else return tempNode->getDat();*/
	
	DbType::const_iterator foundIt = db.find(key);
	if(foundIt == db.end())
		return DAT_NOT_FOUND;
	else return foundIt->second.getDat();
}

