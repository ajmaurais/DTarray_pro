//
//  subCelluarLoc.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "dbase.hpp"

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

DBProtein::DBProtein(string line)
{
	if(line.empty())
		clear();
	else
	{
		vector<string>elems;
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

bool Dbase::readIn(string fname)
{
	utils::File file(fname);
	if(!file.read(fname))
		return false;
	
	string line;
	
	while(!file.end()){
		line = file.getLine_skip_trim();
		DBProtein newDBProtein(line);
		if (newDBProtein.getID() != "ID") //skip line if it is header line
			db->insert(newDBProtein, newDBProtein.getID());
	}
	return true;
}

string Dbase::getDat(string key) const
{
	DBProtein* const tempNode = db->getItem(key);
	if(tempNode == nullptr)
		return DAT_NOT_FOUND;
	else return tempNode->getDat();
}

