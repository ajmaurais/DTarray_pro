//
//  subCelluarLoc.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

DBProtein::DBProtein()
{
	ID = "";
	protein = "";
	description = "";
	loc = "";
}

void DBProtein::clear()
{
	ID = "";
	protein = "";
	description = "";
	loc = "";
}

DBProtein::DBProtein(string line)
{
	if(line.empty())
		clear();
	else
	{
		vector<string>elems;
		util::split(line, '\t', elems);
		assert(elems.size() == 4);
		
		ID = elems[0];
		protein = elems[1];
		description = elems[2];
		loc = elems[3];
	}
}

void DBProtein::operator = (const DBProtein& pget)
{
	ID = pget.ID;
	protein = pget.protein;
	description = pget.description;
	loc = pget.loc;
}
