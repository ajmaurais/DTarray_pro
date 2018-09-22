//
//  subCelluarLoc.hpp
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

#include <dbase.hpp>

DBProtein::DBProtein()
{
	_ID = "";
	_protein = "";
	_description = "";
	dat = "";
}

void DBProtein::clear()
{
	_ID.clear();
	_protein.clear();
	_description.clear();
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
		
		_ID = elems[0];
		_protein = elems[1];
		_description = elems[2];
		dat = elems[3];
	}
}

void DBProtein::operator = (const DBProtein& pget)
{
	_ID = pget._ID;
	_protein = pget._protein;
	_description = pget._description;
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
		
		//skip line if it is header line
		if (newDBProtein.getID() != "ID")
			db[newDBProtein.getID()] = newDBProtein;
	}
	return true;
}

std::string Dbase::getDat(std::string key) const
{
	DbType::const_iterator foundIt = db.find(key);
	if(foundIt == db.end())
		return DAT_NOT_FOUND;
	else return foundIt->second.getDat();
}

