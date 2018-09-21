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

#pragma once

#include <map>
#include <baseClasses.hpp>
#include <utils.hpp>

std::string const DAT_NOT_FOUND = "NOT_FOUND_IN_DB";
size_t const DB_DBASE_SIZE = 10000;

class DBProtein;
class Dbase;

class DBProtein : public ProteinTemplate {
private:
	std::string dat;
	
	//modifers
	void clear();
public:
	//constructor
	DBProtein(std::string);
	DBProtein();
	
	//modifers
	void operator = (const DBProtein&);
	
	//properties
	std::string getDat() const{
		return dat;
	}
};

class Dbase{
private:
	typedef std::map<std::string, DBProtein> DbType;
	DbType db;
	
public:
	
	Dbase() {}
	~Dbase() {}
	
	bool readIn(std::string);
	std::string getDat(std::string) const;
};

/* dbase_hpp */
