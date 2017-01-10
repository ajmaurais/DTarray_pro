//
//  subCelluarLoc.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 11/4/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef dbase_hpp
#define dbase_hpp

#include "baseClasses.hpp"
#include "../lib/hashTable.hpp"
#include "../lib/utils.hpp"

using namespace std;

string const DAT_NOT_FOUND = "NOT_FOUND_IN_DB";

class DBProtein;
class Dbase;

class DBProtein : public ProteinTemplate {
private:
	string dat;
	
	//modifers
	void clear();
public:
	//constructor
	DBProtein(string);
	DBProtein();
	
	//modifers
	void operator = (const DBProtein&);
	
	//properties
	string getDat() const{
		return dat;
	}
};

class Dbase{
private:
	hashTable::HashTable<DBProtein>* db;
	
public:
	
	Dbase() {
		db = new hashTable::HashTable<DBProtein> (10000);
	}
	~Dbase() {
		delete db;
	}
	
	bool readIn(string);
	string getDat(string) const;
};


#endif /* dbase_hpp */
