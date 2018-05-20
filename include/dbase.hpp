//
//  subCelluarLoc.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 11/4/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef dbase_hpp
#define dbase_hpp

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
	//hashTable::HashTable<DBProtein>* db;
	typedef std::map<std::string, DBProtein> DbType;
	DbType db;
	
public:
	
	Dbase() {
		//db = new hashTable::HashTable<DBProtein> (DB_DBASE_SIZE);
	}
	~Dbase() {
		//delete db;
	}
	
	bool readIn(std::string);
	std::string getDat(std::string) const;
};


#endif /* dbase_hpp */
