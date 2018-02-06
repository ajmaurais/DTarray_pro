//
//  saintOutput.hpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 1/11/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#ifndef saintOutput_hpp
#define saintOutput_hpp

#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include <utils.hpp>

namespace saint{

	class BaitFile;

	class BaitFile{
	private:
		//hashTable::HashTable<BaitFileData>* dat;
		typedef std::map<std::string, std::string> DatType;
		DatType dat;
		std::string fname;
		
	public:
		//constructor
		BaitFile(std::string _fname){
			fname = _fname;
		}
		~BaitFile(){}
		
		//modifiers
		bool read();
		
		//properties
		std::string getBaitName(std::string) const;
	};
}

#endif /* saintOutput_hpp */
