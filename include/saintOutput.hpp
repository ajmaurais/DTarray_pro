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
#include <utils.hpp>
#include <hashTable.hpp>

namespace saint{

	class BaitFile;
	class BaitFileData;
	
	class BaitFileData{
	private:
		std::string ipName, baitName, tc;
	public:
		BaitFileData(std::string);
		
		//properties
		inline std::string getBaitName() const{
			return baitName;
		}
		inline std::string getIPname() const{
			return ipName;
		}
		inline bool operator == (std::string comp) const{
			return ipName == comp;
		}
	};

	class BaitFile{
	private:
		hashTable::HashTable<BaitFileData>* dat;
		std::string fname;
		
	public:
		//constructor
		BaitFile(std::string _fname){
			fname = _fname;
			dat = new hashTable::HashTable<BaitFileData> (10);
		}
		~BaitFile(){
			delete dat;
		}
		
		//modifiers
		bool read();
		
		//properties
		std::string getBaitName(std::string) const;
	};
}

#endif /* saintOutput_hpp */
