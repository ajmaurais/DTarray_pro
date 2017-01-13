//
//  saintOutput.hpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 1/11/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#ifndef saintOutput_hpp
#define saintOutput_hpp

#include<iostream>
#include<vector>
#include "../lib/utils.hpp"
#include "../lib/hashTable.hpp"

using namespace std;

namespace saint{

	class BaitFile;
	class BaitFileData;
	
	class BaitFileData{
	private:
		string ipName, baitName, tc;
	public:
		BaitFileData(string);
		
		//properties
		inline string getBaitName() const{
			return baitName;
		}
		inline string getIPname() const{
			return ipName;
		}
		inline bool operator == (string comp) const{
			return ipName == comp;
		}
	};

	class BaitFile{
	private:
		hashTable::HashTable<BaitFileData>* dat;
		string fname;
		
	public:
		//constructor
		BaitFile(string _fname){
			fname = _fname;
			dat = new hashTable::HashTable<BaitFileData> (10);
		}
		~BaitFile(){
			delete dat;
		}
		
		//modifiers
		bool read();
		
		//properties
		string getBaitName(string) const;
	};
}

#endif /* saintOutput_hpp */
