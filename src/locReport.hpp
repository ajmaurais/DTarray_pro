//
//  locReport.hpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 5/31/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#ifndef locReport_hpp
#define locReport_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include "../lib/hashTable.hpp"
#include "../lib/utils.hpp"

using namespace std;

namespace locReport{
	
	/******************************/
	/* namespace scoped constants */
	/******************************/
	
	int const LOC_TABLE_SIZE = 100;
	
	/**********************/
	/* class definitions */
	/*********************/
	
	class LocDat;
	class Loc;
	class LocDB;
	
	class Loc{
	private:
		unsigned int count, specSum, uniqSpecSum, seqCount;
		string sampleName;
		
	public:
		Loc(string _sampleName){
			sampleName = _sampleName;
			count = 0; specSum = 0; uniqSpecSum = 0; seqCount = 0;
		}
		Loc(){
			count = 0; specSum = 0; uniqSpecSum = 0; seqCount = 0;
		};
		~Loc(){}
		
		inline void operator = (const Loc& add){
			count = add.count; specSum = add.specSum;
			uniqSpecSum = add.uniqSpecSum; seqCount = add.seqCount;
		}
		
		inline void consolidate(const Loc&);
		inline void initializeCol(string, unsigned int, unsigned int, unsigned int, unsigned int);
		
		inline string getSampleName() const{
			return sampleName;
		}
		inline int getCount() const{
			return count;
		}
		inline int getSpecSum() const{
			return specSum;
		}
		inline int getUniqSpecSum() const{
			return uniqSpecSum;
		}
		inline int getSeqCount() const{
			return seqCount;
		}
	};
	
	class LocDat{
		//friend class LocDB;
	private:
		string loc;
		string matchDir;
		vector<Loc> col;
		
	public:
		static size_t colSize;
		static params::Params* pars;
		
		LocDat(string _loc, string _matchDir, size_t _colSize, params::Params* _pars){
			matchDir = _matchDir;
			loc = utils::toLower(_loc);
			col.resize(_colSize);
			colSize = _colSize;
			pars = _pars;
		}
		LocDat(){
			loc = "";
		};
		~LocDat(){}
		
		string getKey() const{
			return utils::toLower(loc);
		}
		inline bool operator == (const LocDat& comp) const{
			return getKey() == comp.getKey();
		}
		inline void operator = (const LocDat& add){
			col = add.col;
			loc = add.loc;
		}
		
		inline void initializeCol(size_t colIndex, string _sampleName,
								  unsigned int _count, unsigned int _specSum,
								  unsigned int _uniqSpecSum, unsigned int _seqCount) {
			col[colIndex].initializeCol(_sampleName, _count, _specSum, _uniqSpecSum, _seqCount);
		}
		void consolidate(const LocDat&);
		
		void write(ofstream&, int) const;
		void writeWide(ofstream&) const;
		void writeLong(ofstream&) const;
	};
	
	size_t LocDat::colSize = 0;
	params::Params* LocDat::pars = nullptr;
	
	class LocDB{
	private:
		hashTable::HashTable<LocDat>* locTable;
		
	public:
		LocDB(){
			locTable = new hashTable::HashTable<LocDat>(LOC_TABLE_SIZE);
		}
		~LocDB(){
			delete locTable;
		}
		
		void addLoc(string, LocDat);
		void writeLocReport(ofstream& outF, int fxnNum){
			locTable->write(outF, fxnNum);
		}
	};
}

#endif /* locReport_hpp */
