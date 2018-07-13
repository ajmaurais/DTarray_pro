//
//  locReport.hpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 5/31/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include <dtarray_pro.hpp>
#include <params.hpp>
#include <utils.hpp>

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
		std::string sampleName;
		
	public:
		Loc(std::string _sampleName){
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
		void initializeCol(std::string, unsigned int, unsigned int, unsigned int, unsigned int);
		
		inline std::string getSampleName() const{
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
		std::string loc;
		std::string matchDir;
		std::vector<Loc> col;
		
	public:
		static size_t colSize;
		static params::Params* pars;
		
		LocDat(std::string _loc, std::string _matchDir, size_t _colSize, params::Params* _pars){
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
		
		std::string getKey() const{
			return utils::toLower(loc);
		}
		inline bool operator == (const LocDat& comp) const{
			return getKey() == comp.getKey();
		}
		inline void operator = (const LocDat& add){
			col = add.col;
			loc = add.loc;
		}
		
		inline void initializeCol(size_t colIndex, std::string _sampleName,
								  unsigned int _count, unsigned int _specSum,
								  unsigned int _uniqSpecSum, unsigned int _seqCount) {
			col[colIndex].initializeCol(_sampleName, _count, _specSum, _uniqSpecSum, _seqCount);
		}
		void consolidate(const LocDat&);
		
		void write(std::ofstream&, int) const;
		void writeWide(std::ofstream&) const;
		void writeLong(std::ofstream&) const;
	};
	
	class LocDB{
	private:
		//hashTable::HashTable<LocDat>* locTable;
		typedef std::map<std::string, LocDat> LocTableType;
		LocTableType locTable;
		
	public:
		LocDB(){
			//locTable = new hashTable::HashTable<LocDat>(LOC_TABLE_SIZE);
		}
		~LocDB(){
			//delete locTable;
		}
		
		void addLoc(LocDat newLoc);
		void writeLocReport(std::ofstream& outF, int fxnNum);
	};
}

/* locReport_hpp */
