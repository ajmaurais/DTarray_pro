//
//  locReport.hpp
//  DTarray_pro
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
		void writeLocReport(std::ofstream& outF, int fxnNum) const;
	};
}

/* locReport_hpp */
