//
//  locReport.cpp
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

#include <locReport.hpp>

namespace locReport{
	
	size_t LocDat::colSize = 0;
	params::Params* LocDat::pars = nullptr;
	
	void Loc::consolidate(const Loc& toAdd)
	{
		count += toAdd.count;
		specSum += toAdd.specSum;
		uniqSpecSum += toAdd.uniqSpecSum;
		seqCount += toAdd.seqCount;
	}
	
	void LocDat::consolidate(const LocDat& toAdd)
	{
		unsigned long colSize = col.size();
		assert(colSize == toAdd.col.size());
		
		for(unsigned long i = 0; i < colSize; i++)
			col[i].consolidate(toAdd.col[i]);
	}
	
	void LocDB::addLoc(LocDat newLoc)
	{
		if(!locsOfInterest.empty())
			if(!std::binary_search(locsOfInterest.begin(), locsOfInterest.end(), newLoc.getKey()))
				return;
		
		LocTableType::iterator it = locTable.find(newLoc.getKey());
		if(it == locTable.end()){
			locTable[newLoc.getKey()] = newLoc;
		}
		else{
			locTable[newLoc.getKey()].consolidate(newLoc);
		}
	}
	
	void Loc::initializeCol(std::string _sampleName, unsigned int _count, unsigned int _specSum,
								   unsigned int _uniqSpecSum, unsigned int _seqCount)
	{
		sampleName = _sampleName;
		count = _count;
		specSum = _specSum;
		uniqSpecSum = _uniqSpecSum;
		seqCount = _seqCount;
	}
	
	void LocDB::writeLocReport(std::ofstream& outF, int fxnNum) const
	{
		for(LocTableType::const_iterator it = locTable.begin(); it != locTable.end(); ++it)
		{
			it->second.write(outF, fxnNum);
		}
	}
	
	void LocDat::writeLong(std::ofstream& outF) const
	{
		if(!outF)
			throw std::runtime_error("Bad std::ofstream");
		
		for(size_t i = 0; i < colSize; i++)
		{
			outF << loc << OUT_DELIM << col[i].getSampleName();
			
			if(pars->parseSampleName)
				outF << OUT_DELIM << parseSample(col[i].getSampleName(),
												 pars->sampleNamePrefix,
												 pars->parseSampleName, 1,
												 pars->getMatchRegex())
				<< OUT_DELIM << parseReplicate(col[i].getSampleName());
			
			outF << OUT_DELIM << col[i].getCount()
			<< OUT_DELIM << col[i].getSpecSum();
			
			if(pars->includeUnique)
				outF << OUT_DELIM << col[i].getUniqSpecSum();
			
			outF << OUT_DELIM << col[i].getSeqCount() << std::endl;
		}
	}
	
	void LocDat::writeWide(std::ofstream& outF) const
	{
		if(!outF)
			throw std::runtime_error("Bad std::ofstream");
		
		assert(pars->locSupInfoNum >= 0 && pars->locSupInfoNum <=1);
		if(pars->supInfoOutput == 0)
		{
			outF << loc;
			for(size_t i = 0; i < colSize; i++)
			{
				outF << OUT_DELIM << col[i].getCount()
				<< OUT_DELIM << col[i].getSpecSum();
				
				if(pars->includeUnique)
					outF << OUT_DELIM << col[i].getUniqSpecSum();
				
				outF << OUT_DELIM << col[i].getSeqCount();
			}
			outF << std::endl;
		}
		else if(pars->supInfoOutput == 1)
		{
			outF << loc;
			
			for(size_t i = 0; i < colSize; i++)
				outF << OUT_DELIM << col[i].getCount();
			
			for(size_t i = 0; i < colSize; i++)
				outF << OUT_DELIM << col[i].getSpecSum();
			
			if(pars->includeUnique)
				for(size_t i = 0; i < colSize; i++)
					outF << OUT_DELIM << col[i].getUniqSpecSum();
			
			for(size_t i = 0; i < colSize; i++)
				outF << OUT_DELIM << col[i].getSeqCount();
			
			outF << std::endl;
		}
	}
	
	void LocDat::write(std::ofstream& outF, int fxnNum) const
	{
		if(!pars->includeReverse)
			if(utils::strContains(REVERSE_MATCH, matchDir))
				return;
		
		switch(fxnNum)
		{
			case 0 : writeWide(outF);
				break;
			case 1 : writeLong(outF);
				break;
			default : throw std::runtime_error("function does not exist!");
		}
	}
	
	/**
	 \return reference to locsOfInterest.
	 */
	const LocDB::LocListType& LocDB::get_locsOfInterset() const{
		return locsOfInterest;
	}
	
	/**
	 Set LocDB::locsOfInterest.
	 A deep copy of \p rhs is made.
	 
	 \param rhs list to copy
	 */
	void LocDB::set_locsOfInterest(const LocListType& rhs){
		locsOfInterest.clear();
		
		for(auto it = rhs.begin(); it != rhs.end(); ++it)
			locsOfInterest.push_back(*it);
			
		std::sort(locsOfInterest.begin(), locsOfInterest.end());
	}

	/**
	 Set LocDB::locsOfInterest to locReport::SUMMARY_LOCS.
	 */
	void LocDB::set_summaryLocs(){
		set_locsOfInterest(LocListType(SUMMARY_LOCS, SUMMARY_LOCS + SUMMARY_LOCS_LEN));
	}
}
