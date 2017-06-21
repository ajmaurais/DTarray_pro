//
//  locReport.cpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 5/31/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#include "locReport.hpp"

namespace locReport{
	
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
	
	inline void Loc::initializeCol(string _sampleName, unsigned int _count, unsigned int _specSum,
								   unsigned int _uniqSpecSum, unsigned int _seqCount)
	{
		sampleName = _sampleName;
		count = _count;
		specSum = _specSum;
		uniqSpecSum = _uniqSpecSum;
		seqCount = _seqCount;
	}
	
	void LocDat::writeLong(ofstream& outF) const
	{
		if(!outF)
			throw runtime_error("Bad ofstream");
		
		for(size_t i = 0; i < colSize; i++)
		{
			outF << loc << OUT_DELIM << col[i].getSampleName();
			
			if(pars->parseSampleName)
				outF << OUT_DELIM << parseSample(col[i].getSampleName(),
												 pars->sampleNamePrefix,
												 pars->parseSampleName, 1)
				<< OUT_DELIM << parseReplicate(col[i].getSampleName());
			
			outF << OUT_DELIM << col[i].getCount()
			<< OUT_DELIM << col[i].getSpecSum();
			
			if(pars->includeUnique)
				outF << OUT_DELIM << col[i].getUniqSpecSum();
			
			outF << OUT_DELIM << col[i].getSeqCount() << endl;
		}
	}
	
	void LocDat::writeWide(ofstream& outF) const
	{
		if(!outF)
			throw runtime_error("Bad ofstream");
		
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
			outF << endl;
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
			
			outF << endl;
		}
	}
	
	void LocDat::write(ofstream& outF, int fxnNum) const
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
			default : throw runtime_error("functon does not exist!");
		}
	}
}
