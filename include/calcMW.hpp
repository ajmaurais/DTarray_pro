//
//  calcMW.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef calcMW_h
#define calcMW_h

#include <map>

#include <params.hpp>
#include <molecularFormula.hpp>

namespace mwDB{
	class SeqDB;
	
	/*size_t const SEQ_LIB_SIZE = 10000;
	size_t const AA_DB_SIZE = 20;
	size_t const MAX_PARAM_ITTERATIONS = 100;*/
	
	class SeqDB{
		typedef std::map<std::string, std::string> seqLibraryType;
		seqLibraryType seqLibrary;
	public:
		SeqDB(){}
		~SeqDB(){}
		
		//modifers
		bool readIn(std::string);
		std::string getSequence(std::string) const;
	};
	
	class MWDB_Protein : public molFormula::Residues{
	public:
		SeqDB* seqDB;
		
		MWDB_Protein() : molFormula::Residues(){
			seqDB = new SeqDB;
		}
		~MWDB_Protein(){
			delete seqDB;
		}
		
		bool initalize(std::string wd, const params::Params&);
	};
}

#endif /* calcMW_h */
