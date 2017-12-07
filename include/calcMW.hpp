//
//  calcMW.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef calcMW_h
#define calcMW_h

#include <hashTable.hpp>
#include <params.hpp>

namespace mwDB{
	class MWDB;
	class Peptide;
	class SeqDB;
	
	size_t const SEQ_LIB_SIZE = 10000;
	size_t const AA_DB_SIZE = 20;
	size_t const MAX_PARAM_ITTERATIONS = 100;
	
	class Peptide{
		friend class MWDB;
		friend class SeqDB;
	private:
		string ID, sequence;
	public:
		Peptide();
		
		//modifers
		void operator = (const Peptide&);
		
		//properties
		bool operator == (string) const;
		string getID() const;
		string getSequence() const;
	};
	
	class AminoAcid{
		friend class MWDB;
	private:
		string symbol;
		double avgMass, monoMass;
	public:
		//constructor
		AminoAcid(string);
		AminoAcid(string, double, double);
		AminoAcid(string, double);
		AminoAcid();
		
		//modifer
		void operator += (const AminoAcid&);
		
		//properities
		bool operator == (const string& _comp) const{
			return(symbol == _comp);
		}
	};
	
	class SeqDB{
		hashTable::HashTable<Peptide>* seqLibrary;
	public:
		SeqDB(){
			seqLibrary = new hashTable::HashTable<mwDB::Peptide>(SEQ_LIB_SIZE);
		}
		~SeqDB(){
			delete seqLibrary;
		}
		
		//modifers
		bool readIn(string);
		string getSequence(string) const;
	};
	
	class MWDB{
	public:
		hashTable::HashTable <AminoAcid>* aminoAcidsDB;
		
		//constructor
		MWDB(){
			aminoAcidsDB = new hashTable::HashTable<AminoAcid>(AA_DB_SIZE);
		}
		~MWDB(){
			delete aminoAcidsDB;
		}
		
		//modifers
		bool readIn(string, string);
		
		//properties
		double calcMW(string, int) const;
		double getMW(string, int) const;
		double getMW(char, int) const;
		
	protected:
		//modofers
		bool readInAADB(string);
		bool addStaticMod(const AminoAcid&);
	};
	
	class MWDB_Protein : public MWDB{
	public:
		SeqDB* seqDB;
		
		MWDB_Protein() : MWDB(){
			seqDB = new SeqDB;
		}
		~MWDB_Protein(){
			delete seqDB;
		}
		
		bool readIn(string wd, const params::Params&);
	};
}

#endif /* calcMW_h */
