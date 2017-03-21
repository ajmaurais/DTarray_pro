//
//  calcMW.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef calcMW_h
#define calcMW_h

#include "../lib/hashTable.hpp"
#include "../lib/BinTree.hpp"

class MWDB;
class nwNode;
class Peptide;
class AATree;
class SeqDB;

namespace mwDB{
	size_t const SEQ_LIB_SIZE = 10000;
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
		bool operator < (const AminoAcid&) const;
		bool operator > (const AminoAcid&) const;
		bool operator == (const AminoAcid&) const;
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
		
		void printHistogram() const{
			seqLibrary->printHistogram();
		}
	};
	
	class MWDB{
	public:
		binTree::BinTree <AminoAcid>* aminoAcidsDB;
		
		//constructor
		MWDB();
		~MWDB();
		
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
		
		void printHistogram() const{
			seqDB->printHistogram();
		}
		
		bool readIn(string wd, const params::Params&);
	};
	
}

#endif /* calcMW_h */
