//
//  calcMW.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef calcMW_h
#define calcMW_h

#include "hashTable.hpp"
#include "binTree.cpp"

class MWDB;
class nwNode;
class Peptide;
class AATree;
class SeqDB;

namespace mwDB{
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
		SeqDB();
		~SeqDB();
		
		//modifers
		bool readIn(string);
		
		string getSequence(string) const;
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
		
		bool readIn(string, const FilterFileParams&);
	};
	
}

#endif /* calcMW_h */
