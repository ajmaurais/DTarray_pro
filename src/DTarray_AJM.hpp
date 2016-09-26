//
//  DTarray_AJM.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 3/25/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <algorithm>
#include "BinTree.cpp"
#include "hashTable.cpp"

//#define nullptr NULL //to compile on pleiades this line must be included

using namespace std;

/******************************/
/* globally scoped constants */
/*****************************/

/* dtafilter.cpp */
bool const INCLUDE_FULL_DESCRIPTION = true;
//string const DEFAULT_COL_NAMES [] = {"ID", "Description", "Mass (Da)"};
string const DEFAULT_COL_NAMES [] = {"Protein","ID", "Description", "Mass(Da)", "subcellular_location"};
int const DEFAULT_COL_NAMES_LENGTH = 4;
string const COLUMN_HEADER_LINE_ELEMENTS[] = {"Unique", "FileName", "XCorr", "DeltCN", "Conf%", "M+H+",
	"CalcM+H+", "TotalIntensity", "SpR", "ZScore", "IonProportion", "Redundancy", "Sequence"};
int const COLUMN_HEADER_LINE_ELEMENTS_LENGTH = 13;
string const PARAM_ERROR_MESSAGE = " is an invalid arguement for: ";
int const MAX_PARAM_ITTERATIONS = 100;

//editable params for DB output format
string const DEFAULT_COL_NAMES_DB [] = {"Protein","ID", "Description", "Mass(Da)", "Long_sample_name",
	"Spectral_counts", "Sample", "Replicate"};
int const DEFAULT_COL_NAMES_DB_LENGTH = 8;
string const DEFAULT_COL_NAMES_DB_LOC [] = {"Protein","ID", "Description", "Mass(Da)", "subcellular_location", "Long_sample_name",
	"Spectral_counts", "Sample", "Replicate"};
int const DEFAULT_COL_NAMES_DB_LOC_LENGTH = 9;
string const UNIQUE_PEPTIDE_HEADERS[] = {"SC", "Unique_pep_SC"};
string const MWCALC_HEADERS [] = {"sequence", "avg_mass", "monoisotopic_mass"};
int const MWCALC_HEADERS_LENGTH = 3;
bool const INCLUDE_SEQUENCE = true;

/* utils.cpp */
string const WHITESPACE = " \f\n\r\t\v";
string const COMMENT_SYMBOL = "#"; //if changed, paramsCommentSymbol must also be changed in DTarray_AJM.sh

/* hashTable.cpp */
string const SEQ_NOT_FOUND = "SEQUENCE_NOT_FOUND_IN_DB";

/* subcelluarLoc.cpp */
string const LOC_NOT_FOUND = "NOT_FOUND_IN_DB";


/**********************/
/* class definitions */
/*********************/

class FilterFileParams;
class Protein;
class Proteins;
class DBProtein;
class MWDB;
class nwNode;
class Peptide;
class AATree;

/* #################### subCelluarLoc.cpp #################### */

class DBProtein{
	friend class Protein;
	friend class Proteins;
	string ID, gene, description, loc;
	
	//modifers
	void clear();
	
public:
	//constructor
	DBProtein(string);
	DBProtein();
	
	//modifers
	void operator = (const DBProtein&);
	
	//properties
	bool operator == (string) const;
};

/* #################### dtafilter.cpp #################### */

struct FilterFileParam{
	string path;
	string colname;
};

struct Param {
	string param;
	string value;
	
	//constructor
	Param (string);
};

//stores names and locations for DTA filter files and output paramaters found in
//params file.
class FilterFileParams{
	friend class Proteins;
	vector<FilterFileParam> file;
public:
	int numFiles;
	string outputFormat;
	string sampleNamePrefix;
	bool includeUnique;
	bool getSubCelluarLoc;
	string locDBfname;
	bool calcMW;
	string aaDBfanme, peptideDBfname, staticModsFname, ofname;
	
	//modifiers
	bool readDTParams(string, string);
	bool readFlist(string, string);
};

//stores the data pertaining to a specific filter file (or MS run) for each protein
struct FilterFile{
	string colname, count;
	string uniquePeptides;
	
	//constructor
	FilterFile (string, string, string);
};

//stores data for each protein found in filter file
class Protein{
public:
	//properties
	bool operator == (const Protein&) const;
	bool operator > (const Protein&) const;
	bool operator < (const Protein&) const;

private:
	friend class Proteins;
	vector <FilterFile> col;
	
	//modifier
	void initialize(const vector<string>&);
	bool getProteinData(string, int);
	void consolidate(const Protein&, int);
	DBProtein toDBprotein() const;
	void calcMW(const MWDB&);
	void addLoc(string);
	
	string fullDescription, matchDirrection, ID, description, MW, loc;
	double avgMass, monoMass;
	string calcSequence;
};

//stores data for all proteins found in DTA filter files
class Proteins{
	vector <Protein> proteins;
	int colIndex;
	HashTable<DBProtein> locDB;
	
	//modifers
	bool readIn(string, const FilterFileParam&, bool);
	long insert(const Protein&);
	
	//properties
	//int previousOccurance(const Protein&) const;
	
public:
	vector<string> colNames;
	
	//constructor
	Proteins(const FilterFileParams&);
	Proteins();
	
	bool readInLocDB(string);
	string locSearch(const DBProtein&) const;
	
	//properities
	bool writeOut(string, const FilterFileParams&) const;
	bool writeOutDB(string, const FilterFileParams&) const;
	
	//modifiers
	bool readIn(string, const FilterFileParams&);
	void addSubcelluarLoc();
	void calcMW(const MWDB&);
};

/* #################### calcMW.cpp #################### */

class Peptide{
	friend class MWDB;
private:
	string ID, sequence;
public:
	Peptide();
	//Peptide(string);
	
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

class MWDB{
public:
	HashTable <Peptide> peptideLibrary;
	BinTree <AminoAcid> aminoAcidsDB;
	
	//constructor
	MWDB();
	~MWDB();
	
	//modifers
	bool readIn(string, const FilterFileParams&);
	
	//properties
	double calcMW(string, int) const;
	string getSequence(string) const;
	double getMW(string, int) const;
	double getMW(char, int) const;
	
private:	
	//modofers
	bool readInAAs(string, string);
	bool readInPeptides(string);
	bool readInAADB(string);
	bool addStaticMod(const AminoAcid&);
};

/*************/
/* functions */
/*************/

/* dtafilter.cpp */
bool isColumnHeaderLine(const vector<string>&);
string parseSample(string, string, bool);
int parsePeptideSC(string);
string parseReplicate(string);
string getID(string);

/* utils.cpp */
bool dirExists (string);
bool fileExists (string);
string toString(int);
int toInt(string);
inline bool strContains(string, string);
inline bool strContains(char, string);
void split (const string, char, vector<string> &);
inline string trimTraling(const string&);
inline string trimLeading(const string&);
inline string trim(const string&);
bool isCommentLine(string);
bool isInteger(string);
inline void getLineTrim(ifstream&, string&);
string removeSubstr(string, string);
string toLower(string);
template<class T> long binSearch(const vector<T>&, const T&, long, long);
template<class T> typename vector<T>::iterator insertSorted(vector<T>& vec, const T& item);

