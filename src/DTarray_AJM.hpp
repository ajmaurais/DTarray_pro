
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

//#define nullptr NULL //to compile on pleiades this line must be included

using namespace std;

/******************************/
/* globally scoped constants */
/*****************************/

/* dtafilter.cpp */
string const OF_NAME = "DTarray_AJM.txt";
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

/* utils.cpp */
string const WHITESPACE = " \f\n\r\t\v";
string const COMMENT_SYMBOL = "#"; //if changed, paramsCommentSymbol must also be changed in DTarray_AJM.sh

/* hashTable.cpp */
int const HASH_TABLE_SIZE = 20000;
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
class BinTree;
class MWDB;
class nwNode;
class Peptide;
class LinkedList;
class HashTable;
class AATree;


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
	string aaDBfanme, peptideDBfname, staticModsFname;
	
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
private:
	friend class Proteins;
	vector <FilterFile> col;
	
	//modifier
	void initialize(const vector<string>&);
	bool getProteinData(string, int);
	void consolidate(const Protein&, int);
	DBProtein toDBprotein() const;
	void calcMW(const MWDB&);
	
	string fullDescription, matchDirrection, ID, description, MW, loc;
	double avgMass, monoMass;
	string calcSequence;
};

//stores data for all proteins found in DTA filter files
class Proteins{
	vector <Protein> proteins;
	int colIndex;
	
	//modifers
	bool readIn(string, const FilterFileParam&, bool);
	
	//properties
	int previousOccurance(const Protein&) const;
	
public:
	vector<string> colNames;
	
	//constructor
	Proteins(const FilterFileParams&);
	Proteins();
	
	//properities
	bool writeOut(string, const FilterFileParams&) const;
	bool writeOutDB(string, const FilterFileParams&) const;
	
	//modifiers
	bool readIn(string, const FilterFileParams&);
	void addSubcelluarLoc(const BinTree&);
	void calcMW(const MWDB&);
};

/* #################### subCelluarLoc.cpp #################### */

class DBProtein{
	friend class BinTree;
	friend class Protein;
	string ID, gene, description, loc;
	
	//modifers
	void operator = (const DBProtein&);
	void clear();
	
public:
	//constructor
	DBProtein(string);
	DBProtein();
	
	//properties
	bool operator < (const DBProtein&) const;
	bool operator > (const DBProtein&) const;
	bool operator == (const DBProtein&) const;
};

struct Node{
	DBProtein protein;
	Node();
	
	Node *left;
	Node *right;
};

class BinTree{
public:
	BinTree();
	~BinTree();
	
	void insert(const DBProtein&);
	bool readInProteins(string);
	string locSearch(const DBProtein&) const;
	
private:
	void destroyTree();
	void insert(const DBProtein&, Node *);
	Node *search(const DBProtein&, Node *) const;
	Node *search(const DBProtein&) const;
	void destroyTree(Node *leaf);
	
	Node *root;
};

/* #################### hashTable.cpp #################### */

class Peptide{
public:
	string ID, sequence;
	Peptide();
	Peptide(string);
	
	//modifers
	void operator = (const Peptide&);
	
	//properties
	bool operator == (string) const;
	string getID() const;
};

struct Item
{
	//struct or vars to sture in list
	Peptide val;
	Item* next;
	
	Item();
};

class LinkedList {
private:
	Item * head;
	int length;
	
	void destroyList(Item*);
	void insert(const Peptide&, Item*);
	Item * getItem(string, Item*) const;
	
public:
	//constructor
	LinkedList();
	~LinkedList();
	
	//modifer
	void insert(const Peptide&);
	//bool removeItem(string);
	void destroyList();
	
	//properties
	int getLength();
	Item * getItem(string) const;
};

class HashTable{
private:
	int size;
	LinkedList* array;
	hash<string> str_hash;
	
	int hash(string) const;
	void destroyTable(LinkedList*);
	void insertItem(Item*);
public:
	//coonstructor
	HashTable(int s = HASH_TABLE_SIZE);
	~HashTable();
	
	//modifers
	//bool remove(string);
	void destroyTable();
	void insert(const Peptide&);
	
	//properties
	int getLength();
	void printTable();
	Item * getItem(string) const;
	int getNumItems() const;
	void printHistogram() const;
};

/* #################### calcMW.cpp #################### */

class AminoAcid{
	friend class AATree;
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

class mwNode{
	friend class AATree;
private:
	AminoAcid aa;
	mwNode();
	
	mwNode *left;
	mwNode *right;
};

class AATree{
public:
	AATree();
	~AATree();
	
	bool readInAAs(string, string);
	void insert(const AminoAcid&);
	void destroyTree();
	mwNode *search(const AminoAcid&) const;
	
	double getMW(string, int) const;
	double getMW(char, int) const;
	
private:
	mwNode* root;
	
	bool readInAADB(string);
	bool addStaticMod(const AminoAcid&);
	void destroyTree(mwNode *leaf);
	void insert(const AminoAcid&, mwNode *);
	mwNode *search(const AminoAcid&, mwNode *) const;
};

class MWDB{
public:
	//constructors
	MWDB();
	~MWDB();
	
	//modifers
	bool readIn(string, const FilterFileParams&);
	
	//properties
	double calcMW(string, int) const;
	string getSequence(string) const;
	
	HashTable peptideLibrary;
	AATree aminoAcidsDB;
	
private:	
	//modofers
	bool readInPeptides(string);
	
	//properties
	int search(string, int, int) const;
};


/*************/
/* functions */
/*************/

/* dtafilter.cpp */
bool isColumnHeaderLine(const vector<string>&);
string parseSample(string, string, int);
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
inline int strComp(string, string);
string removeSubstr(string, string);
string toLower(string);
