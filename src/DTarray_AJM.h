
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <functional>

//#define nullptr NULL //to compile on pleiades this line must be included

using namespace std;

/**********************/
/* global constants */
/*********************/

/* dtafilter.cpp */
string const OF_NAME = "DTarray_AJM.txt";
bool const INCLUDE_FULL_DESCRIPTION = true;
//string const DEFAULT_COL_NAMES [] = {"ID", "Description", "Mass (Da)"};
string const DEFAULT_COL_NAMES [] = {"Protein","ID", "Description", "Mass (Da)", "subcellular location"};
int const DEFAULT_COL_NAMES_LENGTH = 4;
string const COLUMN_HEADER_LINE_ELEMENTS[] = {"Unique", "FileName", "XCorr", "DeltCN", "Conf%", "M+H+",
	"CalcM+H+", "TotalIntensity", "SpR", "ZScore", "IonProportion", "Redundancy", "Sequence"};
int const COLUMN_HEADER_LINE_ELEMENTS_LENGTH = 13;
string const PARAM_ERROR_MESSAGE = " is an invalid arguement for: ";
int const MAX_PARAM_ITTERATIONS = 100;

//editable params for DB output format
string const DEFAULT_COL_NAMES_DB [] = {"Protein","ID", "Description", "Mass (Da)", "Long sample name",
	"Spectral counts", "Sample", "Replicate"};
int const DEFAULT_COL_NAMES_DB_LENGTH = 8;
string const DEFAULT_COL_NAMES_DB_LOC [] = {"Protein","ID", "Description", "Mass (Da)", "subcellular location", "Long sample name",
	"Spectral counts", "Sample", "Replicate"};
int const DEFAULT_COL_NAMES_DB_LOC_LENGTH = 9;
string const UNIQUE_PEPTIDE_HEADERS[] = {"SC", "Unique pep. SC"};
string const MWCALC_HEADERS [] = {"sequence", "avg mass", "monoisotopic mass"};
int const MWCALC_HEADERS_LENGTH = 3;

/* utils.cpp */
string const WHITESPACE = " \f\n\r\t\v";
string const COMMENT_SYMBOL = "#"; //if changed, paramsCommentSymbol must also be changed in DTarray_AJM.sh


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
class peptide;

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

/* #################### calcMW.cpp #################### */

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

class Peptide{
	friend class MWDB;
private:
	int ID;
	string sequence;
	
public:
	Peptide(string);
};

class mwNode{
	friend class MWDB;
private:
	AminoAcid aa;
	mwNode();
	
	mwNode *left;
	mwNode *right;
};

class MWDB{
public:
	//constructors
	MWDB();
	~MWDB();
	
	//modifers
	void insert(const AminoAcid&);
	bool readIn(string, const FilterFileParams&);
	
	//properties
	double getMW(string, int) const;
	double getMW(char, int) const;
	mwNode *search(const AminoAcid&) const;
	double calcMW(string, int) const;
	string getSequence(string) const;
	
private:
	vector<Peptide> peptides;
	mwNode *root;
	
	//modofers
	bool readInPeptides(string);
	bool readInAAs(string, string);
	bool readInAADB(string);
	bool addStaticMod(const AminoAcid&);
	void destroyTree(mwNode *leaf);
	void destroyTree();
	void insert(const AminoAcid&, mwNode *);
	
	//properties
	mwNode *search(const AminoAcid&, mwNode *) const;
	int search(string, int, int) const;
};

/*************/
/* functions */
/*************/

/* dtafilter.cpp */
bool isColumnHeaderLine(const vector<string>&);
string parseSample(string, string, string);
int parsePeptideSC(string);
string parseReplicate(string);
string parseIDPrefix(string, string);

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

/* calcMW.cpp */
