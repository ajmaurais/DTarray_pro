/*
 DTarray_AJM reads in a specified number of dtaselect-filter files and writes protein, their molecular
 weights and spectral count data to OF_NAME in the working dirrectory.
 
 Written by Aaron Maurais
 25 May 2016
 */

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include <sys/stat.h>

using namespace std;

//editable paramaters
string const DTA_PARAMS_NAME = "dtarray_ajm.params"; //if value is changed DTarray.sh must also be updated
string const OF_NAME = "DTarray_AJM.txt";
bool const INCLUDE_FULL_DESCRIPTION = true;
//string const DEFAULT_COL_NAMES [] = {"ID", "Description", "Mass (Da)"};
string const DEFAULT_COL_NAMES [] = {"Protein","ID", "Description", "Mass (Da)"};
int const DEFAULT_COL_NAMES_LENGTH = 4;
string const COLUMN_HEADER_LINE_ELEMENTS[] = {"Unique", "FileName", "XCorr", "DeltCN", "Conf%", "M+H+",
	"CalcM+H+", "TotalIntensity", "SpR", "ZScore", "IonProportion", "Redundancy", "Sequence"};
int const COLUMN_HEADER_LINE_ELEMENTS_LENGTH = 13;

//editable params for DB output format
string const DEFAULT_COL_NAMES_DB [] = {"Protein","ID", "Description", "Mass (Da)", "Long sample name",
	"Spectral counts", "Sample", "Replicate"};
int const DEFAULT_COL_NAMES_DB_LENGTH = 8;
string const UNIQUE_PEPTIDE_HEADERS[] = {"SC", "Unique pep. SC"};

//function definitions
bool strContains(char, string);
void split (const string, char, vector<string> &);
bool isColumnHeaderLine(const vector<string>&);
bool dirExists (string);
string parseSample(string, string);
int parsePeptideSC(string);
string parseReplicate(string);
string toString(int);
int toInt(string);

//class definitions and member functions
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

Param::Param(string line)
{
	size_t posStart = line.find("=");
	
	param = line.substr(0, posStart);
	value = line.substr(posStart + 1);
}

struct FilterFileParams{
	vector<FilterFileParam> file;
	int numFiles;
	string sampleNamePrefix;
	
	//modifiers
	bool readDTParams(string, string);
};

bool FilterFileParams::readDTParams(string fname, string path)
{
	ifstream inF ((path + fname).c_str());
	
	int i = 0;
	numFiles = 0;
	string line;
	
	if (!inF)
		return false;
	
	while(!inF.eof())
		{
		getline(inF, line);
		if(strContains('=', line))
			{
			Param param (line);
			if(param.param == "sampleNamePrefix")
				sampleNamePrefix = param.value;
			else return false;
			}
		else {
			vector<string>elems;
			split(line, '\t', elems);
			if(elems.size() == 2)
			{
				FilterFileParam blank;
				file.push_back(blank);
				file[i].colname = elems[0];
				file[i].path = elems[1];
				numFiles++;
				i++;
			}
			else if (elems.size() != 0)
				return false;
		}
		}
	
	return true;
}

struct FilterFile{
	string colname, count;
	string coverage, peptides, uniquePeptides;
	
	//constructor
	FilterFile (string, string, string);
};

FilterFile::FilterFile(string arg1, string arg2, string arg3)
{
	colname = arg1;
	count = arg2;
	uniquePeptides = arg3;
}

struct Protein{
	string fullDescription, matchDirrection, ID, description, MW;
	vector <FilterFile> col;
	
	//modifier
	void initialize(const vector<string>&);
	bool getProteinData(string, int);
	void consolidate(const Protein&, int);
};

void Protein::initialize(const vector<string>& colNames)
{
	int len = int(colNames.size());
	
	for (int i = 0; i < len; i++)
		{
		FilterFile newFilterFile(colNames[i], "0", "0");
		col.push_back(newFilterFile);
		}
}

//parse proten header line and extract desired data
bool Protein::getProteinData(string line, int colIndex)
{
	//split line by tabs
	vector<string> elems;
	split(line, '\t', elems);
	
	if(isColumnHeaderLine(elems))
		return false;
	
	//keep fullDescription but seperate by spaces instead of tabs
	fullDescription = elems[0];
	int len = int(elems.size());
	for (int i = 1; i < len; i++)
		fullDescription += (" " + elems[i]);
	
	//extract matchDirrection
	size_t firstBar = elems[0].find("|");
	matchDirrection = line.substr(0, firstBar);
	
	//Extract uniprotID
	size_t secBar = elems[0].find("|", firstBar+1);
	ID = line.substr(firstBar+1, secBar-firstBar-1);
	if(matchDirrection == "Reverse_sp")
		ID = "reverse_" + ID;
	
	//extract MW
	MW = elems[5];
	
	//extract shortened protein description
	size_t endOfDescription = elems[8].find(" [");
	description = elems[8].substr(0, endOfDescription);
	
	//add spectrum count for *this protein to colname
	col[colIndex].count = elems[2];
	
	return true;
}

void Protein::consolidate(const Protein& toAdd, int colIndex)
{
	col[colIndex] = toAdd.col[colIndex];
}

struct Proteins{
	vector <Protein> proteins;
	vector<string> colNames;
	int colIndex;
	
	//properities
	int previousOccurance(const Protein&) const;
	bool writeOut(string, bool, bool, string) const;
	bool writeOutDB(string, bool, bool, string) const;
	
	//modifiers
	void initialize(const FilterFileParams&);
	bool readIn(string, string, bool);
	void addBlanks();
};

//initialize Proteins.colnames with files contained in params file
void Proteins::initialize(const FilterFileParams& files)
{
	colIndex = 0;
	
	for (int i = 0; i < files.numFiles; i++)
		colNames.push_back(files.file[i].colname);
}

//read in protein headder lines and parse with getProteinData
bool Proteins::readIn(string fname, string colname, bool countUniquePeptides)
{
	ifstream inF(fname.c_str());
	if(!inF)
		return false;
	
	int proteinsIndex = int(proteins.size());
	int previousOccuranceIndex;
	int numUniquePeptides = 0;
	int uniquePeptidesIndex = -1;
	bool inProtein = false;
	bool getNewLine = true;
	Protein blank;
	blank.initialize(colNames);
	string line;
	
	while(!inF.eof()){
		if(getNewLine)
			getline(inF, line);
		getNewLine = true;
		if(strContains('%', line))  //find protein header lines by percent symbol for percent coverage
									//if(strContains('|', line) && strContains('%', line))  //alternativly use both | and % symbols but may loose
									//some uncharacterized proteins
			{
			Protein newProtein;
			newProtein.initialize(colNames);
			if(newProtein.getProteinData(line, colIndex))
				{
				inProtein = true;
				previousOccuranceIndex = previousOccurance(newProtein);
				if (previousOccuranceIndex == -1)
					{
					proteins.push_back(blank);
					proteins[proteinsIndex] = newProtein;
					uniquePeptidesIndex = proteinsIndex;
					proteinsIndex++;
					}
				else {
					proteins[previousOccuranceIndex].consolidate(newProtein, colIndex);
					uniquePeptidesIndex = previousOccuranceIndex;
				}
				}
			if(countUniquePeptides && inProtein)
				{
				do{
					getline(inF, line);
					if(line[0] == '*')
						numUniquePeptides += parsePeptideSC(line);
				} while(!strContains('%', line) && !inF.eof());
				proteins[uniquePeptidesIndex].col[colIndex].uniquePeptides = toString(numUniquePeptides);
				numUniquePeptides = 0;
				getNewLine = false;
				}
			}
	}
	colIndex++;
	return true;
}

//itterate through all previous proteins and return index at which a previous occurance
//of newProtein occures
int Proteins::previousOccurance(const Protein& newProtein) const
{
	int len = int(proteins.size());
	
	for (int i = 0; i < len; i++)
		if (proteins[i].ID == newProtein.ID && proteins[i].matchDirrection == newProtein.matchDirrection)
			return i;
	
	//if newProtein is not found return -1
	return -1;
}

//write out combined protein lists to ofname
bool Proteins::writeOut(string ofname, bool includeUnique, bool parseSampleName, string samplePrefix) const
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	//print header lines
	if(parseSampleName)
		{
		string delim;
		if(includeUnique)
			delim = "\t\t";
		else delim = "\t";
		for (int i = 0; i < DEFAULT_COL_NAMES_LENGTH; i++)
			outF << '\t';
		for (int i = 0; i < colNames.size() ; i++)
			outF << colNames[i] << delim;
		outF << endl;
		}
	if (includeUnique)
		{
		for (int i = 0; i < DEFAULT_COL_NAMES_LENGTH; i++)
			outF << '\t';
		for (int i = 0; i < colNames.size(); i++)
			{
			if (i == 0)
				outF << parseSample(colNames[i], samplePrefix);
			else outF << '\t' <<'\t' << parseSample(colNames[i], samplePrefix);
			}
		outF << endl;
		for (int i = 0; i < DEFAULT_COL_NAMES_LENGTH; i++)
			outF << DEFAULT_COL_NAMES[i] <<'\t';
		for (int i = 0; i < colNames.size(); i++)
			outF << UNIQUE_PEPTIDE_HEADERS[0] << '\t' << UNIQUE_PEPTIDE_HEADERS[1] << '\t';
		outF << endl;
		}
	else
		{
		vector<string> ofColNames;
		for (int i = 0; i < DEFAULT_COL_NAMES_LENGTH; i ++)
			ofColNames.push_back(DEFAULT_COL_NAMES[i]);
		for (int i = 0; i < colNames.size(); i ++)
			ofColNames.push_back(parseSample(colNames[i], samplePrefix));
		int colNamesLen = int(ofColNames.size());
		for (int i = 0; i < colNamesLen; i++)
			outF << ofColNames[i] << '\t';
		outF << endl;
		}
	
	//print proteins and spectral counts
	int proteinsLen = int(proteins.size());
	for (int i = 0; i < proteinsLen; i++)
		{
		if (INCLUDE_FULL_DESCRIPTION)
			outF << proteins[i].fullDescription << '\t';
		
		outF << proteins[i].ID << '\t' <<
		proteins[i].description << '\t' <<
		proteins[i].MW << '\t';
		
		for (int j = 0; j < colIndex; j++)
			{
			outF << proteins[i].col[j].count << '\t';
			if (includeUnique)
				outF << proteins[i].col[j].uniquePeptides << '\t';
			}
		
		outF << endl;
		}
	
	
	return true;
}

bool Proteins::writeOutDB(string ofname, bool includeUnique, bool parseSampleName, string sampleName) const
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	//popuate ofColNames with default col names and unique peptide headers if necissary
	//and print first line of report
	vector<string> ofColNames;
	for (int i = 0; i < DEFAULT_COL_NAMES_DB_LENGTH - (!parseSampleName * 2); i ++)
		ofColNames.push_back(DEFAULT_COL_NAMES_DB[i]);
	if(includeUnique)
		ofColNames.push_back(UNIQUE_PEPTIDE_HEADERS[1]);
	
	for (int i = 0; i < int(ofColNames.size()); i++)
		outF << ofColNames[i] << '\t';
	outF << endl;
	
	//print proteins and spectral counts
	int proteinsLen = int(proteins.size());
	for (int i = 0; i < colNames.size(); i++)
		{
		for (int j = 0; j < proteinsLen; j++)
			{
			if (INCLUDE_FULL_DESCRIPTION)
				outF << proteins[j].fullDescription << '\t';
			
			outF << proteins[j].ID << '\t' <<
			proteins[j].description << '\t' <<
			proteins[j].MW << '\t' <<
			proteins[j].col[i].colname << '\t' <<
			proteins[j].col[i].count;
			
			if (parseSampleName)
				{
				outF << '\t' << parseSample(proteins[j].col[i].colname, sampleName) << '\t' <<
				parseReplicate(proteins[j].col[i].colname);
				}
			
			if (includeUnique)
				{
				outF << '\t' << proteins[j].col[i].uniquePeptides;
				}
			
			outF << endl;
			}
		}
	
	return true;
}
//end class definitions and member functions

//begin main
int main (int argc, char *argv[])
{
	//check paramaters
	string wd = string(argv[1]);
	assert(dirExists(wd));
	string outputFormat = "standard";
	if (argc >= 3)
		{
		outputFormat = argv[2];
		if(outputFormat != "standard" && outputFormat != "DB")
			{
			cout << outputFormat << " is not a valid output format. Exiting..." << endl;
			return 0;
			}
		}
	int includeUnique = 0;
	if (argc >= 4)
		{
		string includeUniqueStr = argv[3];
		if (includeUniqueStr != "0" && includeUniqueStr != "1")
			{
			cout << includeUniqueStr << " is not a valid arguement for includeUnique. Exiting..." << endl;
			return 0;
			}
		includeUnique = toInt(includeUniqueStr);
		}
	int parseSampleName = 0;
	if (argc >= 5)
		{
		string parseSampleNameStr = argv[4];
		if (parseSampleNameStr != "0" && parseSampleNameStr != "1")
			{
			cout << parseSampleNameStr << " is not a valid arguement for parseSampleName. Exiting..." << endl;
			return 0;
			}
		parseSampleName = toInt(parseSampleNameStr);
		}
	
	//read in names of files to combine from params file
	FilterFileParams filterFileParams;
	if (!filterFileParams.readDTParams(DTA_PARAMS_NAME, wd))
		{
		cout <<"Failed to read params file! Exiting..." << endl;
		return 0;
		}
	if(parseSampleName)
		cout << "Parsing colnames by prefix: " << filterFileParams.sampleNamePrefix << endl;
	
	//combine files
	cout << endl;
	Proteins proteins;
	proteins.initialize(filterFileParams);
	for (int i = 0; i < filterFileParams.numFiles; i++)
		{
		if(!proteins.readIn(wd+filterFileParams.file[i].path, filterFileParams.file[i].colname, includeUnique))
			{
			cout <<"Failed to read in " << filterFileParams.file[i].path <<"!" << endl <<
			"Exiting..." << endl;
			return 0;
			}
		cout << "Adding " << filterFileParams.file[i].colname << "..." << endl;
		}
	
	//write out combined data to OF_NAME
	if (outputFormat == "standard")
		{
		if(!proteins.writeOut(wd + OF_NAME, includeUnique, parseSampleName, filterFileParams.sampleNamePrefix))
			{
			cout << "Could not write outFile! Exiting..." << endl;
			return 0;
			}
		}
	else if (outputFormat == "DB")
		{
		if(!proteins.writeOutDB(wd + OF_NAME, includeUnique, parseSampleName, filterFileParams.sampleNamePrefix))
			{
			cout << "Could not write outFile! Exiting..." << endl;
			return 0;
			}
		cout << endl << "Results written in database format." << endl;
		}
	else
		{
		cout << endl << outputFormat << " is not a valid output format! Exiting..." << endl;
		return 0;
		}
	
	//summarize results for user
	cout << proteins.colNames.size() << " files combined." << endl;
	cout << "Results written to: " << OF_NAME << endl;
	
	return 0;
}
//end main

//search through string for char and return true if char is found
bool strContains(char findTxt, string whithinTxt)
{
	int len = int(whithinTxt.length());
	
	for(int i = 0; i < len; i++)
		if(whithinTxt[i] == findTxt)
			return true;
	
	return false;
}

//split str by delim and populate each split into elems
void split (const string str, char delim, vector<string> & elems)
{
	stringstream ss (str);
	string item;
	
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
}

//check if line containing % is a collumn header line instead of a protein header line
bool isColumnHeaderLine(const vector<string>& elems)
{
	int len = int(elems.size());
	assert(len <= COLUMN_HEADER_LINE_ELEMENTS_LENGTH);
	
	for (int i = 0; i < len; i++)
		if(COLUMN_HEADER_LINE_ELEMENTS[i] != elems[i])
			return false;
	
	return true;
}

//returns true if folder at end of path exists and false if it does not
bool dirExists (string path)
{
	struct stat buffer;
	if (stat(path.c_str(), &buffer) == 0 && S_ISDIR(buffer.st_mode))
		return true;
	return false;
}

//optional fxn to parse long sample name
string parseSample(string sampleName, string prefix)
{
	//return unparsed sampleName if prefix is empty string or is not found in sampleName
	if(sampleName.find(prefix) == string::npos || prefix.length() == 0)
		return sampleName;
	
	string sample = sampleName.substr(prefix.length());
	size_t posBegin = sample.find("_");
	size_t posEnd = sample.find_last_of("_");
	sample = sample.substr(0, posBegin) + " " + sample.substr(posBegin+1, sample.length() - posEnd);
	
	return sample;
}

//get replicate number from sample name
string parseReplicate(string str)
{
	return str.substr(str.find_last_of("_")+1);
}

//return number of spectral counts from a
int parsePeptideSC(string line)
{
	//split line by tabs
	vector<string> elems;
	split(line, '\t', elems);
	
	//return SC for peptide as int
	return toInt(elems[11]);
}

//converts int to string because to_string does not work with some c++ compilers
string toString(int num)
{
	string str;
	stringstream convert;
	
	convert << num;
	convert >> str;
	
	return str;
}

//converts string to int because atoi does not work with some c++ compilers
int toInt(string str)
{
	int num;
	stringstream convert;
	
	convert << str;
	convert >> num;
	
	return num;
}
