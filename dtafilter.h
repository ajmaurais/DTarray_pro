
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include "utils.h"
#include "subCelluarLoc.h"

using namespace std;

//editable paramaters
string const OF_NAME = "DTarray_AJM.txt";
bool const INCLUDE_FULL_DESCRIPTION = true;
//string const DEFAULT_COL_NAMES [] = {"ID", "Description", "Mass (Da)"};
string const DEFAULT_COL_NAMES [] = {"Protein","ID", "Description", "Mass (Da)", "subcellular location"};
int const DEFAULT_COL_NAMES_LENGTH = 4;
string const COLUMN_HEADER_LINE_ELEMENTS[] = {"Unique", "FileName", "XCorr", "DeltCN", "Conf%", "M+H+",
	"CalcM+H+", "TotalIntensity", "SpR", "ZScore", "IonProportion", "Redundancy", "Sequence"};
int const COLUMN_HEADER_LINE_ELEMENTS_LENGTH = 13;
string const PARAM_ERROR_MESSAGE = " is an invalid arguement for: ";

//editable params for DB output format
string const DEFAULT_COL_NAMES_DB [] = {"Protein","ID", "Description", "Mass (Da)", "Long sample name",
	"Spectral counts", "Sample", "Replicate"};
int const DEFAULT_COL_NAMES_DB_LENGTH = 8;
string const DEFAULT_COL_NAMES_DB_LOC [] = {"Protein","ID", "Description", "Mass (Da)", "subcellular location", "Long sample name",
	"Spectral counts", "Sample", "Replicate"};
int const DEFAULT_COL_NAMES_DB_LOC_LENGTH = 9;
string const UNIQUE_PEPTIDE_HEADERS[] = {"SC", "Unique pep. SC"};

//function definitions
bool isColumnHeaderLine(const vector<string>&);
string parseSample(string, string, string);
int parsePeptideSC(string);
string parseReplicate(string);

//class definitions and member functions
class Protein;
class Proteins;

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
	
	string outputFormat;
	string sampleNamePrefix;
	bool includeUnique;
	bool getSubCelluarLoc;
	string locDBfname;
	
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
		line = trim(line);
		if(isCommentLine(line) || line == "")
			continue;
		if(strContains('=', line))
		{
			Param param (line);
			if(param.param == "sampleNamePrefix")
			{
				sampleNamePrefix = param.value;
				continue;
			}
			if(param.param == "outputFormat")
			{
				if(param.param != "standard" && param.param == "DB") {
					cout << param.param << PARAM_ERROR_MESSAGE << "outputFormat" << endl;
					return false;
				}
				outputFormat = param.value;
				continue;
			}
			if(param.param == "locDBfname")
			{
				locDBfname = param.value;
				continue;
			}
			if(param.param == "includeUnique")
			{
				assert(param.value == "0" || param.value == "1");
				includeUnique = toInt(param.value);
				continue;
			}
			if(param.param == "getSubCelluarLoc")
			{
				assert(param.value == "0" || param.value == "1");
				getSubCelluarLoc = toInt(param.value);
				continue;
			}
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

class Protein{
	friend class Proteins;
	vector <FilterFile> col;
	
	//modifier
	void initialize(const vector<string>&);
	bool getProteinData(string, int);
	void consolidate(const Protein&, int);
	DBProtein toDBprotein() const;
	
public:
	string fullDescription, matchDirrection, ID, description, MW, loc;
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

DBProtein Protein::toDBprotein() const
{
	DBProtein dbprotein;
	dbprotein.ID = ID;
	
	return dbprotein;
}

void Protein::consolidate(const Protein& toAdd, int colIndex)
{
	col[colIndex] = toAdd.col[colIndex];
}

class Proteins{
	vector <Protein> proteins;
	int colIndex;
	
	//modifers
	bool readIn(string, const FilterFileParam&, bool);
	
	//properties
	int previousOccurance(const Protein&) const;
	
public:
	vector<string> colNames;
	
	//properities
	bool writeOut(string, const FilterFileParams&) const;
	bool writeOutDB(string, const FilterFileParams&) const;
	
	//modifiers
	void initialize(const FilterFileParams&);
	bool readIn(string, const FilterFileParams&);
	void addBlanks();
	void addSubcelluarLoc(const Btree&);
};

//initialize Proteins.colnames with files contained in params file
void Proteins::initialize(const FilterFileParams& files)
{
	colIndex = 0;
	
	for (int i = 0; i < files.numFiles; i++)
		colNames.push_back(files.file[i].colname);
}

//read in protein headder lines and parse with getProteinData
//bool Proteins::readIn(string fname, string colname, bool countUniquePeptides)
bool Proteins::readIn(string wd, const FilterFileParam& filterFile, bool countUniquePeptides)
{
	string fname = wd + filterFile.path;
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

bool Proteins::readIn(string wd, const FilterFileParams& filterFile)
{
	for (int i = 0; i < filterFile.numFiles; i++)
	{
		if(!readIn(wd, filterFile.file[i], filterFile.includeUnique))
		{
			cout <<"Failed to read in " << filterFile.file[i].path <<"!" << endl <<
			"Exiting..." << endl;
			return false;
		}
		cout << "Adding " << filterFile.file[i].colname << "..." << endl;
	}
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
bool Proteins::writeOut(string ofname, const FilterFileParams& filterFileParams) const
{
	ofstream outF (ofname.c_str());
	bool parseSampleName = filterFileParams.sampleNamePrefix != "";
	int colNamesLength = DEFAULT_COL_NAMES_LENGTH;
	if(filterFileParams.getSubCelluarLoc)
		colNamesLength++;
	
	if(!outF)
		return false;
	
	//print header lines
	if(parseSampleName)
	{
		string delim;
		if(filterFileParams.includeUnique)
			delim = "\t\t";
		else delim = "\t";
		for (int i = 0; i < colNamesLength; i++)
			outF << '\t';
		for (int i = 0; i < colNames.size() ; i++)
			outF << colNames[i] << delim;
		outF << endl;
	}
	if (filterFileParams.includeUnique)
	{
		for (int i = 0; i < colNamesLength; i++)
			outF << '\t';
		for (int i = 0; i < colNames.size(); i++)
		{
			if (i == 0)
				outF << parseSample(colNames[i], filterFileParams.sampleNamePrefix, "standard");
			else outF << '\t' <<'\t' << parseSample(colNames[i], filterFileParams.sampleNamePrefix, "standard");
		}
		outF << endl;
		for (int i = 0; i < colNamesLength; i++)
			outF << DEFAULT_COL_NAMES[i] <<'\t';
		for (int i = 0; i < colNames.size(); i++)
			outF << UNIQUE_PEPTIDE_HEADERS[0] << '\t' << UNIQUE_PEPTIDE_HEADERS[1] << '\t';
		outF << endl;
	}
	else
	{
		vector<string> ofColNames;
		for (int i = 0; i < colNamesLength; i ++)
			ofColNames.push_back(DEFAULT_COL_NAMES[i]);
		for (int i = 0; i < colNames.size(); i ++)
			ofColNames.push_back(parseSample(colNames[i], filterFileParams.sampleNamePrefix, "standard"));
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
		
		if(filterFileParams.getSubCelluarLoc)
			outF << proteins[i].loc << '\t';
		
		for (int j = 0; j < colIndex; j++)
		{
			outF << proteins[i].col[j].count << '\t';
			if (filterFileParams.includeUnique)
				outF << proteins[i].col[j].uniquePeptides << '\t';
		}
		
		outF << endl;
	}
	
	
	return true;
}

bool Proteins::writeOutDB(string ofname, const FilterFileParams& filterFileParams) const
{
	ofstream outF (ofname.c_str());
	bool parseSampleName = (filterFileParams.sampleNamePrefix != "");
	
	if(!outF)
		return false;
	
	vector<string> headers;
	int len;
	if(filterFileParams.getSubCelluarLoc)
	{
		for (int i = 0; i < DEFAULT_COL_NAMES_DB_LOC_LENGTH ; i++)
			headers.push_back(DEFAULT_COL_NAMES_DB_LOC[i]);
		len = DEFAULT_COL_NAMES_DB_LOC_LENGTH;
	}
	else {
		for (int i = 0; i < DEFAULT_COL_NAMES_DB_LENGTH ; i++)
			headers.push_back(DEFAULT_COL_NAMES_DB[i]);
		len = DEFAULT_COL_NAMES_DB_LENGTH;
	}
	
	//popuate ofColNames with default col names and unique peptide headers if necissary
	//and print first line of report
	vector<string> ofColNames;
	for (int i = 0; i < len - (!parseSampleName * 2); i ++)
		ofColNames.push_back(headers[i]);
	if(filterFileParams.includeUnique)
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
			proteins[j].MW << '\t';
			
			if(filterFileParams.getSubCelluarLoc)
				outF << proteins[j].loc << '\t';
			
			outF << proteins[j].col[i].colname << '\t' <<
			proteins[j].col[i].count;
			
			if (parseSampleName)
			{
				outF << '\t' << parseSample(proteins[j].col[i].colname, filterFileParams.sampleNamePrefix, "DB") << '\t' <<
				parseReplicate(proteins[j].col[i].colname);
			}
			
			if (filterFileParams.includeUnique)
				outF << '\t' << proteins[j].col[i].uniquePeptides;
			
			outF << endl;
		}
	}
	
	return true;
}

void Proteins::addSubcelluarLoc(const Btree& locDBtree)
{
	int len = int(proteins.size());
	for (int i = 0; i < len; i++)
	{
		DBProtein temp = proteins[i].toDBprotein();
		string loc = locDBtree.locSearch(temp);
		proteins[i].loc = loc;
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


//optional fxn to parse long sample name
string parseSample(string sampleName, string prefix, string outputFormat)
{
	//return unparsed sampleName if prefix is empty string or is not found in sampleName
	if(sampleName.find(prefix) == string::npos || prefix.length() == 0)
		return sampleName;
	
	//old version which adds space instead of _ between sample name
	//size_t posBegin = sample.find("_");
	//sample = sample.substr(0, posBegin) + " " + sample.substr(posBegin+1, sample.length() - posEnd);
	
	string sample = sampleName.substr(prefix.length());
	if(outputFormat == "standard")
		return sample;
	if(outputFormat == "DB")
		return sample.substr(0, sample.find_last_of("_"));
	
	return sampleName;
}

//get replicate number from sample name
string parseReplicate(string str)
{
	return str.substr(str.find_last_of("_")+1);
}

//return number of spectral counts for a peptide from a line in DTA filter file
int parsePeptideSC(string line)
{
	//split line by tabs
	vector<string> elems;
	split(line, '\t', elems);
	
	//return SC for peptide as int
	return toInt(elems[11]);
}
