/*
 DTarray_AJM reads in a specified number of dtaselect-filter files and writes protein their molecular
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
string const DTA_PARAMS_NAME = "dtarray_ajm.params";
string const OF_NAME = "DTarray_AJM.txt";
int const FILTER_FILE_PARAMS_SIZE = 50;
bool const INCLUDE_FULL_DESCRIPTION = true;
//string const DEFAULT_COL_NAMES [] = {"IPI", "Description", "Mass (Da)"};
string const DEFAULT_COL_NAMES [] = {"Protein","IPI", "Description", "Mass (Da)"};
int const DEFAULT_COL_NAMES_LENGTH = 4;
string const COLUMN_HEADER_LINE_ELEMENTS[] = {"Unique", "FileName", "XCorr", "DeltCN", "Conf%", "M+H+",
    "CalcM+H+", "TotalIntensity", "SpR", "ZScore", "IonProportion", "Redundancy", "Sequence"};
int const COLUMN_HEADER_LINE_ELEMENTS_SIZE = 13;

//editable params for DB output format
bool const PARSE_SAMPLE_NAME = true;
string const DEFAULT_COL_NAMES_DB [] = {"Protein","Match dirrection","IPI", "Description", "Mass (Da)", "Long sample name",
	"Spectral counts", "Sample", "Replicate"};
int const DEFAULT_COL_NAMES_DB_LENGTH = 9;
string const SAMPLE_NAME_PREFIX = "Biotin-PG_Tryp_";
string const UNIQUE_PEPTIDE_HEADERS[] = {"SC", "Unique pep. SC"};

//function definitions
bool strContains(char, string);
void split (const string, char, vector<string> &);
bool isColumnHeaderLine(const vector<string>&);
bool dirExists (string);
string parseSample(string);
int parseUniquePeptides(string);
string parseReplicate(string);

//class definitions and functions
struct FilterFileParam{
    string path;
    string colname;
};

struct FilterFileParams{
    FilterFileParam file [FILTER_FILE_PARAMS_SIZE];
    int numFiles;
    
    //modifiers
    bool readDTParams(string, string);
};

bool FilterFileParams::readDTParams(string fname, string path)
{
    ifstream inF ((path + fname).c_str());
    
    numFiles = 0;
    
    if (!inF)
        return false;
    
    for (int i = 0; i < FILTER_FILE_PARAMS_SIZE; i++)
    {
        getline(inF, file[i].colname, '\t');
        getline(inF, file[i].path);
        if (file[i].colname != "")
            numFiles ++;
    }
    
    return true;
}

struct FilterFile{
    string colname, count;
    string coverage, peptides, uniquePeptides;
    
    //modifiers
    FilterFile (string, string, string);
};

FilterFile::FilterFile(string arg1, string arg2, string arg3)
{
    colname = arg1;
    count = arg2;
    uniquePeptides = arg3;
}

struct Protein{
    string fullDescription, matchDirrection, IPI, description, MW;
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
    
    //Extract uniprotID
    size_t firstBar = elems[0].find("|");
    size_t secBar = elems[0].find("|", firstBar+1);
    IPI = line.substr(firstBar+1, secBar-firstBar-1);
    
    //extract matchDirrection
    matchDirrection = line.substr(0, firstBar);
    
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
    bool writeOut(string, bool) const;
    bool writeOutDB(string, bool) const;
    
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
                        numUniquePeptides += parseUniquePeptides(line);
                } while(!strContains('%', line) && !inF.eof());
                proteins[uniquePeptidesIndex].col[colIndex].uniquePeptides = to_string(numUniquePeptides);
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
        if (proteins[i].IPI == newProtein.IPI && proteins[i].matchDirrection == newProtein.matchDirrection)
            return i;
    
    //if newProtein is not found return -1
    return -1;
}

//write out combined protein lists to
bool Proteins::writeOut(string ofname, bool includeUnique) const
{
    ofstream outF (ofname.c_str());
    
    if(!outF)
        return false;
    
    
    //print header lines
    if (includeUnique)
    {
        for (int i = 0; i < DEFAULT_COL_NAMES_LENGTH; i++)
            outF << '\t';
        for (int i = 0; i < colNames.size(); i++)
        {
            if (i == 0)
                outF << colNames[i];
            else outF << '\t' <<'\t' << colNames[i];
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
            ofColNames.push_back(colNames[i]);
        int colNamesLen = int(ofColNames.size());
        for (int i = 0; i < colNamesLen; i++)
        {
            outF << ofColNames[i] << '\t';
        }
        outF << endl;
    }
    
    //print proteins and spectral counts
    int proteinsLen = int(proteins.size());
    for (int i = 0; i < proteinsLen; i++)
    {
        if (INCLUDE_FULL_DESCRIPTION)
            outF << proteins[i].fullDescription << '\t';
        
        outF << proteins[i].IPI << '\t' <<
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

bool Proteins::writeOutDB(string ofname, bool includeUnique) const
{
    ofstream outF (ofname.c_str());
    
    if(!outF)
        return false;
    
    //popuate ofColNames with default col names and unique peptide headers if necissary
    //and print first line of report
    vector<string> ofColNames;
    for (int i = 0; i < DEFAULT_COL_NAMES_DB_LENGTH; i ++)
        ofColNames.push_back(DEFAULT_COL_NAMES_DB[i]);
    if(includeUnique)
        ofColNames.push_back(UNIQUE_PEPTIDE_HEADERS[1]);
    
    for (int i = 0; i < int(ofColNames.size()); i++)
    {
        outF << ofColNames[i] << '\t';
    }
    outF << endl;
    
    //print proteins and spectral counts
    int proteinsLen = int(proteins.size());
    for (int i = 0; i < colNames.size(); i++)
    {
        for (int j = 0; j < proteinsLen; j++)
        {
            if (INCLUDE_FULL_DESCRIPTION)
                outF << proteins[j].fullDescription << '\t';
            
            outF << proteins[j].matchDirrection << '\t' <<
            proteins[j].IPI << '\t' <<
            proteins[j].description << '\t' <<
            proteins[j].MW << '\t' <<
            proteins[j].col[i].colname << '\t' <<
            proteins[j].col[i].count;
            
            if (PARSE_SAMPLE_NAME)
            {
                outF << '\t' << parseSample(proteins[j].col[i].colname) << '\t' <<
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
    if (argc == 4)
    {
        string includeUniqueStr = argv[3];
        if (includeUniqueStr != "0" && includeUniqueStr != "1")
        {
            cout << includeUniqueStr << " is not a valid arguement for includeUnique. Exiting..." << endl;
            return 0;
        }
        includeUnique = stoi(argv[3]);
    }
    
    //read in names of files to combine from params file
    FilterFileParams filterFileParams;
    if (!filterFileParams.readDTParams(DTA_PARAMS_NAME, wd))
    {
        cout <<"Failed to read params file! Exiting..." << endl;
        return 0;
    }
	
	//combine files
    cout << endl;
    Proteins proteins;
    proteins.initialize(filterFileParams);
    for (int i = 0; i < filterFileParams.numFiles; i++)
    {
        if(!proteins.readIn(wd+filterFileParams.file[i].path, filterFileParams.file[i].colname,
                            includeUnique))
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
        if(!proteins.writeOut(wd + OF_NAME, includeUnique))
        {
            cout << "Could not write outFile! Exiting..." << endl;
            return 0;
        }
    }
    else if (outputFormat == "DB")
    {
        if(!proteins.writeOutDB(wd + OF_NAME, includeUnique))
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
} //end main

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
    assert(len <= COLUMN_HEADER_LINE_ELEMENTS_SIZE);
    
    for (int i = 0; i < len; i++)
        if(COLUMN_HEADER_LINE_ELEMENTS[i] != elems[i])
            return false;
    
    return true;
}

//returns true if folder at end of path exists and returns false if it does not
bool dirExists (string path)
{
    struct stat buffer;
    if (stat(path.c_str(), &buffer) == 0 && S_ISDIR(buffer.st_mode))
        return true;
    return false;
}

//optional fxn to parse long sample name
string parseSample(string str)
{
    string sample = str.substr(SAMPLE_NAME_PREFIX.length());
    sample = sample.substr(0, sample.find_last_of("_"));
    size_t posBegin = sample.find("_");
    size_t posEnd = sample.find_last_of("_");
    sample = sample.substr(0, posBegin) + " " + sample.substr(posEnd + 1);
    
    return sample;
}

string parseReplicate(string str)
{
    return str.substr(str.find_last_of("_")+1);
}

int parseUniquePeptides(string line)
{
    //split line by tabs
    vector<string> elems;
    split(line, '\t', elems);
    
    //return SC for peptide as int
    return stoi(elems[11]);
}


