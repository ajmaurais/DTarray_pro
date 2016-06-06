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

//string const PARAMS_LOCATION = "/Users/Aaron/scripts/DTarray_AJM/DTarray_AJM_WD.params";
string const DTA_PARAMS_NAME = "dtarray_ajm.params";
string const OF_NAME = "DTarray_AJM.txt";
int const FILTER_FILE_PARAMS_SIZE = 50;
//string const DEFAULT_COL_NAMES [] = {"IPI", "Description", "Mass (Da)"};
string const DEFAULT_COL_NAMES [] = {"Protein","IPI", "Description", "Mass (Da)"};
int const DEFAULT_COL_NAMES_LENGTH = 4;
string const DEFAULT_COL_NAMES_DB [] = {"Protein","Match dirrection","IPI", "Description", "Mass (Da)", "Long sample name",
    "Spectral counts", "Sample", "replicate"};
int const DEFAULT_COL_NAMES_DB_LENGTH = 9;
string const SAMPLE_NAME_PREFIX = "Biotin-PG_Tryp_";
bool const PARSE_SAMPLE_NAME = true;
bool const INCLUDE_FULL_DESCRIPTION = true;
string const COLUMN_HEADER_LINE_ELEMENTS[] = {"Unique", "FileName", "XCorr", "DeltCN", "Conf%", "M+H+",
    "CalcM+H+", "TotalIntensity", "SpR", "ZScore", "IonProportion", "Redundancy", "Sequence"};
int const COLUMN_HEADER_LINE_ELEMENTS_SIZE = 13;

bool readParams(string&);
bool strContains(char, string);
void split (const string, char, vector<string> &);
bool isColumnHeaderLine(const vector<string>&);
bool dirExists (string);
string parseSample(string);
string parseReplicate(string);

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
    string colname;
    string count;
    
    //modifiers
    FilterFile (string, string);
};

FilterFile::FilterFile(string arg1, string arg2)
{
    colname = arg1;
    count = arg2;
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
        FilterFile newFilterFile(colNames[i], "0");
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

void Protein::consolidate(const Protein& destroy, int colIndex)
{
    col[colIndex] = destroy.col[colIndex];
}

struct Proteins{
    vector <Protein> proteins;
    vector<string> colNames;
    int colIndex;
    
    //properities
    int previousOccurance(const Protein&) const;
    bool writeOut(string) const;
    bool writeOutDB(string) const;
    
    //modifiers
    void initialize(const FilterFileParams&);
    bool readIn(string, string);
    void addBlanks();
};

void Proteins::initialize(const FilterFileParams& files)
{
    colIndex = 0;
    
    for (int i = 0; i < files.numFiles; i++)
        colNames.push_back(files.file[i].colname);
}

//read in protein headder lines and parse with getProteinData
bool Proteins::readIn(string fname, string colname)
{
    ifstream inF(fname.c_str());
    int proteinsIndex = int(proteins.size());
    int previousOccuranceIndex;
    Protein blank;
    blank.initialize(colNames);
    
    if(!inF)
        return false;
    
    string line;
    while(!inF.eof()){
        getline(inF, line);
        if(strContains('%', line))  //find protein header lines by percent symbol for percent coverage
            //if(strContains('|', line) && strContains('%', line))  //alternativly use both | and % symbols but may loose
            //some uncharacterized proteins
        {
            Protein newProtein;
            newProtein.initialize(colNames);
            if(newProtein.getProteinData(line, colIndex))
            {
                previousOccuranceIndex = previousOccurance(newProtein);
                if (previousOccuranceIndex == -1)
                {
                    proteins.push_back(blank);
                    proteins[proteinsIndex] = newProtein;
                    proteinsIndex++;
                }
                else {
                    proteins[previousOccuranceIndex].consolidate(newProtein, colIndex);
                }
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
bool Proteins::writeOut(string ofname) const
{
    ofstream outF (ofname.c_str());
    
    if(!outF)
        return false;
    
    //populate ofColNames with default collumn headings and variable headings
    //and print to first row of ofname
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
            outF << proteins[i].col[j].count << '\t';
        
        outF << endl;
    }
    
    
    return true;
}

bool Proteins::writeOutDB(string ofname) const
{
    ofstream outF (ofname.c_str());
    
    if(!outF)
        return false;
    
    //Print to first row of ofname
    for (int i = 0; i < DEFAULT_COL_NAMES_DB_LENGTH; i++)
    {
        outF << DEFAULT_COL_NAMES_DB[i] << '\t';
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
            
            outF << endl;
        }
    }
    
    return true;
    
}

int main (int argc, char *argv[])
{
    //read in working dirrectory from params file
    /*string wd;
     if (!readParams(wd))
     {
     cout <<"Failed to read params file!_1" << endl;
     return 0;
     }*/
    
    //cout << wd << endl;
    
    //read in paramaters
    //string wd = "/Users/Aaron/Google_Drive/School_Work/Mass_spec_data/Human_RA_and_Healthy_SF/";
    string wd = string(argv[1]);
    assert(dirExists(wd));
    string outputFormat = "standard";
    if (argc == 3)
        outputFormat = argv[2];
    
    //read in names of files to combine from params file
    FilterFileParams filterFileParams;
    if (!filterFileParams.readDTParams(DTA_PARAMS_NAME, wd))
    {
        cout <<"Failed to read params file!" << endl <<
        "Check that DTASelect-filter files exist in subdirectories." << endl;
        return 0;
    }
    
    //combine files
    Proteins proteins;
    proteins.initialize(filterFileParams);
    for (int i = 0; i < filterFileParams.numFiles; i++)
        if(!proteins.readIn(wd+filterFileParams.file[i].path+"DTASelect-filter.txt", filterFileParams.file[i].colname))
        {
            cout <<"Failed to read in " << filterFileParams.file[i].path <<"!" << endl;
            return 0;
        }
    
    //write out combined data to OF_NAME
    if (outputFormat == "standard")
    {
        if(!proteins.writeOut(wd + OF_NAME))
        {
            cout << "Could not write outFile!" << endl;
            return 0;
        }
    }
    else if (outputFormat == "DB")
    {
        if(!proteins.writeOutDB(wd + OF_NAME))
        {
            cout << "Could not write outFile!" << endl;
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
    cout << endl;
    for (int i = 0; i < proteins.colNames.size(); i++)
        cout << "Adding " << proteins.colNames[i] << "..." << endl;
    cout << proteins.colNames.size() << " files combined." << endl;
    cout << "Results written to: " << OF_NAME << endl;
    
    //cout << "Sucess!" << endl;
    return 0;
}

//read working dirrectory from params file in .scripts
/*bool readParams(string& path)
 {
 ifstream inF (PARAMS_LOCATION.c_str());
 
 if (!inF)
 return false;
 
 inF >> path;
 return true;
 }*/

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


