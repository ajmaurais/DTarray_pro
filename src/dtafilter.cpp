
Param::Param(string line)
{
	size_t posStart = line.find("=");
	
	param = line.substr(0, posStart);
	value = line.substr(posStart + 1);
}

//read in files to combine and output paramaters from params file
bool FilterFileParams::readDTParams(string fname, string path)
{
	ifstream inF ((path + fname).c_str());
	
	if (!inF)
		return false;
	
	string line;
	
	do{
		getLineTrim(inF, line);
		if(isCommentLine(line) || line.empty()) //skip line if is comment line
			continue;
		
		if(line == "<params>")
			do{
				getLineTrim(inF, line);
				if(strContains('=', line)) //find lines containing params by = symbol
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
					if(param.param == "calcMW")
					{
						assert(param.value == "0" || param.value == "1");
						calcMW = toInt(param.value);
						continue;
					}
					if(param.param == "peptideDBfname")
					{
						peptideDBfname = param.value;
						continue;
					}
					if(param.param == "aaDBfanme")
					{
						aaDBfanme = param.value;
						continue;
					}
					if(param.param == "staticModsFname")
					{
						staticModsFname = param.value;
						continue;
					}
					else return false;
				}
				else if(isCommentLine(line) || line.empty())
					continue;
				else if(line != "</params>")
					return false;
			} while(line != "</params>");
		
	} while(!inF.eof() && line != "</paramsFile>");
	return true;
}

bool FilterFileParams::readFlist(string fname, string path)
{
	ifstream inF ((path + fname).c_str());
	
	if (!inF)
		return false;
	
	int i = 0;
	numFiles = 0;
	string line;
	
	while(!inF.eof())
	{
		getline(inF, line);
		line = trim(line);
		if(isCommentLine(line) || line.empty()) //skip line if is comment line
			continue;
		else //else line contains data for filter file
		{
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

FilterFile::FilterFile(string arg1, string arg2, string arg3)
{
	colname = arg1;
	count = arg2;
	uniquePeptides = arg3;
}

//fill col element with colNames found in params file
void Protein::initialize(const vector<string>& colNames)
{
	col.clear();
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
	
	//tell readIn function to skip line if it is header line
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

void Protein::calcMW(const MWDB& mwDB)
{
	string sequence = mwDB.getSequence(parseIDPrefix(ID, "KL_"));
	
	avgMass = mwDB.calcMW(sequence, 0);
	monoMass = mwDB.calcMW(sequence, 1);
	calcSequence = sequence;
}

//convert Protein to DBprotein to allow inteface with BinTree of DBproteins
DBProtein Protein::toDBprotein() const
{
	DBProtein dbprotein;
	dbprotein.ID = ID;
	
	return dbprotein;
}

void Protein::consolidate(const Protein& toAdd, int colIndex)
{
	col[colIndex].colname = toAdd.col[colIndex].colname;
	col[colIndex].count = toAdd.col[colIndex].count;
}

Proteins::Proteins()
{
	colIndex = 0;
}

//initialize Proteins.colnames with files contained in params file
Proteins::Proteins(const FilterFileParams& files)
{
	colIndex = 0;
	
	for (int i = 0; i < files.numFiles; i++)
		colNames.push_back(files.file[i].colname);
}

//loop through protein headder and peptide lines in DTA filter file and add data to Proteins
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
		{
			Protein newProtein;
			newProtein.initialize(colNames);
			if(newProtein.getProteinData(line, colIndex)) //if line is not column header line populate protein to Proteins
			{
				inProtein = true; //used to determine if whether it is valid to loop through peptide lines below protein
								  //header line to extract unique peptide spectral counts
				previousOccuranceIndex = previousOccurance(newProtein); //find index at which newProtein occures in Proteins
				if (previousOccuranceIndex == -1) //if protein is not found, add to proteins dataset
				{
					proteins.push_back(blank);
					proteins[proteinsIndex] = newProtein;
					uniquePeptidesIndex = proteinsIndex;
					proteinsIndex++;
				}
				else //if protein is found, add spectral counts to appropiate filterFile
				{
					proteins[previousOccuranceIndex].consolidate(newProtein, colIndex);
					uniquePeptidesIndex = previousOccuranceIndex;
				}
			}
			//extract spectral counts for unique peptides
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
				inProtein = false;
			}
		}
	}
	colIndex++;
	return true;
}

//public Proteins::readIn function which adds all files in filterFileParams to Proteins
//and summarizes progress for user.
bool Proteins::readIn(string wd, const FilterFileParams& filterFileParams)
{
	for (int i = 0; i < filterFileParams.numFiles; i++)
	{
		if(!readIn(wd, filterFileParams.file[i], filterFileParams.includeUnique))
		{
			cout <<"Failed to read in " << filterFileParams.file[i].path <<"!" << endl << "Exiting..." << endl;
			return false;
		}
		cout << "Adding " << filterFileParams.file[i].colname << "..." << endl;
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

//write out combined protein lists to ofname in wide format
bool Proteins::writeOut(string ofname, const FilterFileParams& filterFileParams) const
{
	ofstream outF (ofname.c_str());
	bool parseSampleName = filterFileParams.sampleNamePrefix != "";
	int colNamesLength = DEFAULT_COL_NAMES_LENGTH;
	if(filterFileParams.getSubCelluarLoc)
		colNamesLength++;
	
	if(!outF)
		return false;
	
	//print column headers
	if(filterFileParams.calcMW)
	{
		
	}
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
			else outF << '\t' << '\t' << parseSample(colNames[i], filterFileParams.sampleNamePrefix, "standard");
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
		
		if(filterFileParams.calcMW)
			outF << proteins[i].calcSequence << '\t' <<
			proteins[i].avgMass << '\t' <<
			proteins[i].monoMass << '\t';
			
		
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

//write out combined protein lists to ofname in long format
bool Proteins::writeOutDB(string ofname, const FilterFileParams& filterFileParams) const
{
	ofstream outF (ofname.c_str());
	bool parseSampleName = (filterFileParams.sampleNamePrefix != "");
	
	if(!outF)
		return false;
	
	//print column headers
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

//seearches binary tree containing subcelluar localization data and populates data
//to loc element for each Protein in Proteins
void Proteins::addSubcelluarLoc(const BinTree& locDBinTree)
{
	int len = int(proteins.size());
	for (int i = 0; i < len; i++)
		proteins[i].loc = locDBinTree.locSearch(proteins[i].toDBprotein());
}

void Proteins::calcMW(const MWDB& mwDB)
{
	int len = int(proteins.size());
	for (int i = 0; i < len; i++)
		proteins[i].calcMW(mwDB);
}

//check if line containing % is a collumn header line instead of a protein header line
bool isColumnHeaderLine(const vector<string>& elems)
{
	//int len = int(elems.size());
	//assert(len <= COLUMN_HEADER_LINE_ELEMENTS_LENGTH);
	
	for (int i = 0; i < COLUMN_HEADER_LINE_ELEMENTS_LENGTH; i++)
		if(COLUMN_HEADER_LINE_ELEMENTS[i] != elems[i])
			return false;
	
	return true;
}

//optional fxn to parse long sample name
string parseSample(string sampleName, string prefix, string outputFormat)
{
	//return unparsed sampleName if prefix is empty string or is not found in sampleName
	//if(sampleName.find(prefix) == string::npos || prefix.length() == 0)
	if(strContains(prefix, sampleName) || prefix.length() == 0)
		return sampleName;
	
	//old version which adds space instead of _ between sample name
	//size_t posBegin = sample.find("_");
	//sample = sample.substr(0, posBegin) + " " + sample.substr(posBegin+1, sample.length() - posEnd);
	
	string sample = sampleName.substr(prefix.length()); //remove prefix from begining of sampleName
	if(outputFormat == "standard")
		return sample;
	if(outputFormat == "DB") //remove replicate number from sampleName if using DB outformat
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
	assert(elems.size() > 11);
	
	//return SC for peptide as int
	return toInt(elems[11]);
}

string parseIDPrefix(string str, string prefix)
{
	size_t begin = str.find(prefix);
	string ret = str.substr(begin + prefix.length());
	return ret;
}
