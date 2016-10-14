
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
		util::getLineTrim(inF, line);
		if(util::isCommentLine(line) || line.empty()) //skip line if is comment line
			continue;
		
		if(line == "<params>")
			do{
				util::getLineTrim(inF, line);
				if(util::strContains('=', line)) //find lines containing params by = symbol
				{
					Param param (line);
					if(param.param == "sampleNamePrefix")
					{
						sampleNamePrefix = param.value;
						continue;
					}
					if(param.param == "outputFormat")
					{
						if(param.param != "standard" && param.param == "db") {
							cout << param.param << PARAM_ERROR_MESSAGE << "outputFormat" << endl;
							return false;
						}
						outputFormat = util::toLower(param.value);
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
						includeUnique = util::toInt(param.value);
						continue;
					}
					if(param.param == "getSubCelluarLoc")
					{
						assert(param.value == "0" || param.value == "1");
						getSubCelluarLoc = util::toInt(param.value);
						continue;
					}
					if(param.param == "calcMW")
					{
						assert(param.value == "0" || param.value == "1");
						calcMW = util::toInt(param.value);
						continue;
					}
					if(param.param == "mwDBFname")
					{
						mwDBFname = param.value;
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
					if(param.param == "ofname")
					{
						ofname = param.value;
						continue;
					}
					if(param.param == "seqDBfname")
					{
						seqDBfname = param.value;
						continue;
					}
					if(param.param == "includeSeq")
					{
						assert(param.value == "0" || param.value == "1");
						includeSeq = util::toInt(param.value);
						continue;
					}
					else return false;
				}
				else if(util::isCommentLine(line) || line.empty())
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
		line = util::trim(line);
		if(util::isCommentLine(line) || line.empty()) //skip line if is comment line
			continue;
		else //else line contains data for filter file
		{
			vector<string>elems;
			util::split(line, '\t', elems);
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

bool Protein::operator == (const Protein& comp) const
{
	return comp.ID == ID;
}

bool Protein::operator > (const Protein& comp) const
{
	return comp.ID > ID;
}

bool Protein::operator < (const Protein& comp) const
{
	return comp.ID < ID;
}


//parse proten header line and extract desired data
bool Protein::getProteinData(string line, int colIndex)
{
	//split line by tabs
	vector<string> elems;
	util::split(line, '\t', elems);
	
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
	ID = getID(elems[0]);
	if(matchDirrection == "Reverse_sp")
		ID = "reverse_" + ID;
	
	//extract MW
	MW = elems[5];
	
	//extract shortened protein name and description
	size_t endOfDescription = elems[8].find(" [");
	string descriptionTemp = elems[8].substr(0, endOfDescription);
	getProteinAndDescr(descriptionTemp);
	
	
	//add spectrum count for *this protein to colname
	col[colIndex].count = elems[2];
	
	return true;
}

inline void Protein::getProteinAndDescr(string str)
{
	size_t firstSpace = str.find(" ");
	
	if(firstSpace == string::npos)
		protein = str;
	else{
		protein = str.substr(0, firstSpace);
		description = str.substr(firstSpace + 1);
	}
}

void Protein::addLoc(string newLoc)
{
	loc = newLoc;
}

void Protein::calcMW(const MWDB& mwDB)
{
	string sequence = mwDB.seqDB->getSequence(ID);

	if(sequence == SEQ_NOT_FOUND)
	{
		avgMass = -1;
		monoMass = -1;
		calcSequence = SEQ_NOT_FOUND;
		return;
	}
	
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
	locDB = new hashTable::HashTable<DBProtein>;
	proteins = new vector<Protein>;
}

Proteins::~Proteins()
{
	delete locDB;
	delete proteins;
}

//initialize Proteins.colnames with files contained in params file
Proteins::Proteins(const FilterFileParams& files)
{
	colIndex = 0;
	locDB = new hashTable::HashTable<DBProtein>;
	proteins = new vector<Protein>;
	
	for (int i = 0; i < files.numFiles; i++)
		colNames.push_back(files.file[i].colname);
}

long Proteins::insert(const Protein& newProtein)
{
	if(proteins->size() == 0)
	{
		proteins->push_back(newProtein);
		return 0;
	}
	else {
		long index = util::binSearch(proteins, newProtein, 0, proteins->size() - 1);
		if(index == -1)
		{
			vector<Protein>::iterator it = util::insertSorted(proteins, newProtein);
			return int(distance(proteins->begin(), it));
		}
		else {
			proteins->at(index).consolidate(newProtein, colIndex);
			return index;
		}
	}
}

//loop through protein headder and peptide lines in DTA filter file and add data to Proteins
bool Proteins::readIn(string wd, const FilterFileParam& filterFile, bool countUniquePeptides)
{
	string fname = wd + filterFile.path;
	ifstream inF(fname.c_str());
	if(!inF)
		return false;
	
	int numUniquePeptides = 0;
	long uniquePeptidesIndex = -1;
	bool inProtein = false;
	bool getNewLine = true;
	Protein blank, newProtein;
	blank.initialize(colNames);
	string line;
	
	while(!inF.eof()){
		if(getNewLine)
			getline(inF, line);
		getNewLine = true;
		if(util::strContains('%', line))  //find protein header lines by percent symbol for percent coverage
		{
			newProtein.initialize(colNames);
			if(newProtein.getProteinData(line, colIndex)) //if line is not column header line populate protein to Proteins
			{
				inProtein = true; //used to determine if whether it is valid to loop through peptide lines below protein
								  //header line to extract unique peptide spectral counts
				
				uniquePeptidesIndex = insert(newProtein); //insert new protein into proteins list
			}
			//extract spectral counts for unique peptides
			if(countUniquePeptides && inProtein)
			{
				do{
					getline(inF, line);
					if(line[0] == '*')
						numUniquePeptides += parsePeptideSC(line);
				} while(!util::strContains('%', line) && !inF.eof());
				proteins->at(uniquePeptidesIndex).col[colIndex].uniquePeptides = util::toString(numUniquePeptides);
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
	reverse(proteins->begin(), proteins->end()); //reverse proteins list so that reverse matches won't be in begining of dataset
	return true;
}

bool Proteins::readInLocDB(string fname)
{
	ifstream inF (fname.c_str());
	
	if(!inF)
		return false;
	
	string line;
	
	while(!inF.eof()){
		getline(inF, line);
		DBProtein newDBProtein(line);
		if (newDBProtein.ID != "ID") //skip line if it is header line
			locDB->insert(newDBProtein, newDBProtein.ID);
	}
	return true;
}

//write out combined protein lists to ofname in wide format
bool Proteins::writeOut(string ofname, const FilterFileParams& filterFileParams) const
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	bool parseSampleName = filterFileParams.sampleNamePrefix != "";
	int colNamesLength = DEFAULT_COL_NAMES_LENGTH;
	if(filterFileParams.getSubCelluarLoc)
		colNamesLength++;
	vector<string> headers;
	vector<string>::iterator it = headers.begin();
	headers.insert(it, DEFAULT_COL_NAMES, DEFAULT_COL_NAMES + colNamesLength);
	
	//print column headers
	if(filterFileParams.calcMW)
	{
		int mwHeadersLen = MWCALC_HEADERS_LENGTH;
		if(filterFileParams.includeSeq)
			mwHeadersLen++;
		
		it = headers.begin();
		for(int i = 0; i < colNamesLength; i++)
			if(headers[i] == "Mass(Da)")
			{
				headers.insert(it + i + 1, MWCALC_HEADERS, MWCALC_HEADERS + mwHeadersLen);
				colNamesLength = int(headers.size());
				break;
			}
	}
	else if(filterFileParams.includeSeq && !filterFileParams.calcMW)
	{
		it = headers.begin();
		for(int i = 0; i < headers.size(); i++)
			if(headers[i] == "Mass(Da)")
			{
				headers.insert(it + i + 1, MWCALC_HEADERS[2]);
				break;
			}
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
				outF << parseSample(colNames[i], filterFileParams.sampleNamePrefix, 0);
			else outF << '\t' << '\t' << parseSample(colNames[i], filterFileParams.sampleNamePrefix, 0);
		}
		outF << endl;
		for (int i = 0; i < colNamesLength; i++)
			outF << headers[i] <<'\t';
		for (int i = 0; i < colNames.size(); i++)
			outF << UNIQUE_PEPTIDE_HEADERS[0] << '\t' << UNIQUE_PEPTIDE_HEADERS[1] << '\t';
		outF << endl;
	}
	else
	{
		vector<string> ofColNames;
		int len = int(headers.size());
		for (int i = 0; i < len; i ++)
			ofColNames.push_back(headers[i]);
		for (int i = 0; i < colNames.size(); i ++)
			ofColNames.push_back(parseSample(colNames[i], filterFileParams.sampleNamePrefix, 0));
		int colNamesLen = int(ofColNames.size());
		for (int i = 0; i < colNamesLen; i++)
			outF << ofColNames[i] << '\t';
		outF << endl;
	}
	
	//print proteins and spectral counts
	int proteinsLen = int(proteins->size());
	for (int i = 0; i < proteinsLen; i++)
	{
		if (INCLUDE_FULL_DESCRIPTION)
			outF << proteins->at(i).fullDescription << '\t';
		
		outF << proteins->at(i).ID << '\t' <<
		proteins->at(i).protein << '\t' <<
		proteins->at(i).description << '\t' <<
		proteins->at(i).MW << '\t';
		
		if(filterFileParams.calcMW)
			outF << proteins->at(i).avgMass << '\t' <<
			proteins->at(i).monoMass << '\t';
			
		if(filterFileParams.includeSeq)
			outF << proteins->at(i).calcSequence << '\t';
		
		
		if(filterFileParams.getSubCelluarLoc)
			outF << proteins->at(i).loc << '\t';
		
		for (int j = 0; j < colIndex; j++)
		{
			outF << proteins->at(i).col[j].count << '\t';
			if (filterFileParams.includeUnique)
				outF << proteins->at(i).col[j].uniquePeptides << '\t';
		}
		
		outF << endl;
	}
	
	return true;
}

//write out combined protein lists to ofname in long format
bool Proteins::writeOutDB(string ofname, const FilterFileParams& filterFileParams) const
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	bool parseSampleName = (filterFileParams.sampleNamePrefix != "");
	
	//print column headers
	vector<string> headers;
	vector<string>::iterator it = headers.begin();
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
	
	if(filterFileParams.calcMW)
	{
		int mwHeadersLen = MWCALC_HEADERS_LENGTH;
		if(filterFileParams.includeSeq)
			mwHeadersLen++;
		
		it = headers.begin();
		for(int i = 0; i < headers.size(); i++)
			if(headers[i] == "Mass(Da)")
			{
				headers.insert(it + i + 1, MWCALC_HEADERS, MWCALC_HEADERS + mwHeadersLen);
				len = int(headers.size());
				break;
			}
	}
	else if(filterFileParams.includeSeq && !filterFileParams.calcMW)
	{
		it = headers.begin();
		for(int i = 0; i < headers.size(); i++)
			if(headers[i] == "Mass(Da)")
			{
				headers.insert(it + i + 1, MWCALC_HEADERS[2]);
				len = int(headers.size());
				break;
			}
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
	int proteinsLen = int(proteins->size());
	for (int i = 0; i < colNames.size(); i++)
	{
		for (int j = 0; j < proteinsLen; j++)
		{
			if (INCLUDE_FULL_DESCRIPTION)
				outF << proteins->at(j).fullDescription << '\t';
			
			outF << proteins->at(j).ID << '\t' <<
			proteins->at(j).protein << '\t' <<
			proteins->at(j).description << '\t' <<
			proteins->at(j).MW << '\t';
			
			if(filterFileParams.calcMW)
				outF << proteins->at(j).avgMass << '\t' <<
				proteins->at(j).monoMass << '\t';
			
			if(filterFileParams.includeSeq)
				outF << proteins->at(j).calcSequence << '\t';
			
			if(filterFileParams.getSubCelluarLoc)
				outF << proteins->at(j).loc << '\t';
			
			outF << proteins->at(j).col[i].colname << '\t' <<
			proteins->at(j).col[i].count;
			
			if (parseSampleName)
			{
				outF << '\t' << parseSample(proteins->at(j).col[i].colname, filterFileParams.sampleNamePrefix, 1) << '\t' <<
				parseReplicate(proteins->at(j).col[i].colname);
			}
			
			if (filterFileParams.includeUnique)
				outF << '\t' << proteins->at(j).col[i].uniquePeptides;
			
			outF << endl;
		}
	}
	return true;
}

//seearches binary tree containing subcelluar localization data and populates data
//to loc element for each Protein in Proteins
void Proteins::addSubcelluarLoc()
{
	hashTable::Node<DBProtein>* nodeTemp;
	
	int len = int(proteins->size());
	for (int i = 0; i < len; i++)
	{
		nodeTemp = locDB->getItem(proteins->at(i).ID);
		nodeTemp == nullptr ? proteins->at(i).loc = string(LOC_NOT_FOUND) : proteins->at(i).loc = string(nodeTemp->val.loc);
	}
}

void Proteins::calcMW(const MWDB& mwDB)
{
	int len = int(proteins->size());
	for (int i = 0; i < len; i++)
		proteins->at(i).calcMW(mwDB);
}

void Proteins::addSeq(const SeqDB& seqDB)
{
	int len = int(proteins->size());
	for(int i = 0; i < len; i++)
		proteins->at(i).calcSequence = seqDB.getSequence(proteins->at(i).ID);
}

//check if line containing % is a collumn header line instead of a protein header line
bool isColumnHeaderLine(const vector<string>& elems)
{
	for (int i = 0; i < COLUMN_HEADER_LINE_ELEMENTS_LENGTH; i++)
		if(COLUMN_HEADER_LINE_ELEMENTS[i] != elems[i])
			return false;
	
	return true;
}

//optional fxn to parse long sample name
string parseSample(string sampleName, string prefix, bool outputFormat)
{
	//return unparsed sampleName if prefix is empty string or is not found in sampleName
	if(!util::strContains(prefix, sampleName) || prefix.length() == 0)
		return sampleName;
	
	string sample = util::removeSubstr(prefix, sampleName); //remove prefix from sampleName
	
	return outputFormat ? sample.substr(0, sample.find_last_of("_")) : sample;
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
	util::split(line, '\t', elems);
	assert(elems.size() > 11);
	
	//return SC for peptide as int
	return util::toInt(elems[11]);
}

string getID(string str)
{
	size_t firstBar = str.find("|");
	size_t secBar = str.find("|", firstBar+1);
	return str.substr(firstBar+1, secBar-firstBar-1);
}
