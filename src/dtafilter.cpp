
#include "DTarray_AJM.hpp"

inline bool ProteinTemplate::operator == (const ProteinTemplate& comp) const
{
	return comp.ID == ID;
}

inline bool ProteinTemplate::operator == (string comp) const
{
	return comp == ID;
}

inline bool ProteinTemplate::operator > (const ProteinTemplate& comp) const
{
	return comp.ID > ID;
}

inline bool ProteinTemplate::operator < (const ProteinTemplate& comp) const
{
	return comp.ID < ID;
}

template<class T>
long DBTemplate<T>::insert(const T& newData)
{
	if(data->size() == 0)
	{
		data->push_back(newData);
		return 0;
	}
	else {
		long index = util::binSearch(data, newData, 0, data->size() - 1);
		if(index == -1)
		{
			vector<Protein>::iterator it = util::insertSorted(data, newData);
			return int(distance(data->begin(), it));
		}
		else {
			data->at(index).consolidate(newData, colIndex);
			return index;
		}
	}
}

template <class T>
void DBTemplate<T>::calcMW(const mwDB::MWDB& mwDB)
{
	int len = int(data->size());
	for (int i = 0; i < len; i++)
		data->at(i).calcMW(mwDB, data->at(i).get_ID());
}

void ProteinDataTemplate::consolidate(const ProteinDataTemplate& toAdd, int colIndex)
{
	col[colIndex].colname = toAdd.col[colIndex].colname;
	col[colIndex].count = toAdd.col[colIndex].count;
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

void ProteinDataTemplate::calcMW(const mwDB::MWDB& mwDB, string ID)
{
	string tempSequence = mwDB.seqDB->getSequence(ID);

	if(sequence == SEQ_NOT_FOUND)
	{
		sequence = SEQ_NOT_FOUND;
		avgMass = -1;
		monoMass = -1;
		return;
	}
	
	avgMass = mwDB.calcMW(sequence, 0);
	monoMass = mwDB.calcMW(sequence, 1);
	sequence = tempSequence;
}

//convert Protein to DBprotein to allow inteface with BinTree of DBproteins
DBProtein Protein::toDBprotein() const
{
	DBProtein dbprotein;
	dbprotein.ID = ID;
	
	return dbprotein;
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
	
	bool extractPeptides = true;
	
	while(!inF.eof()){
		if(getNewLine)
			getline(inF, line);
		getNewLine = true;
		if(util::strContains('%', line))  //find protein header lines by percent symbol for percent coverage
		{
			newProtein.initialize(colNames);
			if(newProtein.getProteinData(line, colIndex)) //if line is not column header line populate protein to Proteins
			{
				inProtein = true; //used to determine if it is valid to loop through peptide lines below protein
								  //header line to extract unique peptide spectral counts
				
				uniquePeptidesIndex = insert(newProtein); //insert new protein into proteins list
			}
			//extract spectral counts for unique peptides
			if((countUniquePeptides || extractPeptides ) && inProtein)
			{
				if(countUniquePeptides)
				{
				   do{
					getline(inF, line);
					if(line[0] == '*')
						numUniquePeptides += parsePeptideSC(line);
				} while(!util::strContains('%', line) && !inF.eof());
				data->at(uniquePeptidesIndex).col[colIndex].uniquePeptides = util::toString(numUniquePeptides);
				numUniquePeptides = 0;
				getNewLine = false;
				inProtein = false;
				}
				
				/*if(extractPeptides)
				{
					do{
						getline(inF, line);
						
						string temp = line;
						
						parsePeptide(line);
						
					} while(!util::strContains('%', line) && !inF.eof());*/
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
	reverse(data->begin(), data->end()); //reverse proteins list so that reverse matches won't be in begining of dataset
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
	int dataLen = int(data->size());
	for (int i = 0; i < dataLen; i++)
	{
		if (INCLUDE_FULL_DESCRIPTION)
			outF << data->at(i).fullDescription << '\t';
		
		outF << data->at(i).ID << '\t' <<
		data->at(i).protein << '\t' <<
		data->at(i).description << '\t' <<
		data->at(i).MW << '\t';
		
		if(filterFileParams.calcMW)
			outF << data->at(i).avgMass << '\t' <<
			data->at(i).monoMass << '\t';
			
		if(filterFileParams.includeSeq)
			outF << data->at(i).sequence << '\t';
		
		
		if(filterFileParams.getSubCelluarLoc)
			outF << data->at(i).loc << '\t';
		
		for (int j = 0; j < colIndex; j++)
		{
			outF << data->at(i).col[j].count << '\t';
			if (filterFileParams.includeUnique)
				outF << data->at(i).col[j].uniquePeptides << '\t';
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
	int dataLen = int(data->size());
	for (int i = 0; i < colNames.size(); i++)
	{
		for (int j = 0; j < dataLen; j++)
		{
			if (INCLUDE_FULL_DESCRIPTION)
				outF << data->at(j).fullDescription << '\t';
			
			outF << data->at(j).ID << '\t' <<
			data->at(j).protein << '\t' <<
			data->at(j).description << '\t' <<
			data->at(j).MW;
			
			if(filterFileParams.calcMW)
				outF << '\t' << data->at(j).avgMass << '\t' <<
				data->at(j).monoMass << '\t';
			
			if(filterFileParams.includeSeq)
				outF << '\t' << data->at(j).sequence << '\t';
			
			if(filterFileParams.getSubCelluarLoc)
				outF << '\t' << data->at(j).loc << '\t';
			
			outF << '\t' << data->at(j).col[i].colname << '\t' <<
			data->at(j).col[i].count;
			
			if (parseSampleName)
			{
				outF << '\t' << parseSample(data->at(j).col[i].colname, filterFileParams.sampleNamePrefix, 1) << '\t' <<
				parseReplicate(data->at(j).col[i].colname);
			}
			
			if (filterFileParams.includeUnique)
				outF << '\t' << data->at(j).col[i].uniquePeptides;
			
			outF << endl;
		}
	}
	return true;
}

//seearches binary tree containing subcelluar localization data and populates data
//to loc element for each Protein in data
void Proteins::addSubcelluarLoc()
{
	hashTable::Node<DBProtein>* nodeTemp;
	
	int len = int(data->size());
	for (int i = 0; i < len; i++)
	{
		nodeTemp = locDB->getItem(data->at(i).ID);
		nodeTemp == nullptr ? data->at(i).loc = string(LOC_NOT_FOUND) : data->at(i).loc = string(nodeTemp->val.loc);
	}
}

void Proteins::addSeq(const mwDB::SeqDB& seqDB)
{
	int len = int(data->size());
	for(int i = 0; i < len; i++)
		data->at(i).sequence = seqDB.getSequence(data->at(i).ID);
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
