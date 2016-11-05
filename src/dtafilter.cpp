
#include "DTarray_AJM.hpp"
#include "baseClasses.hpp"

void Protein::consolidate(const Protein& toAdd)
{
	col[*colIndex].colname = toAdd.col[*colIndex].colname;
	col[*colIndex].count = toAdd.col[*colIndex].count;
	col[*colIndex].coverage = toAdd.col[*colIndex].coverage;
}

void Peptide::consolidate(const Peptide& toAdd)
{
	col[*colIndex].colname = toAdd.col[*colIndex].colname;
	col[*colIndex].count = toAdd.col[*colIndex].count;
	col[*colIndex].scan = toAdd.col[*colIndex].scan;
	col[*colIndex].parentFile = toAdd.col[*colIndex].parentFile;
	col[*colIndex].obsMH = toAdd.col[*colIndex].obsMH;
}

//fill col element with colNames found in params file
void Protein::initialize(FilterFileData_protein* const tempColNames, size_t colNamesLen, size_t* _colIndex)
{
	//free(col);
	col = (FilterFileData_protein*) calloc(colNamesLen, colNamesLen * sizeof(FilterFileData_protein));
	memcpy(col, tempColNames, colNamesLen * sizeof(FilterFileData_protein));
	colIndex = _colIndex;
	colSize = colNamesLen;
}

void Peptide::initialize(FilterFileData_peptide* const tempColNames, size_t colNamesLen, size_t* _colIndex)
{
	//free(col);
	col = (FilterFileData_peptide*) calloc(colNamesLen, colNamesLen * sizeof(FilterFileData_peptide));
	memcpy(col, tempColNames, colNamesLen * sizeof(FilterFileData_peptide));
	colIndex = _colIndex;
	colSize = colNamesLen;
}

void Protein::addSeq()
{
	assert(seqDB != nullptr);
	sequence = seqDB->getSequence(ID);
}

void Protein::addLoc()
{
	assert(locDB != nullptr);
	loc = locDB->getDat(ID);
}

void Protein::addFxn()
{
	assert(fxnDB != nullptr);
	fxn = fxnDB->getDat(ID);
}

//parse proten header line and extract desired data
bool Protein::getProteinData(string line, size_t colIndex)
{
	//split line by tabs
	vector<string> elems;
	util::split(line, IN_DELIM, elems);
	
	//tell readIn function to skip line if it is header line
	if(isColumnHeaderLine(elems))
		return false;
	
	//keep fullDescription but seperate by spaces instead of tabs
	fullDescription = elems[0];
	int len = int(elems.size());
	for (int i = 1; i < len; i++)
		fullDescription += (" " + elems[i]);
	
	//extract sequence count
	sequenceCount = util::toInt(elems[1]);
	
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
	
	//add spectrum count and coverage for *this protein to colname
	col[colIndex].count = elems[2];
	string coverageTemp = elems[3];
	coverageTemp = coverageTemp.substr(0, coverageTemp.find("%"));
	col[colIndex].coverage = coverageTemp;
	
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

void Protein::calcMW()
{
	string tmp_sequence = mwdb->seqDB->getSequence(ID);
	
	if(tmp_sequence == SEQ_NOT_FOUND)
	{
		sequence = SEQ_NOT_FOUND;
		avgMass = -1;
		monoMass = -1;
		return;
	}
	
	sequence = tmp_sequence;
	avgMass = mwdb->calcMW(sequence, 0);
	monoMass = mwdb->calcMW(sequence, 1);
}

void Peptide::calcMW()
{
	avgMass = mwdb->calcMW(calcSequence, 0);
	monoMass = mwdb->calcMW(calcSequence, 1);
}

//loop through protein headder and peptide lines in DTA filter file and add data to Proteins
bool Proteins::readIn(string wd, FilterFileParams* const pars, const FilterFileParam& filterFile,
					  FilterFileData_protein* const colNamesTemp, FilterFileData_peptide* const pColNamesTemp, Peptides * const peptides)
{
	string fname = wd + filterFile.path;
	ifstream inF(fname.c_str());
	if(!inF)
		return false;
	
	int numUniquePeptides = 0;
	hashTable::Node<Protein>* uniquePeptidesIndex = nullptr;
	bool inProtein = false;
	bool getNewLine = true;
	Protein newProtein(pars, locDB, fxnDB, mwdb, seqDB);
	Peptide newPeptide(pars, mwdb);
	string line;
	hashTable::Node<Peptide>* peptidesIndex = nullptr;
	int end = 0;
	size_t colNamesLen = pars->numFiles;
	
	while(!inF.eof()){
		if(getNewLine)
			getline(inF, line);
		getNewLine = true;
		if(util::strContains('%', line))  //find protein header lines by percent symbol for percent coverage
		{
			newProtein.clear();
			newProtein.initialize(colNamesTemp, colNamesLen, &colIndex);
			if(newProtein.getProteinData(line, colIndex)) //if line is not column header line populate protein to Proteins
			{
				end = newProtein.sequenceCount;
				inProtein = true; //used to determine if it is valid to loop through peptide lines below protein
								  //header line to extract unique peptide spectral counts
				
				uniquePeptidesIndex = data->consolidate(newProtein, newProtein.ID); //insert new protein into proteins list
			}
			//extract spectral counts for unique peptides
			if((pars->includeUnique || pars->includePeptides ) && inProtein)
			{
				getline(inF, line);
				getNewLine = false;
				for(int i = 0; i < end && !inF.eof() && !util::strContains('%', line); i++)
				{
					if(getNewLine)
						getline(inF, line);
					getNewLine = true;
					
					if(util::strContains('%', line))
						throw runtime_error("Bad filter file: " + filterFile.colname);

					if(pars->includePeptides)
					{
						newPeptide.initialize(pColNamesTemp, colNamesLen, &peptides->colIndex);
						
						newPeptide.proteinID = newProtein.ID;
						newPeptide.protein = newProtein.protein;
						newPeptide.description = newProtein.description;
						newPeptide.parsePeptide(line);
						peptidesIndex = peptides->data->consolidate(newPeptide, newPeptide.key);
					}
					if(pars->includeUnique)
					{
						if(line[0] == '*')
							numUniquePeptides += parsePeptideSC(line);
					}
				}
				uniquePeptidesIndex->val.col[colIndex].uniquePeptides = util::toString(numUniquePeptides);
				numUniquePeptides = 0;
				getNewLine = false;
				inProtein = false;
			}
		}
	}
	colIndex++;
	if(pars->includePeptides)
		peptides->colIndex++;
	return true;
}

void Peptide::parsePeptide(const string& line)
{
	vector<string> elems;
	util::split(line, IN_DELIM, elems);
	
	unique = elems[0] == "*";
	calcMH = elems[5];
	
	calcSequence = parseSequence(elems[12]);
	sequence = elems[12];
	col[*colIndex].count = elems[11];
	col[*colIndex].obsMH = elems[5];
	
	calcMH = elems[6];
	
	string temp = elems[1];
	elems.clear();
	util::split(temp, '.', elems);
	scan = elems[2];
	charge = elems[3];
	key = makeKey();
	col[*colIndex].scan = elems[2];
	col[*colIndex].parentFile = elems[0];
}

inline void Protein::clear()
{
	MW.clear();
	fullDescription.clear();
	matchDirrection.clear();
	sequenceCount = 0;
	ID.clear();
	protein.clear();
	description.clear();
	loc.clear();
	avgMass = 0;
	monoMass = 0;
	sequence.clear();
}

inline void Peptide::clear()
{
	unique = 0;
	calcMH.clear();
	key.clear();
	scan.clear();
	avgMass = 0;
	monoMass = 0;
	proteinID.clear();
}

//public Proteins::readIn function which adds all files in filterFileParams to Proteins
//and summarizes progress for user.
bool Proteins::readIn(string wd, FilterFileParams& filterFileParams, Peptides * const peptides)
{
	size_t colNamesLen = filterFileParams.numFiles;
	FilterFileData_protein* colNamesTemp = (FilterFileData_protein*) calloc(colNamesLen, sizeof(FilterFileData_protein));
	FilterFileData_peptide* pColNamesTemp = (FilterFileData_peptide*) calloc(colNamesLen, sizeof(FilterFileData_peptide));
	
	for(int i = 0; i < colNamesLen; i++)
	{
		FilterFileData_protein temp (colNames[i]);
		*(colNamesTemp + i) = temp;
		FilterFileData_peptide pTemp (colNames[i]);
		*(pColNamesTemp + i) = pTemp;
	}
	
	for (int i = 0; i < filterFileParams.numFiles; i++)
	{
		if(!readIn(wd, &filterFileParams, filterFileParams.file[i], colNamesTemp, pColNamesTemp, peptides))
		{
			cout <<"Failed to read in " << filterFileParams.file[i].path <<"!" << endl << "Exiting..." << endl;
			return false;
		}
		cout << "Adding " << filterFileParams.file[i].colname << "..." << endl;
	}
	return true;
}

bool Proteins::readInMWdb(string wd, const FilterFileParams& par)
{
	mwdb = new mwDB::MWDB_Protein;
	return mwdb->readIn(wd, par);
}

bool Proteins::readInSeqDB(string fname)
{
	seqDB = new mwDB::SeqDB;
	return seqDB->readIn(fname);
}

bool Proteins::readInFxnDB(string fname)
{
	fxnDB = new Dbase;
	return fxnDB->readIn(fname);
}

bool Proteins::readInLocDB(string fname)
{
	locDB = new Dbase;
	return locDB->readIn(fname);
}

void Protein::write(ofstream& outF)
{
	if(!outF)
		throw runtime_error("Bad ofstream");
	
	if(!supDataAdded)
	{
		if(par->calcMW)
			calcMW();
		if(par->getSubCelluarLoc)
			addLoc();
		if(par->getSeq && !par->calcMW)
			addSeq();
		if(par->getFxn)
			addFxn();
		supDataAdded = true;
	}
	
	if (INCLUDE_FULL_DESCRIPTION)
		outF << fullDescription << OUT_DELIM;
	
	outF << ID <<
	OUT_DELIM << protein <<
	OUT_DELIM << description <<
	OUT_DELIM << MW;
	
	if(par->getFxn)
		outF << OUT_DELIM << fxn;
	
	if(par->calcMW)
		outF << OUT_DELIM << avgMass <<
		OUT_DELIM << monoMass;
	
	if(par->getSeq)
		outF << OUT_DELIM << sequence;
	
	if(par->getSubCelluarLoc)
		outF << OUT_DELIM <<  loc;
	
	if(par->outputFormat == 1)
	{
		for(int i = 0; i < colSize; i++)
		{
			outF << OUT_DELIM << col[i].count;
			if (par->includeUnique)
				outF << OUT_DELIM << col[i].uniquePeptides;
			
			if(par->includeCoverage)
				outF << OUT_DELIM << col[i].coverage;
		}
	}
	else if(par->outputFormat == 2)
	{
		if(par->includeCoverage)
			outF << OUT_DELIM << col[*colIndex].coverage;
		
		outF << OUT_DELIM << col[*colIndex].colname <<
		OUT_DELIM << col[*colIndex].count;
		
		if (!par->sampleNamePrefix.empty())
		{
			outF << OUT_DELIM << parseSample(col[*colIndex].colname, par->sampleNamePrefix, 1) << OUT_DELIM <<
			parseReplicate(col[*colIndex].colname);
		}
		
		if (par->includeUnique)
			outF << OUT_DELIM << col[*colIndex].uniquePeptides;
	}
	outF << endl;
}

//write out combined protein lists to ofname in wide format
bool Proteins::writeOut(string ofname, const FilterFileParams& filterFileParams)
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	int outputFormat = filterFileParams.outputFormat;
	filterFileParams.outputFormat = 1;
	
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
		if(filterFileParams.getSeq)
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
	else if(filterFileParams.getSeq && !filterFileParams.calcMW)
	{
		it = headers.begin();
		for(int i = 0; i < headers.size(); i++)
			if(headers[i] == "Mass(Da)")
			{
				headers.insert(it + i + 1, MWCALC_HEADERS[2]);
				break;
			}
	}
	if(filterFileParams.getFxn)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Mass(Da)")
			{
				headers.insert(it + 1, "Function");
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
			outF << OUT_DELIM;
		for (int i = 0; i < filterFileParams.numFiles ; i++)
			outF << colNames[i] << delim;
		outF << endl;
	}
	if (filterFileParams.includeUnique || filterFileParams.includeCoverage)
	{
		int numCols = filterFileParams.includeUnique + filterFileParams.includeCoverage;
		string tabs = "\t";
		
		for(int i = 0; i < numCols; i++)
			tabs += OUT_DELIM;
		
		for (int i = 0; i < colNamesLength; i++)
			outF << OUT_DELIM;
		for (int i = 0; i < filterFileParams.numFiles; i++)
		{
			if (i == 0)
				outF << parseSample(colNames[i], filterFileParams.sampleNamePrefix, 0);
			else outF << tabs << parseSample(colNames[i], filterFileParams.sampleNamePrefix, 0);
		}
		outF << endl;
		for (int i = 0; i < colNamesLength; i++)
			outF << headers[i] << OUT_DELIM;
		
		string repeatedHeaders = "";
		if(filterFileParams.includeUnique)
			repeatedHeaders += (SUP_INFO_HEADERS[0] + OUT_DELIM + SUP_INFO_HEADERS[1]);
		if(filterFileParams.includeCoverage && filterFileParams.includeUnique)
			repeatedHeaders += (OUT_DELIM + SUP_INFO_HEADERS[2]);
		else if(filterFileParams.includeCoverage)
			repeatedHeaders += (SUP_INFO_HEADERS[0] + OUT_DELIM + SUP_INFO_HEADERS[2]);
		
		for (int i = 0; i < filterFileParams.numFiles; i++)
		{
			if(i == 0)
				outF << repeatedHeaders;
			else outF << OUT_DELIM << repeatedHeaders;
		}
		outF << endl;
	}
	else
	{
		vector<string> ofColNames;
		int len = int(headers.size());
		for (int i = 0; i < len; i ++)
			ofColNames.push_back(headers[i]);
		for (int i = 0; i < filterFileParams.numFiles; i ++)
			ofColNames.push_back(parseSample(colNames[i], filterFileParams.sampleNamePrefix, 0));
		int colNamesLen = int(ofColNames.size());
		for (int i = 0; i < colNamesLen; i++)
			outF << ofColNames[i] << OUT_DELIM;
		outF << endl;
	}
	
	//print proteins and spectral counts
	data->write(outF);
	
	filterFileParams.outputFormat = outputFormat;
	
	return true;
}

//write out combined protein lists to ofname in long format
bool Proteins::writeOutDB(string ofname, const FilterFileParams& filterFileParams)
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	int outputFormat = filterFileParams.outputFormat;
	filterFileParams.outputFormat = 2;

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
		if(filterFileParams.getSeq)
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
	else if(filterFileParams.getSeq && !filterFileParams.calcMW)
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
	if(filterFileParams.getFxn)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Mass(Da)")
			{
				headers.insert(it + 1, "Function");
				len = int(headers.size());
				break;
			}
	}
	
	if(filterFileParams.includeCoverage)
	{
		for(it = headers.begin(); it != headers.end(); it++)
		{
			if(*it == "Long_sample_name")
			{
				headers.insert(it, "percent_coverage");
				len = int(headers.size());
				break;
			}
		}
	}
	
	vector<string> ofColNames;
	for (int i = 0; i < len - (!parseSampleName * 2); i ++)
		ofColNames.push_back(headers[i]);
	if(filterFileParams.includeUnique)
		ofColNames.push_back(SUP_INFO_HEADERS[1]);
	
	for (int i = 0; i < int(ofColNames.size()); i++)
		outF << ofColNames[i] << OUT_DELIM;
	outF << endl;
	
	//print proteins and spectral counts
	for (colIndex = 0; colIndex < filterFileParams.numFiles; colIndex++)
	{
		data->write(outF);
	}
	
	filterFileParams.outputFormat = outputFormat;
	
	return true;
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
	util::split(line, IN_DELIM, elems);
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

inline string parseSequence(string str)
{
	size_t firstP = str.find(".");
	size_t secP = str.find_last_of(".");

	return (firstP == string::npos || secP == string::npos) ? str : str.substr(firstP + 1, secP - (str.length() - secP));
}

void Peptide::write(ofstream& outF)
{
	if(!outF)
		throw runtime_error("Bad ofstream!");
	
	if(!supDataAdded)
	{
		if(par->calcMW)
			calcMW();
		supDataAdded = true;
	}
	
	outF << proteinID << OUT_DELIM <<
	protein << OUT_DELIM <<
	description << OUT_DELIM <<
	sequence << OUT_DELIM <<
	charge << OUT_DELIM <<
	unique << OUT_DELIM <<
	calcMH << OUT_DELIM;
	
	if(par->calcMW)
		outF << avgMass << OUT_DELIM <<
		monoMass << OUT_DELIM;
	
	if(par->outputFormat == 1)
	{
		for (int i = 0; i < colSize; i++)
		{
			if(i == 0)
				outF << col[i].count;
			else outF << OUT_DELIM << col[i].count;
		}
	}
	else if(par->outputFormat == 2)
	{
		outF << col[*colIndex].obsMH << OUT_DELIM <<
		col[*colIndex].scan << OUT_DELIM <<
		col[*colIndex].parentFile << OUT_DELIM <<
		col[*colIndex].colname << OUT_DELIM <<
		col[*colIndex].count;
		
		if(!par->sampleNamePrefix.empty())
		{
			outF << OUT_DELIM << parseSample(col[*colIndex].colname, par->sampleNamePrefix, 1)
			<< OUT_DELIM <<	parseReplicate(col[*colIndex].colname);
		}
	}
	
	outF << endl;
}

bool Peptides::writeOut(string ofname, const FilterFileParams& pars)
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	int outputFormat = pars.outputFormat;
	pars.outputFormat = 1;
	
	//generate headers based off params
	vector<string> headers;
	vector<string>::iterator it;
	for(int i = 0; i < DEFALUT_PEPTIDE_COLNAMES_LEN; i++)
		headers.push_back(DEFALUT_PEPTIDE_COLNAMES[i]);
	if(pars.calcMW)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "calcMH")
			{
				headers.insert(it + 1, MWCALC_HEADERS, MWCALC_HEADERS + 2);
				break;
			}
	}
	headers.insert(headers.end(), colNames, colNames + (pars.numFiles));
	
	//print headers
	for(it = headers.begin(); it != headers.end(); it++)
	{
		if(it == headers.begin())
			outF << *it;
		else outF << OUT_DELIM << *it;
	}
	outF << endl;
	
	//print peptides and spectral counts
	data->write(outF);
	
	pars.outputFormat = outputFormat;
	
	return true;
}


bool Peptides::writeOutDB(string ofname, const FilterFileParams& pars)
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	int outputFormat = pars.outputFormat;
	pars.outputFormat = 2;
	
	//generate headers based off params
	vector<string> headers;
	vector<string>::iterator it;
	int defaultColNamesLen = DEFALUT_PEPTIDE_DB_COLNAMES_LEN;
	
	if(!pars.sampleNamePrefix.empty())
		defaultColNamesLen += 2;
	
	for(int i = 0; i < defaultColNamesLen; i++)
		headers.push_back(DEFALUT_PEPTIDE_DB_COLNAMES[i]);
	if(pars.calcMW)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "sequence")
			{
				headers.insert(it, MWCALC_HEADERS, MWCALC_HEADERS + 2);
				break;
			}
	}
	
	//print headers
	for(it = headers.begin(); it != headers.end(); it++)
	{
		if(it == headers.begin())
			outF << *it;
		else outF << OUT_DELIM << *it;
	}
	outF << endl;

	for(colIndex = 0; colIndex < pars.numFiles; colIndex++)
		data->write(outF);
	
	pars.outputFormat = outputFormat;
	
	return true;
}


