//
//  dtafilter.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 3/25/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "dtafilter.hpp"

void Protein::consolidate(const Protein& toAdd)
{
	col[*toAdd.colIndex].colname = toAdd.col[*toAdd.colIndex].colname;
	col[*toAdd.colIndex].count = toAdd.col[*toAdd.colIndex].count;
	col[*toAdd.colIndex].coverage = toAdd.col[*toAdd.colIndex].coverage;
	col[*toAdd.colIndex].sequenceCount = toAdd.col[*toAdd.colIndex].sequenceCount;
}

void Peptide::consolidate(const Peptide& toAdd)
{
	col[*toAdd.colIndex].colname = toAdd.col[*toAdd.colIndex].colname;
	col[*toAdd.colIndex].count += toAdd.col[*toAdd.colIndex].count;
	col[*toAdd.colIndex].scan = toAdd.col[*toAdd.colIndex].scan;
	col[*toAdd.colIndex].parentFile = toAdd.col[*toAdd.colIndex].parentFile;
	col[*toAdd.colIndex].obsMH = toAdd.col[*toAdd.colIndex].obsMH;
}

template <class _Tp>
void ProteinDataTemplate<_Tp>::initialize(const vector<_Tp>& tempColNames, size_t colNamesLen, size_t* _colIndex)
{
	col.clear();
	col.insert(col.begin(), tempColNames.begin(), tempColNames.end());
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

inline void Protein::operator = (const Protein& _p)
{
	ID = _p.ID;
	protein = _p.protein;
	description = _p.description;
	MW = _p.MW;
	pI = _p.pI;
	loc = _p.loc;
	fxn = _p.fxn;
	fullDescription = _p.fullDescription;
	matchDirrection = _p.matchDirrection;
	col.insert(col.begin(), _p.col.begin(), _p.col.end());
	avgMass = _p.avgMass;
	monoMass = _p.monoMass;
	sequence = _p.sequence;
	supDataAdded = _p.supDataAdded;
}

inline void Peptide::operator = (const Peptide& _p)
{
	key = _p.key;
	calcSequence = _p.calcSequence;
	proteinID = _p.proteinID;
	calcMH = _p.calcMH;
	fileName = _p.fileName;
	protein = _p.protein;
	description = _p.description;
	charge = _p.charge;
	unique = _p.unique;
	col.insert(col.begin(), _p.col.begin(), _p.col.end());
	avgMass = _p.avgMass;
	monoMass = _p.monoMass;
	sequence = _p.sequence;
	supDataAdded = _p.supDataAdded;
}

//parse proten header line and extract desired data
void Protein::getProteinData(string line)
{
	//split line by tabs
	vector<string> elems;
	utils::split(line, IN_DELIM, elems);
	assert(elems.size() == 9);
	
	//keep fullDescription but seperate by spaces instead of tabs
	fullDescription = elems[0];
	int len = int(elems.size());
	for (int i = 1; i < len; i++)
		fullDescription += (" " + elems[i]);
	
	//extract pI
	pI = elems[6];
	
	//extract protein length
	length = elems[4];
	
	//extract matchDirrection
	size_t firstBar = elems[0].find("|");
	matchDirrection = line.substr(0, firstBar);
	
	//Extract uniprotID
	ID = ::getID(elems[0]);
	if(matchDirrection == "Reverse_sp")
		ID = "reverse_" + ID;
	
	//extract MW
	MW = elems[5];
	
	//extract shortened protein name and description
	size_t endOfDescription = elems[8].find(" [");
	description = elems[8].substr(0, endOfDescription);
	getProtein(description);
	
	//add spectrum count, coverage and sequence count for *this protein to colname
	col[*colIndex].count = utils::toInt(elems[2]);
	string coverageTemp = elems[3];
	coverageTemp = coverageTemp.substr(0, coverageTemp.find("%"));
	col[*colIndex].coverage = coverageTemp;
	col[*colIndex].sequenceCount = elems[1];
}

inline void Protein::getProtein(string str)
{
	size_t firstSpace = str.find(" ");
	
	if(firstSpace == string::npos)
		protein = str;
	else{
		protein = str.substr(0, firstSpace);
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
bool Proteins::readIn(string wd, filterFile::FilterFileParams* const pars,
					  const vector <filterFile::FilterFileData_protein>& colNamesTemp,
					  const vector<filterFile::FilterFileData_peptide>& pColNamesTemp,
					  Peptides * const peptides)
{
	string fname = wd + pars->getFilePath(colIndex);
	
	utils::File file;
	if(!file.read(fname))
		return false;
	
	int numUniquePeptides = 0;
	Protein* uniquePeptidesIndex = nullptr;
	bool inProtein = false;
	bool getNewLine = true;
	string line;
	Peptide* peptidesIndex = nullptr;
	size_t colNamesLen = pars->numFiles;
	
	while(!file.end()){
		if(getNewLine)
			line = file.getLine();
		getNewLine = true;
		if(utils::strContains('%', line))  //find protein header lines by percent symbol for percent coverage
		{
			if(isColumnHeaderLine(line)) //skip if header line
				continue;
			else{
				Protein newProtein(pars, locDB, fxnDB, mwdb, seqDB);
				newProtein.initialize(colNamesTemp, colNamesLen, &colIndex);
				newProtein.getProteinData(line);
				inProtein = true; //used to determine if it is valid to loop through peptide lines below protein
								  //header line to extract unique peptide spectral counts
					
				uniquePeptidesIndex = data->consolidate(newProtein, newProtein.ID); //insert new protein into proteins list
				
				//get peptide data if applicable
				if((pars->includeUnique || pars->includePeptides ) && inProtein)
				{
					numUniquePeptides = 0;
					line = file.getLine();
					getNewLine = false;
					while(!file.end())
					{
						Peptide newPeptide(pars, mwdb);
						if(getNewLine)
							line = file.getLine();
						getNewLine = true;
						
						//break if starting new protein or end of file
						if(utils::strContains('%', line) || line == "\tProteins\tPeptide IDs\tSpectra")
							break;
						
						if(pars->includePeptides)
						{
							newPeptide.initialize(pColNamesTemp, colNamesLen, &peptides->colIndex);
							newPeptide.proteinID = newProtein.ID;
							newPeptide.protein = newProtein.protein;
							newPeptide.description = newProtein.description;
							newPeptide.parsePeptide(line);
							
							if(pars->peptideGroupMethod == filterFile::byScan)
								peptides->data->insert(newPeptide, newPeptide.key);
							else peptidesIndex = peptides->data->consolidate(newPeptide, newPeptide.key);
						}
						if(pars->includeUnique)
						{
							if(line[0] == '*')
								numUniquePeptides += parsePeptideSC(line);
						}
					}
					uniquePeptidesIndex->col[colIndex].uniquePeptides = utils::toString(numUniquePeptides);
					numUniquePeptides = 0;
					getNewLine = false;
					inProtein = false;
				}
			}
		}
	}
	return true;
}

void Peptide::parsePeptide(const string& line)
{
	vector<string> elems;
	utils::split(line, IN_DELIM, elems);
	
	assert(elems.size() == 13);
	
	unique = elems[0] == "*";
	calcMH = elems[5];
	
	calcSequence = parseSequence(elems[12]);
	length = utils::toString(calcSequence.length());
	sequence = elems[12];
	col[*colIndex].count = utils::toInt(elems[11]);
	col[*colIndex].obsMH = elems[5];
	
	calcMH = elems[6];
	fileName = elems[1];
	
	string temp = elems[1];
	elems.clear();
	utils::split(temp, '.', elems);
	charge = elems[3];
	key = makeKey();
	col[*colIndex].scan = elems[2];
	col[*colIndex].parentFile = elems[0];
}

inline string Peptide::makeKey() const {
	switch(par->peptideGroupMethod){
		case filterFile::byScan : return fileName;
			break;
		case filterFile::byProtein : return proteinID + "_" + sequence + "_" + charge;
			break;
		case filterFile::byCharge : return proteinID + "_" + sequence;
			break;
		default : throw runtime_error("Invalid type");
	}
}

inline void Protein::clear()
{
	col.clear();
	MW.clear();
	fullDescription.clear();
	matchDirrection.clear();
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
	col.clear();
	unique = 0;
	calcMH.clear();
	key.clear();
	fileName.clear();
	avgMass = 0;
	monoMass = 0;
	proteinID.clear();
}

//public Proteins::readIn function which adds all files in filterFileParams to Proteins
//and summarizes progress for user.
bool Proteins::readIn(string wd, filterFile::FilterFileParams& par, Peptides * const peptides)
{
	size_t colNamesLen = par.numFiles;
	
	vector<filterFile::FilterFileData_protein> colNamesTemp;
	vector<filterFile::FilterFileData_peptide> pColNamesTemp;
	
	for(int i = 0; i < colNamesLen; i++)
	{
		filterFile::FilterFileData_protein temp (colNames[i]);
		colNamesTemp.push_back(temp);
		filterFile::FilterFileData_peptide pTemp (colNames[i]);
		pColNamesTemp.push_back(pTemp);
	}
	
	for(DBTemplate::colIndex = 0; DBTemplate::colIndex < par.numFiles; DBTemplate::colIndex++)
	{
		Peptides::colIndex = Proteins::colIndex;
		if(!readIn(wd, &par, colNamesTemp, pColNamesTemp, peptides))
		{
			cout <<"Failed to read in " << par.getFilePath(DBTemplate::colIndex) <<"!" << endl << "Exiting..." << endl;
			return false;
		}
		cout << "Adding " << par.getFileColname(DBTemplate::colIndex) << "..." << endl;
	}
	return true;
}

bool Proteins::readInMWdb(string wd, const filterFile::FilterFileParams& par)
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
	OUT_DELIM << pI <<
	OUT_DELIM << length <<
	OUT_DELIM << MW;
	
	if(par->calcMW)
		outF << OUT_DELIM << avgMass <<
		OUT_DELIM << monoMass;
	
	if(par->getFxn)
		outF << OUT_DELIM << fxn;
	
	if(par->getSeq)
		outF << OUT_DELIM << sequence;
	
	if(par->getSubCelluarLoc)
		outF << OUT_DELIM <<  loc;
	
	if(par->outputFormat == 1)
	{
		if(par->supInfoOutput == 0)
		{
			for(int i = 0; i < colSize; i++)
			{
				outF << OUT_DELIM << col[i].count;
				if (par->includeUnique)
					outF << OUT_DELIM << col[i].uniquePeptides;
				
				if(par->includeCoverage)
					outF << OUT_DELIM << col[i].coverage;
				
				if(par->includeSequenceCount)
					outF << OUT_DELIM << col[i].sequenceCount;
			}
		}
		else if(par->supInfoOutput == 1)
		{
			assert(par->supInfoNum > 0);
			if(par->supInfoNum == 1)
			{
				writeCount(outF);
				if(par->includeCoverage)
					writeCoverage(outF);
				else if(par->includeUnique)
					writeUnique(outF);
				else if(par->includeSequenceCount)
					writeSequenceCount(outF);
				else throw runtime_error("Bad pars!");
			}
			if(par->supInfoNum == 2)
			{
				writeCount(outF);
				if(par->includeUnique && par->includeCoverage)
				{
					writeUnique(outF);
					writeCoverage(outF);
				}
				else if(par->includeUnique && par->includeSequenceCount)
				{
					writeUnique(outF);
					writeSequenceCount(outF);
				}
				else if(par->includeCoverage && par->includeSequenceCount)
				{
					writeCoverage(outF);
					writeSequenceCount(outF);
				}
				else throw runtime_error("Bad pars!");
			}
			if(par->supInfoNum == 3)
			{
				writeCount(outF);
				writeUnique(outF);
				writeCoverage(outF);
				writeSequenceCount(outF);
			}
		}
	}
	else if(par->outputFormat == 2)
	{
		outF << OUT_DELIM << col[*colIndex].colname;
		
		if (!par->sampleNamePrefix.empty())
		{
			outF << OUT_DELIM << parseSample(col[*colIndex].colname, par->sampleNamePrefix, 1) << OUT_DELIM <<
			parseReplicate(col[*colIndex].colname);
		}
		
		outF << OUT_DELIM << col[*colIndex].count;
		
		if (par->includeUnique)
			outF << OUT_DELIM << col[*colIndex].uniquePeptides;
		
		if(par->includeCoverage)
			outF << OUT_DELIM << col[*colIndex].coverage;
		
		if(par->includeSequenceCount)
			outF << OUT_DELIM << col[*colIndex].sequenceCount;
	}
	outF << endl;
}

void Protein::writeCount(ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < colSize; i++)
		outF << OUT_DELIM << col[i].count;
}
void Protein::writeUnique(ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < colSize; i++)
		outF << OUT_DELIM << col[i].uniquePeptides;
}
void Protein::writeCoverage(ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < colSize; i++)
		outF << OUT_DELIM << col[i].coverage;
}
void Protein::writeSequenceCount(ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < colSize; i++)
		outF << OUT_DELIM << col[i].sequenceCount;
}

//write out combined protein lists to ofname in wide format
bool Proteins::writeOut(string ofname, const filterFile::FilterFileParams& par)
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	filterFile::OutputFormat outputFormat = par.outputFormat;
	par.outputFormat = filterFile::wideFormat;
	
	bool parseSampleName = par.sampleNamePrefix != "";
	int colNamesLength = DEFAULT_COL_NAMES_LENGTH;
	vector<string> headers;
	vector<string>::iterator it = headers.begin();
	headers.insert(it, DEFAULT_COL_NAMES, DEFAULT_COL_NAMES + colNamesLength);
	
	//print column headers
	if(par.getSubCelluarLoc)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Mass(Da)")
			{
				headers.insert(it + 1, "subcelluar_loc");
				break;
			}
	}
	if(par.getFxn)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Mass(Da)")
			{
				headers.insert(it + 1, "Function");
				break;
			}
	}
	if(par.calcMW)
	{
		int mwHeadersLen = MWCALC_HEADERS_LENGTH;
		if(par.getSeq)
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
	else if(par.getSeq && !par.calcMW)
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
		if(par.supInfoOutput == 0)
			 delim = utils::repeat(string(1, OUT_DELIM), par.supInfoNum + 1);
		else delim = string(1, OUT_DELIM);
		for (int i = 0; i < colNamesLength; i++)
			outF << OUT_DELIM;
		for (int i = 0; i < par.numFiles ; i++)
			outF << colNames[i] << delim;
		outF << endl;
	}
	if((par.includeUnique || par.includeCoverage || par.includeSequenceCount) && (par.supInfoOutput == 0))
	{
		string tabs = utils::repeat(string(1, OUT_DELIM), par.supInfoNum + 1);
		
		for (int i = 0; i < colNamesLength; i++)
			outF << OUT_DELIM;
		for (int i = 0; i < par.numFiles; i++)
		{
			if (i == 0)
				outF << parseSample(colNames[i], par.sampleNamePrefix, 0);
			else outF << tabs << parseSample(colNames[i], par.sampleNamePrefix, 0);
		}
		outF << endl;
		for (int i = 0; i < colNamesLength; i++)
			outF << headers[i] << OUT_DELIM;
		
		string repeatedHeaders = SUP_INFO_HEADERS[0];
		if(par.supInfoNum == 1)
		{
			if(par.includeUnique)
				repeatedHeaders += (OUT_DELIM + SUP_INFO_HEADERS[1]);
			else if(par.includeCoverage)
				repeatedHeaders += (OUT_DELIM + SUP_INFO_HEADERS[2]);
			else if(par.includeSequenceCount)
				repeatedHeaders += ( OUT_DELIM + SUP_INFO_HEADERS[3]);
			else return false;
		}
		else if(par.supInfoNum == 2)
		{
			if(par.includeUnique && par.includeCoverage)
				repeatedHeaders += (OUT_DELIM + SUP_INFO_HEADERS[1] + OUT_DELIM + SUP_INFO_HEADERS[2]);
			else if(par.includeUnique && par.includeSequenceCount)
				repeatedHeaders += (OUT_DELIM + SUP_INFO_HEADERS[1] + OUT_DELIM + SUP_INFO_HEADERS[3]);
			else if(par.includeCoverage && par.includeSequenceCount)
				repeatedHeaders += (OUT_DELIM + SUP_INFO_HEADERS[2] + OUT_DELIM + SUP_INFO_HEADERS[3]);
			else return false;
		}
		else if(par.supInfoNum == 3)
		{
			if(par.includeUnique && par.includeCoverage && par.includeSequenceCount)
				repeatedHeaders += (OUT_DELIM + SUP_INFO_HEADERS[1] + OUT_DELIM + SUP_INFO_HEADERS[2] + OUT_DELIM + SUP_INFO_HEADERS[3]);
			else return false;
		}
		
		for (int i = 0; i < par.numFiles; i++)
		{
			if(i == 0)
				outF << repeatedHeaders;
			else outF << OUT_DELIM << repeatedHeaders;
		}
		outF << endl;
	}
	else {
		int len = int(headers.size());
		if(par.supInfoOutput == 1)
		{
			assert(par.supInfoNum > 0);
			string preBuffer = utils::repeat(string(1, OUT_DELIM), len);
			string postBuffer = utils::repeat(string(1, OUT_DELIM), colNames.size());
			
			outF << preBuffer << "Spectral counts";
			
			if(par.supInfoNum == 1)
			{
				if(par.includeCoverage)
					outF << postBuffer << "Percent_coverage" << postBuffer;
				else if(par.includeUnique)
					outF << postBuffer << "Unique_peptide_SC" << postBuffer;
				else if(par.includeSequenceCount)
					outF << preBuffer << "Sequence_count" << postBuffer;
				else return false;
			}
			else if(par.supInfoNum == 2)
			{
				if(par.includeCoverage && par.includeUnique)
					outF << postBuffer << "Unique_peptide_SC" <<
					postBuffer << "Percent_coverage"; // << postBuffer;
				else if(par.includeCoverage && par.includeSequenceCount)
					outF << postBuffer << "Percent_coverage" <<
					postBuffer << "Sequence_count"; // << postBuffer;
				else if(par.includeSequenceCount && par.includeCoverage)
					outF << postBuffer << "Unique_peptide_SC" <<
					postBuffer << "Sequence_count"; // << postBuffer;
				else return false;
			}
			else if(par.supInfoNum == 3)
			{
				outF << postBuffer << "Unique_peptide_SC" <<
				postBuffer << "Percent_coverage" <<
				postBuffer << "Sequence_count";
				
			}
			else return false;
			outF << endl;
		}

		vector<string> ofColNames;
		for (int i = 0; i < len; i ++)
			ofColNames.push_back(headers[i]);
		for(int i = 0; i <= par.supInfoNum; i++)
		{
			for (int i = 0; i < par.numFiles; i ++)
				ofColNames.push_back(parseSample(colNames[i], par.sampleNamePrefix, 0));
		}
		int colNamesLen = int(ofColNames.size());
		for (int i = 0; i < colNamesLen; i++)
		{
			if(i == 0)
				outF << ofColNames[i];
			else outF << OUT_DELIM << ofColNames[i];
		}
		outF << endl;
	}
	
	//print proteins and spectral counts
	data->write(outF);
	
	par.outputFormat = outputFormat;
	
	return true;
}

//write out combined protein lists to ofname in long format
bool Proteins::writeOutDB(string ofname, const filterFile::FilterFileParams& par)
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	filterFile::OutputFormat outputFormat = par.outputFormat;
	par.outputFormat = filterFile::longFormat;

	bool parseSampleName = !(par.sampleNamePrefix.empty());
	
	//print column headers
	vector<string> headers;
	vector<string>::iterator it = headers.begin();
	int len = DEFAULT_COL_NAMES_DB_LENGTH;
	
	for (int i = 0; i < DEFAULT_COL_NAMES_DB_LENGTH ; i++)
		headers.push_back(DEFAULT_COL_NAMES_DB[i]);
	
	if(parseSampleName)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Long_sample_name")
			{
				headers.insert(it + 1, PARSE_SAMPLE_NAME_HEADERS, PARSE_SAMPLE_NAME_HEADERS + PARSE_SAMPLE_NAME_HEADERS_LEN);
				len = int(headers.size());
				break;
			}
	}
	if(par.getSubCelluarLoc)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Mass(Da)")
			{
				headers.insert(it + 1, "subcelluar_loc");
				len = int(headers.size());
				break;
			}
	}
	if(par.getFxn)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Mass(Da)")
			{
				headers.insert(it + 1, "Function");
				len = int(headers.size());
				break;
			}
	}
	if(par.calcMW)
	{
		int mwHeadersLen = MWCALC_HEADERS_LENGTH;
		if(par.getSeq)
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
	else if(par.getSeq && !par.calcMW)
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
	if(par.includeSequenceCount)
	{
		for(it = headers.begin(); it != headers.end(); it++)
		{
			if(*it == "Spectral_counts")
			{
				headers.insert(it + 1, "sequence_count");
				len = int(headers.size());
				break;
			}
		}
	}
	if(par.includeCoverage)
	{
		for(it = headers.begin(); it != headers.end(); it++)
		{
			if(*it == "Spectral_counts")
			{
				headers.insert(it + 1, "percent_coverage");
				len = int(headers.size());
				break;
			}
		}
	}
	if(par.includeUnique)
	{
		for(it = headers.begin(); it != headers.end(); it++)
		{
			if(*it == "Spectral_counts")
			{
				headers.insert(it + 1, SUP_INFO_HEADERS[1]);
				len = int(headers.size());
				break;
			}
		}
	}

	for (it = headers.begin(); it != headers.end(); it++)
	{
		if(it == headers.begin())
			outF << *it;
		else outF << OUT_DELIM << *it;
	}
	outF << endl;
	
	Protein::colIndex = &colIndex;
	//print proteins and spectral counts
	for (colIndex = 0; colIndex < par.numFiles; colIndex++)
		data->write(outF);
	
	par.outputFormat = outputFormat;
	
	return true;
}

//check if line containing % is a collumn header line instead of a protein header line
inline bool isColumnHeaderLine(const string& line)
{
	return utils::strContains("Conf%", line);
}

//optional fxn to parse long sample name
string parseSample(string sampleName, string prefix, bool outputFormat)
{
	//return unparsed sampleName if prefix is empty string or is not found in sampleName
	if(!utils::strContains(prefix, sampleName) || prefix.length() == 0)
		return sampleName;
	
	string sample = utils::removeSubstr(prefix, sampleName); //remove prefix from sampleName
	
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
	utils::split(line, IN_DELIM, elems);
	assert(elems.size() > 11);
	
	//return SC for peptide as int
	return utils::toInt(elems[11]);
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

	return (firstP == string::npos || secP == string::npos) ?
		str : str.substr(firstP + 1, secP - (str.length() - secP));
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
	
	if(par->outputFormat == 2)
		if(!par->includeNullPeptides && col[*colIndex].isNull())
			return;
	
	outF << proteinID << OUT_DELIM <<
	protein << OUT_DELIM <<
	description << OUT_DELIM <<
	sequence << OUT_DELIM <<
	length << OUT_DELIM;
	
	if(par->peptideGroupMethod != filterFile::byCharge)
		outF << charge << OUT_DELIM;
	
	outF << unique << OUT_DELIM <<
	calcMH << OUT_DELIM;
	
	if(par->calcMW)
		outF << avgMass << OUT_DELIM <<
		monoMass << OUT_DELIM;
	
	if(par->outputFormat == 1)
	{
		for(int i = 0; i < colSize; i++)
		{
			if(i == 0)
				outF << col[i].count;
			else outF << OUT_DELIM << col[i].count;
		}
	}
	else if(par->outputFormat == 2)
	{
		if(par->peptideGroupMethod != filterFile::byCharge)
			outF << col[*colIndex].obsMH << OUT_DELIM <<
			col[*colIndex].scan << OUT_DELIM <<
			col[*colIndex].parentFile << OUT_DELIM;
		
		outF << col[*colIndex].colname << OUT_DELIM <<
		col[*colIndex].count;
		
		if(!par->sampleNamePrefix.empty())
		{
			outF << OUT_DELIM << parseSample(col[*colIndex].colname, par->sampleNamePrefix, 1)
			<< OUT_DELIM <<	parseReplicate(col[*colIndex].colname);
		}
	}
	
	outF << endl;
}

bool Peptides::writeOut(string ofname, const filterFile::FilterFileParams& pars)
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	filterFile::OutputFormat outputFormat = pars.outputFormat;
	pars.outputFormat = filterFile::wideFormat;
	
	//generate headers based off params
	vector<string> headers;
	vector<string>::iterator it;
	for(int i = 0; i < DEFALUT_PEPTIDE_COLNAMES_LEN; i++)
		headers.push_back(DEFALUT_PEPTIDE_COLNAMES[i]);
	if(pars.peptideGroupMethod != filterFile::byCharge)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "length(aa)")
			{
				headers.insert(it + 1, "charge");
				break;
			}
	}
	if(pars.calcMW)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "calcMH")
			{
				headers.insert(it + 1, MWCALC_HEADERS, MWCALC_HEADERS + 2);
				break;
			}
	}
	headers.insert(headers.end(), colNames.begin(), colNames.end());
	
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


bool Peptides::writeOutDB(string ofname, const filterFile::FilterFileParams& pars)
{
	ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	filterFile::OutputFormat outputFormat = pars.outputFormat;
	pars.outputFormat = filterFile::longFormat;
	
	//generate headers based off params
	vector<string> headers;
	vector<string>::iterator it;
	int defaultColNamesLen = DEFALUT_PEPTIDE_DB_COLNAMES_LEN;
	
	if(!pars.sampleNamePrefix.empty())
		defaultColNamesLen += 2;
	
	for(int i = 0; i < defaultColNamesLen; i++)
		headers.push_back(DEFALUT_PEPTIDE_DB_COLNAMES[i]);
	if(pars.peptideGroupMethod != filterFile::byCharge)
	{
		//headers to add
		string add [] = {"obsMH", "Scan", "parent_file"};
		for(it = headers.begin(); it != headers.end(); it++)
		{
			if(*it == "length(aa)")
			{
				headers.insert(it + 1, "charge");
				break;
			}
		}
		for(it = headers.begin(); it != headers.end(); it++)
		{
			if(*it == "calcMH")
			{
				headers.insert(it + 1, add, add + 3);
				break;
			}
		}
	}
	if(pars.calcMW)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "unique")
			{
				headers.insert(it + 1, MWCALC_HEADERS, MWCALC_HEADERS + 2);
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
	
	Peptide::colIndex = &colIndex;
	for(colIndex = 0; colIndex < pars.numFiles; colIndex++)
		data->write(outF);
	
	pars.outputFormat = outputFormat;
	
	return true;
}
