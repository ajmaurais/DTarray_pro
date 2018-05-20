//
//  dtafilter.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 3/25/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include <dtafilter.hpp>

void Protein::consolidate(const Protein& toAdd)
{
	col[*toAdd.colIndex].colname = toAdd.col[*toAdd.colIndex].colname;
	col[*toAdd.colIndex].count = toAdd.col[*toAdd.colIndex].count;
	col[*toAdd.colIndex].uniquePeptides = toAdd.col[*toAdd.colIndex].uniquePeptides;
	col[*toAdd.colIndex].coverage = toAdd.col[*toAdd.colIndex].coverage;
	col[*toAdd.colIndex].sequenceCount = toAdd.col[*toAdd.colIndex].sequenceCount;
	col[*toAdd.colIndex].modPeptides = toAdd.col[*toAdd.colIndex].modPeptides;
	col[*toAdd.colIndex].modPeptidesSC = toAdd.col[*toAdd.colIndex].modPeptidesSC;
}

void Peptide::consolidate(const Peptide& toAdd)
{
	col[*toAdd.colIndex].colname = toAdd.col[*toAdd.colIndex].colname;
	col[*toAdd.colIndex].count += toAdd.col[*toAdd.colIndex].count;
	col[*toAdd.colIndex].scan = toAdd.col[*toAdd.colIndex].scan;
	col[*toAdd.colIndex].parentFile = toAdd.col[*toAdd.colIndex].parentFile;
	col[*toAdd.colIndex].obsMH = toAdd.col[*toAdd.colIndex].obsMH;
	col[*toAdd.colIndex].modPeptidesSC += toAdd.col[*toAdd.colIndex].modPeptidesSC;
	col[*toAdd.colIndex].xCorr = toAdd.col[*toAdd.colIndex].xCorr;
}

template <class _Tp>
void ProteinDataTemplate<_Tp>::initialize(const std::vector<_Tp>& tempColNames, size_t colNamesLen, size_t* _colIndex)
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

void Protein::operator = (const Protein& _p)
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

void Peptide::operator = (const Peptide& _p)
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
void Protein::getProteinData(std::string line)
{
	//split line by tabs
	std::vector<std::string> elems;
	utils::split(line, IN_DELIM, elems);
	assert(elems.size() == 9);
	
	//keep fullDescription but seperate by spaces instead of tabs
	fullDescription = elems[0];
	for(int i = 1; i < 9; i++)
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
	if(utils::strContains(REVERSE_MATCH, matchDirrection))
		ID = "reverse_" + ID;
	
	//extract MW
	MW = elems[5];
	
	//extract shortened protein name and description
	size_t endOfDescription = elems[8].find(" [");
	description = elems[8].substr(0, endOfDescription);
	getProtein(description);
	
	//add spectrum count, coverage and sequence count for *this protein to colname
	col[*colIndex].count = utils::toInt(elems[2]);
	std::string coverageTemp = elems[3];
	coverageTemp = coverageTemp.substr(0, coverageTemp.find("%"));
	col[*colIndex].coverage = coverageTemp;
	col[*colIndex].sequenceCount = elems[1];
}

void Protein::getProtein(std::string str)
{
	size_t firstSpace = str.find(" ");
	
	if(firstSpace == std::string::npos)
		protein = str;
	else{
		protein = str.substr(0, firstSpace);
	}
}

void Protein::calcMW()
{
	std::string tmp_sequence = mwdb->seqDB->getSequence(ID);
	
	if(tmp_sequence == SEQ_NOT_FOUND)
	{
		sequence = SEQ_NOT_FOUND;
		avgMass = -1;
		monoMass = -1;
		formula = "NA";
		return;
	}
	
	sequence = tmp_sequence;
	avgMass = mwdb->calcMW(sequence);
	monoMass = mwdb->calcMono(sequence);
	formula = mwdb->calcFormula(sequence, par->unicode);
}

void Peptide::calcMW()
{
	avgMass = mwdb->calcMW(calcSequence);
	monoMass = mwdb->calcMono(calcSequence);
	formula = mwdb->calcFormula(calcSequence, par->unicode);
}

//loop through protein header and peptide lines in DTA filter file and add data to Proteins
bool Proteins::readIn(params::Params* const pars,
					  const std::vector <SampleData_protein>& colNamesTemp,
					  const std::vector<SampleData_peptide>& pColNamesTemp,
					  Peptides * const peptides)
{
	std::string fname = pars->getwd() + pars->getFilePath(colIndex);
	utils::File file;
	if(!file.read(fname))
		return false;
	
	DataType::iterator proteinIndex;
	bool inProtein = false;
	bool getNewLine = true;
	bool foundHeader = false;
	std::string line;
	std::map<std::string, Peptide>::iterator peptidesIndex;
	size_t colNamesLen = pars->getNumFiles();
	
	while(!file.end()){
		if(getNewLine)
			line = file.getLine();
		getNewLine = true;
		if(utils::strContains('%', line)) //find protein header lines by percent symbol for percent coverage
		{
			if(!foundHeader)
			{
				if(utils::strContains("Conf%", line)) //skip if header line
				{
					foundHeader = true;
					continue;
				}
			}
			else{
				//initalize Protein to hold data for current line
				Protein newProtein(pars, locDB, fxnDB, mwdb, seqDB, baitFile, locTable);
				newProtein.initialize(colNamesTemp, colNamesLen, &colIndex);
				newProtein.getProteinData(line);
				inProtein = true; //used to determine if it is valid to loop through peptide lines below protein
								  //header line to extract unique peptide spectral counts
					
				//proteinIndex = data->consolidate(newProtein, newProtein.ID); //insert new protein into proteins list
				proteinIndex = data.find(newProtein.getID());
				if(proteinIndex == data.end()){
					data[newProtein.getID()] = newProtein;
				}
				else{
					data[newProtein.getID()].consolidate(newProtein);
				}
				proteinIndex = data.find(newProtein.getID());
				
				//get peptide data if applicable
				if((pars->includeUnique || pars->includePeptides || pars->includeModStat) && inProtein)
				{
					line = file.getLine();
					getNewLine = false;
					while(!file.end())
					{
						if(getNewLine)
							line = file.getLine();
						getNewLine = true;
						
						//break if starting new protein or end of file
						if(utils::strContains('%', line) || line == "\tProteins\tPeptide IDs\tSpectra")
							break;
						
						if(pars->includePeptides)
						{
							Peptide newPeptide(pars, mwdb);
							newPeptide.initialize(pColNamesTemp, colNamesLen, &peptides->colIndex);
							newPeptide.proteinID = newProtein.ID;
							newPeptide.protein = newProtein.protein;
							newPeptide.description = newProtein.description;
							newPeptide.matchDirrection = newProtein.matchDirrection;
							newPeptide.parsePeptide(line);
							
							if(pars->peptideGroupMethod == params::byScan){
								peptides->data[newPeptide.key] = newPeptide;
							}
							else{
								peptidesIndex = peptides->data.find(newPeptide.key);
								if(peptidesIndex == peptides->data.end()){
									peptides->data[newPeptide.key] = newPeptide;
								}
								else{
									peptides->data[newPeptide.key].consolidate(newPeptide);
								}
							}
						}
						if(pars->includeUnique)
						{
							if(line[0] == '*')
								proteinIndex->second.col[colIndex].uniquePeptides += parsePeptideSC(line);
						}
						if(pars->includeModStat)
						{
							int modPeptideSC = parseModPeptide(line);
							if(modPeptideSC > 0)
							{
								proteinIndex->second.col[colIndex].modPeptidesSC += modPeptideSC;
								proteinIndex->second.col[colIndex].modPeptides++;
							}
							
						}
					}//end of while
					getNewLine = false;
					inProtein = false;
				}//end of if
			}//end of else
		}//end of if
	}//end of while
	return true;
}

void Peptide::parsePeptide(const std::string& line)
{
	std::vector<std::string> elems;
	utils::split(line, IN_DELIM, elems);
	
	assert(elems.size() == 13);
	
	unique = elems[0] == "*";
	parseSequence(elems[12]); //get calcSequence
	length = utils::toString(calcSequence.length());
	sequence = elems[12];
	col[*colIndex].obsMH = elems[5];
	col[*colIndex].count = utils::toInt(elems[11]);
	col[*colIndex].xCorr = elems[2];
	
	if(par->includeModStat)
		col[*colIndex].modPeptidesSC += parseModPeptide(line);
	
	calcMH = elems[6];
	fileName = elems[1];
	
	std::string temp = elems[1];
	elems.clear();
	utils::split(temp, '.', elems);
	charge = elems[3];
	key = makeKey();
	col[*colIndex].scan = elems[2];
	col[*colIndex].parentFile = elems[0];
}

std::string Peptide::makeKey() const {
	switch(par->peptideGroupMethod){
		case params::byScan : return fileName;
			break;
		case params::byProtein : return proteinID + "_" + calcSequence + "_" + charge;
			break;
		case params::byCharge : return proteinID + "_" + calcSequence;
			break;
		default : throw std::runtime_error("Invalid type");
	}
}

void Protein::clear()
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

void Peptide::clear()
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
bool Proteins::readIn(params::Params* const par, Peptides* const peptides)
{
	size_t colNamesLen = par->getNumFiles();
	
	std::vector<SampleData_protein> colNamesTemp;
	std::vector<SampleData_peptide> pColNamesTemp;
	
	for(int i = 0; i < colNamesLen; i++)
	{
		SampleData_protein temp (colNames[i]);
		colNamesTemp.push_back(temp);
		SampleData_peptide pTemp (colNames[i]);
		pColNamesTemp.push_back(pTemp);
	}
	
	for(Proteins::colIndex = 0; Proteins::colIndex < par->getNumFiles(); Proteins::colIndex++)
	{
		Peptides::colIndex = Proteins::colIndex;
		if(!readIn(par, colNamesTemp, pColNamesTemp, peptides))
		{
			std::cerr <<"Failed to read in " << par->getFilePath(Proteins::colIndex) <<"!" << std::endl << "Exiting..." << std::endl;
			return false;
		}
		std::cerr << "Adding " << par->getFileColname(Proteins::colIndex) << "..." << std::endl;
	}
	return true;
}

bool Proteins::readInMWdb(const params::Params& par)
{
	mwdb = new mwDB::MWDB_Protein;
	return mwdb->initalize(par);
}

bool Proteins::readInSeqDB(std::string fname)
{
	seqDB = new mwDB::SeqDB;
	return seqDB->readIn(fname);
}

bool Proteins::readInFxnDB(std::string fname)
{
	fxnDB = new Dbase;
	return fxnDB->readIn(fname);
}

bool Proteins::readInLocDB(std::string fname)
{
	locDB = new Dbase;
	return locDB->readIn(fname);
}

bool Proteins::readBaitFile(std::string fname)
{
	baitFile = new saint::BaitFile(fname);
	return baitFile->read();
}

locReport::LocDat Protein::toLocDat(const std::string& _loc) const
{
	locReport::LocDat newLocDat(_loc, matchDirrection, colSize, par);
	
	for(size_t i = 0; i < colSize; i++)
		newLocDat.initializeCol(i, col[i].colname, (col[i].count > 0), col[i].count, col[i].uniquePeptides, utils::toInt(col[i].sequenceCount));
	
	return newLocDat;
}

void Protein::addLocToTable()
{
	if(!supDataAdded)
		addSupData();
	
	std::vector<std::string> elems;
	utils::split(loc, DB_DELIM, elems);
	
	for(std::vector<std::string>::iterator it = elems.begin(); it != elems.end(); ++it)
		locTable->addLoc(toLocDat(*it));
}

/*void Protein::apply(int fxnNum)
{
	switch(fxnNum){
		case 0 : addLocToTable();
			break;
		default : throw std::runtime_error("Function does not exist!");
	}
}*/

void Protein::writePrey(std::ofstream& outF) const
{
	assert(outF);
	outF << ID << OUT_DELIM <<
	length << OUT_DELIM << description << std::endl;
}

void Protein::writeInteractions(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < colSize; i++)
	{
		if(col[i].count > 0)
		{
			outF << col[i].colname << OUT_DELIM <<
			baitFile->getBaitName(col[i].colname) << OUT_DELIM <<
			ID << OUT_DELIM <<
			col[i].count << std::endl;
		}
	}
}

void Protein::addSupData()
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

void Protein::writeProtein(std::ofstream& outF)
{
	if(!outF)
		throw std::runtime_error("Bad std::ofstream");
	
	if(!par->includeReverse)
		if(utils::strContains(REVERSE_MATCH, matchDirrection))
			return;
	
	if(!supDataAdded)
		addSupData();
	
	if(INCLUDE_FULL_DESCRIPTION)
		outF << fullDescription << OUT_DELIM;
	
	outF << ID <<
	OUT_DELIM << protein <<
	OUT_DELIM << description <<
	OUT_DELIM << pI <<
	OUT_DELIM << length <<
	OUT_DELIM << MW;
	
	std::streamsize ss = outF.precision();
	outF.precision(5); //write floating point nums with 5 digits
	if(par->calcMW)
		outF << OUT_DELIM << std::fixed << avgMass <<
		OUT_DELIM << std::fixed << monoMass <<
		OUT_DELIM << formula;
	outF.precision(ss);
	
	if(par->getSeq)
		outF << OUT_DELIM << sequence;
	
	if(par->getFxn)
		outF << OUT_DELIM << fxn;
	
	if(par->getSubCelluarLoc)
		outF << OUT_DELIM << loc;
	
	if(par->outputFormat == 1)
	{
		assert(par->supInfoNum >= 0 && par->supInfoNum <= 5);
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
				
				if(par->includeModStat)
					outF << OUT_DELIM << col[i].modPeptides
						<< OUT_DELIM << col[i].modPeptidesSC;
			}
		}
		else if(par->supInfoOutput == 1)
		{
			writeCount(outF);
			
			if(par->includeUnique)
				writeUnique(outF);
			
			if(par->includeCoverage)
				writeCoverage(outF);
			
			if(par->includeSequenceCount)
				writeSequenceCount(outF);
			
			if(par->includeModStat)
				writeModStat(outF);
		}
	}
	else if(par->outputFormat == 2)
	{
		outF << OUT_DELIM << col[*colIndex].colname;
		
		if (par->parseSampleName)
		{
			outF << OUT_DELIM << parseSample(col[*colIndex].colname, par->sampleNamePrefix, par->parseSampleName, 1)
			<< OUT_DELIM << parseReplicate(col[*colIndex].colname);
		}
		
		outF << OUT_DELIM << col[*colIndex].count;
		
		if (par->includeUnique)
			outF << OUT_DELIM << col[*colIndex].uniquePeptides;
		
		if(par->includeCoverage)
			outF << OUT_DELIM << col[*colIndex].coverage;
		
		if(par->includeSequenceCount)
			outF << OUT_DELIM << col[*colIndex].sequenceCount;
		
		if(par->includeModStat)
			outF << OUT_DELIM << col[*colIndex].modPeptides
			<< OUT_DELIM << col[*colIndex].modPeptidesSC;
	}
	outF << std::endl;
}

void Protein::writeCount(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < colSize; i++)
		outF << OUT_DELIM << col[i].count;
}
void Protein::writeUnique(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < colSize; i++)
		outF << OUT_DELIM << col[i].uniquePeptides;
}
void Protein::writeCoverage(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < colSize; i++)
		outF << OUT_DELIM << col[i].coverage;
}
void Protein::writeSequenceCount(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < colSize; i++)
		outF << OUT_DELIM << col[i].sequenceCount;
}

void Protein::writeModStat(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < colSize; i++)
		outF << OUT_DELIM << col[i].modPeptides
			<< OUT_DELIM << col[i].modPeptidesSC;
}

//write out combined protein lists to ofname in wide format
bool Proteins::writeOut(std::string ofname, const params::Params& par)
{
	std::ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	params::OutputFormat outputFormat = par.outputFormat;
	par.outputFormat = params::wideFormat;
	
	int colNamesLength = DEFAULT_COL_NAMES_LENGTH;
	
	bool supInfoS [] = {true, par.includeUnique, par.includeCoverage, par.includeSequenceCount,
		par.includeModStat, par.includeModStat};
	bool supInfo = false;
	std::vector<std::string> supInfoHeaders;
	for(int i = 0; i < SUP_INFO_HEADERS_LEN; i++)
	{
		if(supInfoS[i])
		{
			supInfoHeaders.push_back(SUP_INFO_HEADERS[i]);
			if(i > 0)
				supInfo = true;
		}
	}
	
	std::vector<std::string> headers;
	std::vector<std::string>::iterator it = headers.begin();
	headers.insert(it, DEFAULT_COL_NAMES, DEFAULT_COL_NAMES + colNamesLength);
	
	//print column headers
	if(par.getSubCelluarLoc)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Mass(Da)")
			{
				headers.insert(it + 1, "subcelluar_loc");
				colNamesLength = int(headers.size());
				break;
			}
	}
	if(par.getFxn)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Mass(Da)")
			{
				headers.insert(it + 1, "Function");
				colNamesLength = int(headers.size());
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
				colNamesLength = int(headers.size());
				break;
			}
	}
	if(!par.sampleNamePrefix.empty())
	{
		std::string delim;
		if(par.supInfoOutput == 0)
			 delim = utils::repeat(std::string(1, OUT_DELIM), par.supInfoNum + 1);
		else delim = std::string(1, OUT_DELIM);
		for(int i = 0; i < colNamesLength; i++)
			outF << OUT_DELIM;
		for(int i = 0; i < par.getNumFiles(); i++)
			outF << colNames[i] << delim;
		outF << std::endl;
	}
	if(supInfo && (par.supInfoOutput == 0))
	{
		std::string tabs = utils::repeat(std::string(1, OUT_DELIM), par.supInfoNum + 1);
		
		std::string repeatHeaders;
		for(std::vector<std::string>::iterator it = supInfoHeaders.begin(); it != supInfoHeaders.end(); ++it)
		{
			if(it == supInfoHeaders.begin())
				repeatHeaders = *it;
			else repeatHeaders += (OUT_DELIM + *it);
		}
		
		for(int i = 0; i < colNamesLength; i++)
			outF << OUT_DELIM;
		for(int i = 0; i < par.getNumFiles(); i++)
		{
			if (i == 0)
				outF << parseSample(colNames[i], par.sampleNamePrefix, false, 0);
			else outF << tabs << parseSample(colNames[i], par.sampleNamePrefix, false, 0);
		}
		outF << std::endl;
		for(int i = 0; i < colNamesLength; i++)
			outF << headers[i] << OUT_DELIM;
		
		for (int i = 0; i < par.getNumFiles(); i++)
		{
			if(i == 0)
				outF << repeatHeaders;
			else outF << OUT_DELIM << repeatHeaders;
		}
		outF << std::endl;
	}
	else {
		int len = int(headers.size());
		if(par.supInfoOutput == 1 && supInfo)
		{
			assert(par.supInfoNum >= 1 && par.supInfoNum <= 5);
			std::string preBuffer = utils::repeat(std::string(1, OUT_DELIM), len);
			std::string postBuffer = utils::repeat(std::string(1, OUT_DELIM), colNames.size());
			
			outF << preBuffer;
			for(std::vector<std::string>::iterator it = supInfoHeaders.begin(); it != supInfoHeaders.end(); ++it)
				outF << *it << postBuffer;
			outF << std::endl;
		}

		std::vector<std::string> ofColNames;
		for(int i = 0; i < len; i++)
			ofColNames.push_back(headers[i]);
		for(int i = 0; i <= par.supInfoNum; i++)
		{
			for(int i = 0; i < par.getNumFiles(); i ++)
				ofColNames.push_back(parseSample(colNames[i], par.sampleNamePrefix, false, 0));
		}
		int colNamesLen = int(ofColNames.size());
		for(int i = 0; i < colNamesLen; i++)
		{
			if(i == 0)
				outF << ofColNames[i];
			else outF << OUT_DELIM << ofColNames[i];
		}
		outF << std::endl;
	}
	
	//print proteins and spectral counts
	//data->write(outF);
	for(DataType::iterator it = data.begin(); it != data.end(); ++it){
		it->second.writeProtein(outF);
	}
	
	par.outputFormat = outputFormat;
	
	return true;
}

//write out combined protein lists to ofname in long format
bool Proteins::writeOutDB(std::string ofname, const params::Params& par)
{
	std::ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	params::OutputFormat outputFormat = par.outputFormat;
	par.outputFormat = params::longFormat;
	
	//print column headers
	std::vector<std::string> headers;
	std::vector<std::string>::iterator it = headers.begin();
	int len = DEFAULT_COL_NAMES_DB_LENGTH;
	
	for (int i = 0; i < DEFAULT_COL_NAMES_DB_LENGTH ; i++)
		headers.push_back(DEFAULT_COL_NAMES_DB[i]);
	
	if(par.parseSampleName)
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
				headers.insert(it + 1, "Subcelluar_loc");
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
	if(par.includeModStat)
	{
		for(it = headers.begin(); it != headers.end(); it++)
		{
			if(*it == "Spectral_counts")
			{
				headers.insert(it + 1, SUP_INFO_HEADERS[5]);
				headers.insert(it + 1, SUP_INFO_HEADERS[4]);
				len = int(headers.size());
				break;
			}
		}
	}
	if(par.includeSequenceCount)
	{
		for(it = headers.begin(); it != headers.end(); it++)
		{
			if(*it == "Spectral_counts")
			{
				headers.insert(it + 1, SUP_INFO_HEADERS[3]);
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
				headers.insert(it + 1, SUP_INFO_HEADERS[2]);
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
	outF << std::endl;
	
	Protein::colIndex = &colIndex;
	//print proteins and spectral counts
	for(DataType::iterator it = data.begin(); it != data.end(); ++it){
		for(colIndex = 0; colIndex < par.getNumFiles(); colIndex++){
			it->second.writeProtein(outF);
		}
	}
	
	par.outputFormat = outputFormat;
	
	return true;
}

bool Proteins::writeSaint(std::string fname, OutputFiles file) const
{
	std::ofstream outF(fname.c_str());
	
	if(!outF)
		return false;
	
	if(file == preyFile){
		for(DataType::const_iterator it = data.begin(); it != data.end(); ++it)
			it->second.writePrey(outF);
	}
	else if(file == interactionFile){
		for(DataType::const_iterator it = data.begin(); it != data.end(); ++it)
			it->second.writeInteractions(outF);
	}
	
	return true;
}

void Proteins::buildLocTable()
{
	//apply loc build fxn across data hash table
	locTable = new locReport::LocDB;
	Protein::locTable = locTable;
	//data->apply(0);
	for(DataType::iterator it = data.begin(); it != data.end(); ++it)
		it->second.addLocToTable();
}

bool Proteins::writeLongLocTable(std::string fname, const params::Params& pars) const
{
	std::ofstream outF(fname.c_str());
	
	if(!outF)
		return false;
	
	std::vector<std::string> headers;
	std::vector<std::string>::iterator it;
	headers.push_back("Location");
	
	if(pars.parseSampleName)
	{
		headers.push_back("Long_sample_name");
		headers.push_back("Sample");
		headers.push_back("Replicate");
	}
	else{
		headers.push_back("Sample");
	}
	
	headers.insert(headers.end(), utils::begin(LOC_REPORT_HEADERS), utils::end(LOC_REPORT_HEADERS));
	
	if(pars.includeUnique)
	{
		for(it = headers.begin(); it != headers.end(); ++it)
		{
			if(*it == "Sum_SC")
			{
				headers.insert(it + 1, "Sum_uniq_SC");
				break;
			}
		}
	}
	
	for(it = headers.begin(); it != headers.end(); ++it)
	{
		if(it == headers.begin())
			outF << *it;
		else outF << OUT_DELIM << *it;
	}
	outF << std::endl;
	
	locTable->writeLocReport(outF, 1);
	
	return true;
}

bool Proteins::writeWideLocTable(std::string fname, const params::Params& pars) const
{
	std::ofstream outF(fname.c_str());
	
	if(!outF)
		return false;
	
	std::vector<std::string> headers;
	std::vector<std::string> repeatedHeadersV;
	std::vector<std::string>::iterator it;
	headers.push_back("Location");
	
	int supInfoNum = pars.includeUnique + 3;
	
	repeatedHeadersV.insert(repeatedHeadersV.end(), utils::begin(LOC_REPORT_HEADERS), utils::end(LOC_REPORT_HEADERS));
	
	if(pars.includeUnique)
	{
		for(it = repeatedHeadersV.begin(); it != repeatedHeadersV.end(); ++it)
		{
			if(*it == "Sum_SC")
			{
				repeatedHeadersV.insert(it + 1, "Sum_uniq_SC");
				break;
			}
		}
	}
	
	if(!pars.sampleNamePrefix.empty())
	{
		std::string delim;
		if(pars.supInfoOutput == 0)
			delim = utils::repeat(std::string(1, OUT_DELIM), supInfoNum);
		else delim = std::string(1, OUT_DELIM);
		for(int i = 0; i < headers.size(); i++)
			outF << OUT_DELIM;
		for(int i = 0; i < pars.getNumFiles() ; i++)
			outF << colNames[i] << delim;
		outF << std::endl;
	}
	
	assert(pars.locSupInfoNum >= 0 && pars.locSupInfoNum <=1);
	if(pars.supInfoOutput == 0)
	{
		std::string tabs = utils::repeat(std::string(1, OUT_DELIM), supInfoNum);
		
		std::string repeatHeaders;
		for(it = repeatedHeadersV.begin(); it != repeatedHeadersV.end(); ++it)
		{
			if(it == repeatedHeadersV.begin())
				repeatHeaders = *it;
			else repeatHeaders += (OUT_DELIM + *it);
		}
		
		for(int i = 0; i < headers.size(); i++)
			outF << OUT_DELIM;
		for(int i = 0; i < pars.getNumFiles(); i++)
		{
			if (i == 0)
			outF << parseSample(colNames[i], pars.sampleNamePrefix, false, 0);
			else outF << tabs << parseSample(colNames[i], pars.sampleNamePrefix, false, 0);
		}
		outF << std::endl;
		for(int i = 0; i < headers.size(); i++)
		outF << headers[i] << OUT_DELIM;
		
		for (int i = 0; i < pars.getNumFiles(); i++)
		{
			if(i == 0)
			outF << repeatHeaders;
			else outF << OUT_DELIM << repeatHeaders;
		}
		outF << std::endl;
	}
	else if(pars.supInfoOutput == 1)
	{
		std::string preBuffer = utils::repeat(std::string(1, OUT_DELIM), headers.size());
		std::string postBuffer = utils::repeat(std::string(1, OUT_DELIM), colNames.size());
		
		outF << preBuffer;
		for(std::vector<std::string>::iterator it = repeatedHeadersV.begin(); it != repeatedHeadersV.end(); ++it)
			outF << *it << postBuffer;
		outF << std::endl;
		
		std::vector<std::string> ofColNames;
		for(int i = 0; i < headers.size(); i ++)
			ofColNames.push_back(headers[i]);
		for(int i = 0; i <= supInfoNum - 1; i++)
		{
			for(int i = 0; i < pars.getNumFiles(); i ++)
				ofColNames.push_back(parseSample(colNames[i], pars.sampleNamePrefix, false, 0));
		}
		int colNamesLen = int(ofColNames.size());
		for(int i = 0; i < colNamesLen; i++)
		{
			if(i == 0)
				outF << ofColNames[i];
			else outF << OUT_DELIM << ofColNames[i];
		}
		outF << std::endl;
	}
	
	locTable->writeLocReport(outF, 0);
	
	return true;
}

//optional fxn to parse long sample name
std::string parseSample(std::string sampleName, std::string prefix, bool parseSampleName, bool outputFormat)
{
	//return unparsed sampleName if prefix is empty std::string or is not found in sampleName
	if(!parseSampleName && (!utils::strContains(prefix, sampleName) || prefix.length() == 0)){
		return sampleName;
	}
	else if(parseSampleName && prefix.empty()){
		return sampleName.substr(0, sampleName.find_last_of("_"));
	}
	else {
		std::string sample = utils::removeSubstr(prefix, sampleName); //remove prefix from sampleName
		return outputFormat ? sample.substr(0, sample.find_last_of("_")) : sample;
	}
}

//get replicate number from sample name
std::string parseReplicate(std::string str)
{
	return str.substr(str.find_last_of("_")+1);
}

//return number of spectral counts for a peptide from a line in DTA filter file
int parsePeptideSC(std::string line)
{
	//split line by tabs
	std::vector<std::string> elems;
	utils::split(line, IN_DELIM, elems);
	assert(elems.size() == 13);
	
	//return SC for peptide as int
	return utils::toInt(elems[11]);
}

int parseModPeptide(std::string line)
{
	std::vector<std::string> elems;
	utils::split(line, IN_DELIM, elems);
	assert(elems.size() == 13);
	
	std::string sequence = elems[12];
	size_t firstP = sequence.find(".");
	size_t secP = sequence.find_last_of(".");
	sequence = sequence.substr(firstP + 1, secP - (sequence.length() - secP));
	
	//return true if any DIFMODS are found in peptide std::string.
	for(const char* p = DIFFMODS; *p; p++)
		if(utils::strContains(*p, sequence))
			return utils::toInt(elems[11]);
	
	return 0;
}

std::string getID(std::string str)
{
	size_t firstBar = str.find("|");
	size_t secBar = str.find("|", firstBar+1);
	return str.substr(firstBar+1, secBar-firstBar-1);
}

inline void Peptide::parseSequence(const std::string& str)
{
	size_t firstP = str.find(".");
	size_t secP = str.find_last_of(".");

	if(firstP == std::string::npos || secP == std::string::npos)
		calcSequence = str;
	else calcSequence = str.substr(firstP + 1, secP - (str.length() - secP));
	
	if(par->modGroupMethod == 1)
		for(const char* p = DIFFMODS; *p; p++)
			calcSequence = utils::removeChars(*p, calcSequence);
}

void Peptide::write(std::ofstream& outF)
{
	if(!outF)
		throw std::runtime_error("Bad std::ofstream!");
	
	if(!par->includeReverse)
		if(utils::strContains(REVERSE_MATCH, matchDirrection))
			return;
	
	if(!supDataAdded)
	{
		if(par->calcMW)
			calcMW();
		supDataAdded = true;
	}
	
	if(par->outputFormat == 2)
		if(!par->includeNullPeptides && col[*colIndex].isNull())
			return;
	
	outF << proteinID
	<< OUT_DELIM <<	protein
	<< OUT_DELIM <<	description
	<< OUT_DELIM <<	calcSequence
	<< OUT_DELIM <<	length;
	
	if(par->peptideGroupMethod != params::byCharge)
		outF << OUT_DELIM << charge;
	
	outF << OUT_DELIM << unique;
	
	std::streamsize ss = outF.precision();
	outF.precision(6); //write floating point nums with 5 digits
	if(par->calcMW)
		outF << OUT_DELIM << std::fixed << avgMass
		<< OUT_DELIM << std::fixed << monoMass
		<< OUT_DELIM << formula;
	outF.precision(ss);
	
	outF << OUT_DELIM << calcMH;
	
	if(par->outputFormat == 1)
	{
		assert(par->peptideSupInfoNum >= 0 && par->peptideSupInfoNum <= 2);
		if(par->supInfoOutput == 0)
		{
			for(int i = 0; i < colSize; i++)
			{
				outF << OUT_DELIM << col[i].count;
				
				if(par->includeModStat)
					outF << OUT_DELIM << col[i].modPeptidesSC;
			}
		}
		else if(par->supInfoOutput == 1)
		{
			for(int i = 0; i < colSize; i++)
				outF << OUT_DELIM << col[i].count;
			
			if(par->includeModStat)
				for(int i = 0; i < colSize; i++)
					outF << OUT_DELIM << col[i].modPeptidesSC;
		}
		
	}
	else if(par->outputFormat == 2)
	{
		if(par->peptideGroupMethod != params::byCharge)
			outF << OUT_DELIM << col[*colIndex].obsMH
			<< OUT_DELIM << col[*colIndex].xCorr
			<< OUT_DELIM << col[*colIndex].scan
			<< OUT_DELIM <<	col[*colIndex].parentFile;
		
		outF << OUT_DELIM << col[*colIndex].colname
		<< OUT_DELIM << col[*colIndex].count;
		
		if(par->includeModStat)
			outF << OUT_DELIM << col[*colIndex].modPeptidesSC;
		
		if(par->parseSampleName)
		{
			outF << OUT_DELIM << parseSample(col[*colIndex].colname, par->sampleNamePrefix, par->parseSampleName, 1)
			<< OUT_DELIM <<	parseReplicate(col[*colIndex].colname);
		}
	}
	outF << std::endl;
}

bool Peptides::writeOut(std::string ofname, const params::Params& pars)
{
	std::ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	params::OutputFormat outputFormat = pars.outputFormat;
	pars.outputFormat = params::wideFormat;
	
	bool supInfoS [] = {true, pars.includeModStat};
	bool supInfo = false;
	std::vector<std::string> supInfoHeaders;
	for(int i = 0; i < PEP_SUP_INFO_HEADERS_LEN; i++)
	{
		if(supInfoS[i])
		{
			supInfoHeaders.push_back(PEP_SUP_INFO_HEADERS[i]);
			if(i > 0)
				supInfo = true;
		}
	}
	
	//generate headers based off params
	std::vector<std::string> headers;
	std::vector<std::string>::iterator it;
	headers.insert(headers.begin(), DEFALUT_PEPTIDE_COLNAMES, utils::end(DEFALUT_PEPTIDE_COLNAMES));
	
	if(pars.peptideGroupMethod != params::byCharge)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Length(aa)")
			{
				headers.insert(it + 1, "Charge");
				break;
			}
	}
	if(pars.calcMW)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Unique")
			{
				headers.insert(it + 1, MWCALC_HEADERS, MWCALC_HEADERS + 3);
				break;
			}
	}
	
	if(!pars.sampleNamePrefix.empty())
	{
		std::string delim;
		if(pars.supInfoOutput == 0)
			delim = utils::repeat(std::string(1, OUT_DELIM), pars.peptideSupInfoNum + 1);
		else delim = std::string(1, OUT_DELIM);
		for(int i = 0; i < headers.size(); i++)
			outF << OUT_DELIM;
		for(int i = 0; i < pars.getNumFiles() ; i++)
			outF << colNames[i] << delim;
		outF << std::endl;
	}
	if(supInfo && (pars.supInfoOutput == 0))
	{
		std::string tabs = utils::repeat(std::string(1, OUT_DELIM), pars.peptideSupInfoNum + 1);
		
		std::string repeatHeaders;
		for(std::vector<std::string>::iterator it = supInfoHeaders.begin(); it != supInfoHeaders.end(); ++it)
		{
			if(it == supInfoHeaders.begin())
				repeatHeaders = *it;
			else repeatHeaders += (OUT_DELIM + *it);
		}
		
		for(int i = 0; i < headers.size(); i++)
			outF << OUT_DELIM;
		for(int i = 0; i < pars.getNumFiles(); i++)
		{
			if (i == 0)
				outF << parseSample(colNames[i], pars.sampleNamePrefix, false, 0);
			else outF << tabs << parseSample(colNames[i], pars.sampleNamePrefix, false, 0);
		}
		outF << std::endl;
		for(int i = 0; i < headers.size(); i++)
			outF << headers[i] << OUT_DELIM;
		
		for (int i = 0; i < pars.getNumFiles(); i++)
		{
			if(i == 0)
				outF << repeatHeaders;
			else outF << OUT_DELIM << repeatHeaders;
		}
		outF << std::endl;
	}
	else {
		int len = int(headers.size());
		if(pars.supInfoOutput == 1)
		{
			//assert(pars.peptideSupInfoNum >= 1 && pars.peptideSupInfoNum <= 1);
			std::string preBuffer = utils::repeat(std::string(1, OUT_DELIM), len);
			std::string postBuffer = utils::repeat(std::string(1, OUT_DELIM), colNames.size());
			
			outF << preBuffer;
			for(std::vector<std::string>::iterator it = supInfoHeaders.begin(); it != supInfoHeaders.end(); ++it)
				outF << *it << postBuffer;
			outF << std::endl;
		}
		
		std::vector<std::string> ofColNames;
		for(int i = 0; i < len; i ++)
			ofColNames.push_back(headers[i]);
		for(int i = 0; i <= pars.peptideSupInfoNum; i++)
		{
			for(int i = 0; i < pars.getNumFiles(); i ++)
				ofColNames.push_back(parseSample(colNames[i], pars.sampleNamePrefix, false, 0));
		}
		int colNamesLen = int(ofColNames.size());
		for(int i = 0; i < colNamesLen; i++)
		{
			if(i == 0)
				outF << ofColNames[i];
			else outF << OUT_DELIM << ofColNames[i];
		}
		outF << std::endl;
	}
	
	//print peptides and spectral counts
	for(DataType::iterator it = data.begin(); it != data.end(); ++it){
		it->second.write(outF);
	}
	
	pars.outputFormat = outputFormat;
	
	return true;
}


bool Peptides::writeOutDB(std::string ofname, const params::Params& pars)
{
	std::ofstream outF (ofname.c_str());
	
	if(!outF)
		return false;
	
	params::OutputFormat outputFormat = pars.outputFormat;
	pars.outputFormat = params::longFormat;
	
	//generate headers based off params
	std::vector<std::string> headers;
	std::vector<std::string>::iterator it;
	int defaultColNamesLen = DEFALUT_PEPTIDE_DB_COLNAMES_LEN;
	
	if(pars.parseSampleName)
		defaultColNamesLen += 2;
	
	for(int i = 0; i < defaultColNamesLen; i++)
		headers.push_back(DEFALUT_PEPTIDE_DB_COLNAMES[i]);
	if(pars.peptideGroupMethod != params::byCharge)
	{
		//headers to add
		std::string add [] = {"ObsMH", "xCorr", "Scan", "Parent_file"};
		for(it = headers.begin(); it != headers.end(); it++)
		{
			if(*it == "Length(aa)")
			{
				headers.insert(it + 1, "Charge");
				break;
			}
		}
		for(it = headers.begin(); it != headers.end(); it++)
		{
			if(*it == "CalcMH")
			{
				headers.insert(it + 1, utils::begin(add), utils::end(add));
				break;
			}
		}
	}
	if(pars.calcMW)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Unique")
			{
				headers.insert(it + 1, MWCALC_HEADERS, MWCALC_HEADERS + 3);
				break;
			}
	}
	if(pars.includeModStat)
	{
		for(it = headers.begin(); it != headers.end(); it++)
			if(*it == "Spectral_counts")
			{
				headers.insert(it+1, "Mod_pep_SC");
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
	outF << std::endl;
	
	Peptide::colIndex = &colIndex;
	for(DataType::iterator it = data.begin(); it != data.end(); ++it){
		for(colIndex = 0; colIndex < pars.getNumFiles(); colIndex++){
			it->second.write(outF);
		}
	}
	
	pars.outputFormat = outputFormat;
	
	return true;
}
