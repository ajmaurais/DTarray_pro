//
//  dtafilter.cpp
//  DTarray_AJM
// -----------------------------------------------------------------------------
// Copyright 2018 Aaron maurais
// -----------------------------------------------------------------------------
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
// -----------------------------------------------------------------------------
//

#include <dtafilter.hpp>

Dbase* Protein::_locDB = nullptr;
Dbase* Protein::_fxnDB = nullptr;
mwDB::MWDB_Protein* Protein::_mwdb = nullptr;
mwDB::SeqDB* Protein::_seqDB = nullptr;
saint::BaitFile* Protein::_baitFile = nullptr;
locReport::LocDB* Protein::_locTable = nullptr;

molFormula::Residues* Peptide::mwdb = nullptr;

//diffmod symbols to search for
const char* DIFFMODS = "*";

void Protein::consolidate(const Protein& toAdd)
{
	_col[*toAdd._colIndex]._colname = toAdd._col[*toAdd._colIndex]._colname;
	_col[*toAdd._colIndex]._count = toAdd._col[*toAdd._colIndex]._count;
	_col[*toAdd._colIndex]._uniquePeptides = toAdd._col[*toAdd._colIndex]._uniquePeptides;
	_col[*toAdd._colIndex]._coverage = toAdd._col[*toAdd._colIndex]._coverage;
	_col[*toAdd._colIndex].sequenceCount = toAdd._col[*toAdd._colIndex].sequenceCount;
	_col[*toAdd._colIndex]._modPeptides = toAdd._col[*toAdd._colIndex]._modPeptides;
	_col[*toAdd._colIndex]._modPeptidesSC = toAdd._col[*toAdd._colIndex]._modPeptidesSC;
}

void Peptide::consolidate(const Peptide& toAdd)
{
	_col[*toAdd._colIndex]._colname = toAdd._col[*toAdd._colIndex]._colname;
	_col[*toAdd._colIndex]._count += toAdd._col[*toAdd._colIndex]._count;
	_col[*toAdd._colIndex]._scan = toAdd._col[*toAdd._colIndex]._scan;
	_col[*toAdd._colIndex]._parentFile = toAdd._col[*toAdd._colIndex]._parentFile;
	_col[*toAdd._colIndex]._obsMH = toAdd._col[*toAdd._colIndex]._obsMH;
	_col[*toAdd._colIndex]._modPeptidesSC += toAdd._col[*toAdd._colIndex]._modPeptidesSC;
	_col[*toAdd._colIndex]._xCorr = toAdd._col[*toAdd._colIndex]._xCorr;
}

template <class _Tp>
void ProteinDataTemplate<_Tp>::initialize(const std::vector<_Tp>& tempColNames, size_t colNamesLen, size_t* colIndex)
{
	_col.clear();
	_col.insert(_col.begin(), tempColNames.begin(), tempColNames.end());
	_colIndex = colIndex;
	_colSize = colNamesLen;
}

void Protein::addSeq()
{
	assert(_seqDB != nullptr);
	_sequence = _seqDB->getSequence(_ID);
}

void Protein::addLoc()
{
	assert(_locDB != nullptr);
	loc = _locDB->getDat(_ID);
}

void Protein::addFxn()
{
	assert(_fxnDB != nullptr);
	fxn = _fxnDB->getDat(_ID);
}

void Protein::operator = (const Protein& p)
{
	_ID = p._ID;
	_protein = p._protein;
	_description = p._description;
	MW = p.MW;
	pI = p.pI;
	_length = p._length;
	loc = p.loc;
	fxn = p.fxn;
	fullDescription = p.fullDescription;
	_matchDirrection = p._matchDirrection;
	_col.insert(_col.begin(), p._col.begin(), p._col.end());
	_avgMass = p._avgMass;
	_monoMass = p._monoMass;
	_sequence = p._sequence;
	_supDataAdded = p._supDataAdded;
}

void Peptide::operator = (const Peptide& p)
{
	key = p.key;
	calcSequence = p.calcSequence;
	_length = p._length;
	_proteinID = p._proteinID;
	_calcMH = p._calcMH;
	_fileName = p._fileName;
	_protein = p._protein;
	_description = p._description;
	_charge = p._charge;
	unique = p.unique;
	_col.insert(_col.begin(), p._col.begin(), p._col.end());
	_avgMass = p._avgMass;
	_monoMass = p._monoMass;
	_sequence = p._sequence;
	_supDataAdded = p._supDataAdded;
}

//parse proten header line and extract desired data
bool Protein::getProteinData(std::string line)
{
	//split line by tabs
	std::vector<std::string> elems;
	utils::split(line, IN_DELIM, elems);
	assert(elems.size() == 9);
	
	//parse elem[0]
	if(!parse_matchDir_ID_Protein(elems[0])) return false;
	
	//keep fullDescription but seperate by spaces instead of tabs
	fullDescription = elems[0];
	for(int i = 1; i < 9; i++)
		fullDescription += (" " + elems[i]);
	
	//extract pI
	pI = elems[6];
	
	//extract protein length
	_length = elems[4];
	
	//extract MW
	MW = elems[5];
	
	//extract shortened protein name and description
	size_t endOfDescription = elems[8].find(" [");
	_description = elems[8].substr(0, endOfDescription);
	
	//add spectrum count, coverage and sequence count for *this protein to colname
	_col[*_colIndex]._count = std::stoi(elems[2]);
	std::string coverageTemp = elems[3];
	coverageTemp = coverageTemp.substr(0, coverageTemp.find("%"));
	_col[*_colIndex]._coverage = coverageTemp;
	_col[*_colIndex].sequenceCount = elems[1];
	
	return true;
}

bool Protein::parse_matchDir_ID_Protein(std::string str)
{
	std::vector<std::string> elems;
	utils::split(str, '|', elems);
	
	try{
		if(elems.size() == 3){
			_matchDirrection = elems.at(0);
			_ID = elems.at(1);
			
			size_t underScoreI = elems.at(2).find_last_of("_");
			_protein = elems.at(2).substr(0, underScoreI);
		}
		else{
			_matchDirrection = elems.at(0);
			_ID = elems.at(0);
			_protein = elems.at(0);
		}
	} catch(std::out_of_range& e){
		std::cerr << "\n Error parsing protein id for " << str <<"\n Skipping...\n";
		return false;
	}
	
	if(utils::strContains(REVERSE_MATCH, _matchDirrection))
		_ID = "reverse_" + _ID;
	return true;
}

void Protein::calcMW()
{
	std::string tmp_sequence = _mwdb->seqDB->getSequence(_ID);
	
	if(tmp_sequence == mwDB::SEQ_NOT_FOUND)
	{
		_sequence = mwDB::SEQ_NOT_FOUND;
		_avgMass = -1;
		_monoMass = -1;
		_formula = "NA";
		return;
	}
	
	_sequence = tmp_sequence;
	_avgMass = _mwdb->calcMW(_sequence);
	_monoMass = _mwdb->calcMono(_sequence);
	_formula = _mwdb->calcFormula(_sequence, _par->unicode);
}

void Peptide::calcMW()
{
	_avgMass = mwdb->calcMW(calcSequence);
	_monoMass = mwdb->calcMono(calcSequence);
	_formula = mwdb->calcFormula(calcSequence, _par->unicode);
}

//loop through protein header and peptide lines in DTA filter file and add data to Proteins
bool Proteins::readIn(params::Params* const pars,
					  const std::vector <SampleData_protein>& colNamesTemp,
					  const std::vector<SampleData_peptide>& pColNamesTemp,
					  Peptides * const peptides)
{
	std::string fname = pars->getwd() + pars->getFilePath(_colIndex);
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
				Protein newProtein(pars, &_locDB, &_fxnDB, &_mwdb, &_seqDB, &_baitFile, &_locTable);
				newProtein.initialize(colNamesTemp, colNamesLen, &_colIndex);
				
				//parse protein header line and skip if error
				if(!newProtein.getProteinData(line)) continue;
				inProtein = true; //used to determine if it is valid to loop through peptide lines below protein
								  //header line to extract unique peptide spectral counts
					
				//insert new protein into proteins list
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
							Peptide newPeptide(pars, &peptides->_mwdb);
							newPeptide.initialize(pColNamesTemp, colNamesLen, &peptides->_colIndex);
							newPeptide._proteinID = newProtein._ID;
							newPeptide._protein = newProtein._protein;
							newPeptide._description = newProtein._description;
							newPeptide._matchDirrection = newProtein._matchDirrection;
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
								proteinIndex->second._col[_colIndex]._uniquePeptides += parsePeptideSC(line);
						}
						if(pars->includeModStat)
						{
							int modPeptideSC = parseModPeptide(line);
							if(modPeptideSC > 0)
							{
								proteinIndex->second._col[_colIndex]._modPeptidesSC += modPeptideSC;
								proteinIndex->second._col[_colIndex]._modPeptides++;
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
	_length = std::to_string(calcSequence.length());
	_sequence = elems[12];
	_col[*_colIndex]._obsMH = elems[5];
	_col[*_colIndex]._count = std::stoi(elems[11]);
	_col[*_colIndex]._xCorr = elems[2];
	
	if(_par->includeModStat)
		_col[*_colIndex]._modPeptidesSC += parseModPeptide(line);
	
	_calcMH = elems[6];
	_fileName = elems[1];
	
	std::string temp = elems[1];
	elems.clear();
	utils::split(temp, '.', elems);
	_charge = elems[3];
	key = makeKey();
	_col[*_colIndex]._scan = elems[2];
	_col[*_colIndex]._parentFile = elems[0];
}

std::string Peptide::makeKey() const {
	switch(_par->peptideGroupMethod){
		case params::byScan : return _fileName;
			break;
		case params::byProtein : return _proteinID + "_" + calcSequence + "_" + _charge;
			break;
		case params::byCharge : return _proteinID + "_" + calcSequence;
			break;
		default : throw std::runtime_error("Invalid type");
	}
}

void Protein::clear()
{
	_col.clear();
	MW.clear();
	fullDescription.clear();
	_matchDirrection.clear();
	_ID.clear();
	_protein.clear();
	_description.clear();
	loc.clear();
	_avgMass = 0;
	_monoMass = 0;
	_sequence.clear();
}

void Peptide::clear()
{
	_col.clear();
	unique = 0;
	_calcMH.clear();
	key.clear();
	_fileName.clear();
	_avgMass = 0;
	_monoMass = 0;
	_proteinID.clear();
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
		SampleData_protein temp (_colNames[i]);
		colNamesTemp.push_back(temp);
		SampleData_peptide pTemp (_colNames[i]);
		pColNamesTemp.push_back(pTemp);
	}
	
	for(Proteins::_colIndex = 0; Proteins::_colIndex < par->getNumFiles(); Proteins::_colIndex++)
	{
		Peptides::_colIndex = Proteins::_colIndex;
		if(!readIn(par, colNamesTemp, pColNamesTemp, peptides))
		{
			std::cerr <<"Failed to read in " << par->getFilePath(Proteins::_colIndex) <<"!" << std::endl << "Exiting..." << std::endl;
			return false;
		}
		std::cerr << "Adding " << par->getFileColname(Proteins::_colIndex) << "..." << std::endl;
	}
	return true;
}

bool Proteins::readInMWdb(const params::Params& par)
{
	return _mwdb.initalize(par);
}

bool Proteins::readInSeqDB(std::string fname)
{
	return _seqDB.readIn(fname);
}

bool Proteins::readInFxnDB(std::string fname)
{
	return _fxnDB.readIn(fname);
}

bool Proteins::readInLocDB(std::string fname)
{
	return _locDB.readIn(fname);
}

bool Proteins::readBaitFile(std::string fname)
{
	return _baitFile.read();
}

bool Peptides::readInMWdb(const params::Params& par)
{
	return _mwdb.initalize(par.atomCountTableFname,
							par.atomMassTableFname);
}

locReport::LocDat Protein::toLocDat(const std::string& _loc) const
{
	locReport::LocDat newLocDat(_loc, _matchDirrection, _colSize, _par);
	
	for(size_t i = 0; i < _colSize; i++)
		newLocDat.initializeCol(i, _col[i]._colname, (_col[i]._count > 0), _col[i]._count,
							_col[i]._uniquePeptides, std::stoi(_col[i].sequenceCount));
	
	return newLocDat;
}

void Protein::addLocToTable()
{
	if(!_supDataAdded)
		addSupData();
	
	std::vector<std::string> elems;
	utils::split(loc, DB_DELIM, elems);
	
	for(std::vector<std::string>::iterator it = elems.begin(); it != elems.end(); ++it)
		_locTable->addLoc(toLocDat(*it));
}

void Protein::writePrey(std::ofstream& outF) const
{
	assert(outF);
	outF << _ID << OUT_DELIM <<
	_length << OUT_DELIM << _description << std::endl;
}

void Protein::writeInteractions(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < _colSize; i++)
	{
		if(_col[i]._count > 0)
		{
			outF << _col[i]._colname << OUT_DELIM <<
			_baitFile->getBaitName(_col[i]._colname) << OUT_DELIM <<
			_ID << OUT_DELIM <<
			_col[i]._count << std::endl;
		}
	}
}

void Protein::addSupData()
{
	if(_par->calcMW)
		calcMW();
	if(_par->getSubCelluarLoc)
		addLoc();
	if(_par->getSeq && !_par->calcMW)
		addSeq();
	if(_par->getFxn)
		addFxn();
	_supDataAdded = true;
}

void Protein::writeProtein(std::ofstream& outF)
{
	if(!outF)
		throw std::runtime_error("Bad std::ofstream");
	
	if(!_par->includeReverse)
		if(utils::strContains(REVERSE_MATCH, _matchDirrection))
			return;
	
	if(!_supDataAdded)
		addSupData();
	
	if(INCLUDE_FULL_DESCRIPTION)
		outF << fullDescription << OUT_DELIM;
	
	outF << _ID <<
	OUT_DELIM << _protein <<
	OUT_DELIM << _description <<
	OUT_DELIM << pI <<
	OUT_DELIM << _length <<
	OUT_DELIM << MW;
	
	std::streamsize ss = outF.precision();
	outF.precision(5); //write floating point nums with 5 digits
	if(_par->calcMW)
		outF << OUT_DELIM << std::fixed << _avgMass <<
		OUT_DELIM << std::fixed << _monoMass <<
		OUT_DELIM << _formula;
	outF.precision(ss);
	
	if(_par->getSeq)
		outF << OUT_DELIM << _sequence;
	
	if(_par->getFxn)
		outF << OUT_DELIM << fxn;
	
	if(_par->getSubCelluarLoc)
		outF << OUT_DELIM << loc;
	
	if(_par->outputFormat == 1)
	{
		assert(_par->supInfoNum >= 0 && _par->supInfoNum <= 5);
		if(_par->supInfoOutput == 0)
		{
			for(int i = 0; i < _colSize; i++)
			{
				outF << OUT_DELIM << _col[i]._count;
				if (_par->includeUnique)
					outF << OUT_DELIM << _col[i]._uniquePeptides;
				
				if(_par->includeCoverage)
					outF << OUT_DELIM << _col[i]._coverage;
				
				if(_par->includeSequenceCount)
					outF << OUT_DELIM << _col[i].sequenceCount;
				
				if(_par->includeModStat)
					outF << OUT_DELIM << _col[i]._modPeptides
						<< OUT_DELIM << _col[i]._modPeptidesSC;
			}
		}
		else if(_par->supInfoOutput == 1)
		{
			writeCount(outF);
			
			if(_par->includeUnique)
				writeUnique(outF);
			
			if(_par->includeCoverage)
				writeCoverage(outF);
			
			if(_par->includeSequenceCount)
				writeSequenceCount(outF);
			
			if(_par->includeModStat)
				writeModStat(outF);
		}
	}
	else if(_par->outputFormat == 2)
	{
		outF << OUT_DELIM << _col[*_colIndex]._colname;
		
		if (_par->parseSampleName)
		{
			outF << OUT_DELIM << parseSample(_col[*_colIndex]._colname, _par->sampleNamePrefix, _par->parseSampleName, 1)
			<< OUT_DELIM << parseReplicate(_col[*_colIndex]._colname);
		}
		
		outF << OUT_DELIM << _col[*_colIndex]._count;
		
		if (_par->includeUnique)
			outF << OUT_DELIM << _col[*_colIndex]._uniquePeptides;
		
		if(_par->includeCoverage)
			outF << OUT_DELIM << _col[*_colIndex]._coverage;
		
		if(_par->includeSequenceCount)
			outF << OUT_DELIM << _col[*_colIndex].sequenceCount;
		
		if(_par->includeModStat)
			outF << OUT_DELIM << _col[*_colIndex]._modPeptides
			<< OUT_DELIM << _col[*_colIndex]._modPeptidesSC;
	}
	outF << std::endl;
}

void Protein::writeCount(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < _colSize; i++)
		outF << OUT_DELIM << _col[i]._count;
}
void Protein::writeUnique(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < _colSize; i++)
		outF << OUT_DELIM << _col[i]._uniquePeptides;
}
void Protein::writeCoverage(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < _colSize; i++)
		outF << OUT_DELIM << _col[i]._coverage;
}
void Protein::writeSequenceCount(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < _colSize; i++)
		outF << OUT_DELIM << _col[i].sequenceCount;
}

void Protein::writeModStat(std::ofstream& outF) const
{
	assert(outF);
	for(int i = 0; i < _colSize; i++)
		outF << OUT_DELIM << _col[i]._modPeptides
			<< OUT_DELIM << _col[i]._modPeptidesSC;
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
			outF << _colNames[i] << delim;
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
				outF << parseSample(_colNames[i], par.sampleNamePrefix, false, 0);
			else outF << tabs << parseSample(_colNames[i], par.sampleNamePrefix, false, 0);
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
			std::string postBuffer = utils::repeat(std::string(1, OUT_DELIM), _colNames.size());
			
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
				ofColNames.push_back(parseSample(_colNames[i], par.sampleNamePrefix, false, 0));
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
	
	Protein::_colIndex = &_colIndex;
	//print proteins and spectral counts
	for(DataType::iterator it = data.begin(); it != data.end(); ++it){
		for(_colIndex = 0; _colIndex < par.getNumFiles(); _colIndex++){
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
	Protein::_locTable = &_locTable;
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
	
	_locTable.writeLocReport(outF, 1);
	
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
			outF << _colNames[i] << delim;
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
			outF << parseSample(_colNames[i], pars.sampleNamePrefix, false, 0);
			else outF << tabs << parseSample(_colNames[i], pars.sampleNamePrefix, false, 0);
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
		std::string postBuffer = utils::repeat(std::string(1, OUT_DELIM), _colNames.size());
		
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
				ofColNames.push_back(parseSample(_colNames[i], pars.sampleNamePrefix, false, 0));
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
	
	_locTable.writeLocReport(outF, 0);
	
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
	return std::stoi(elems[11]);
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
			return std::stoi(elems[11]);
	
	return 0;
}

inline void Peptide::parseSequence(const std::string& str)
{
	size_t firstP = str.find(".");
	size_t secP = str.find_last_of(".");

	if(firstP == std::string::npos || secP == std::string::npos)
		calcSequence = str;
	else calcSequence = str.substr(firstP + 1, secP - (str.length() - secP));
	
	if(_par->modGroupMethod == 1)
		for(const char* p = DIFFMODS; *p; p++)
			calcSequence = utils::removeChars(*p, calcSequence);
}

void Peptide::write(std::ofstream& outF)
{
	if(!outF)
		throw std::runtime_error("Bad std::ofstream!");
	
	if(!_par->includeReverse)
		if(utils::strContains(REVERSE_MATCH, _matchDirrection))
			return;
	
	if(!_supDataAdded)
	{
		if(_par->calcMW)
			calcMW();
		_supDataAdded = true;
	}
	
	if(_par->outputFormat == 2)
		if(!_par->includeNullPeptides && _col[*_colIndex].isNull())
			return;
	
	outF << _proteinID
	<< OUT_DELIM <<	_protein
	<< OUT_DELIM <<	_description
	<< OUT_DELIM << _sequence
	<< OUT_DELIM <<	calcSequence
	<< OUT_DELIM <<	_length;
	
	if(_par->peptideGroupMethod != params::byCharge)
		outF << OUT_DELIM << _charge;
	
	outF << OUT_DELIM << unique;
	
	std::streamsize ss = outF.precision();
	outF.precision(6); //write floating point nums with 5 digits
	if(_par->calcMW)
		outF << OUT_DELIM << std::fixed << _avgMass
		<< OUT_DELIM << std::fixed << _monoMass
		<< OUT_DELIM << _formula;
	outF.precision(ss);
	
	outF << OUT_DELIM << _calcMH;
	
	if(_par->outputFormat == 1)
	{
		assert(_par->peptideSupInfoNum >= 0 && _par->peptideSupInfoNum <= 2);
		if(_par->supInfoOutput == 0)
		{
			for(int i = 0; i < _colSize; i++)
			{
				outF << OUT_DELIM << _col[i]._count;
				
				if(_par->includeModStat)
					outF << OUT_DELIM << _col[i]._modPeptidesSC;
			}
		}
		else if(_par->supInfoOutput == 1)
		{
			for(int i = 0; i < _colSize; i++)
				outF << OUT_DELIM << _col[i]._count;
			
			if(_par->includeModStat)
				for(int i = 0; i < _colSize; i++)
					outF << OUT_DELIM << _col[i]._modPeptidesSC;
		}
		
	}
	else if(_par->outputFormat == 2)
	{
		if(_par->peptideGroupMethod != params::byCharge)
			outF << OUT_DELIM << _col[*_colIndex]._obsMH
			<< OUT_DELIM << _col[*_colIndex]._xCorr
			<< OUT_DELIM << _col[*_colIndex]._scan
			<< OUT_DELIM <<	_col[*_colIndex]._parentFile;
		
		outF << OUT_DELIM << _col[*_colIndex]._colname
		<< OUT_DELIM << _col[*_colIndex]._count;
		
		if(_par->includeModStat)
			outF << OUT_DELIM << _col[*_colIndex]._modPeptidesSC;
		
		if(_par->parseSampleName)
		{
			outF << OUT_DELIM << parseSample(_col[*_colIndex]._colname, _par->sampleNamePrefix, _par->parseSampleName, 1)
			<< OUT_DELIM <<	parseReplicate(_col[*_colIndex]._colname);
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
			outF << _colNames[i] << delim;
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
				outF << parseSample(_colNames[i], pars.sampleNamePrefix, false, 0);
			else outF << tabs << parseSample(_colNames[i], pars.sampleNamePrefix, false, 0);
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
			std::string postBuffer = utils::repeat(std::string(1, OUT_DELIM), _colNames.size());
			
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
				ofColNames.push_back(parseSample(_colNames[i], pars.sampleNamePrefix, false, 0));
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
	
	Peptide::_colIndex = &_colIndex;
	for(DataType::iterator it = data.begin(); it != data.end(); ++it){
		for(_colIndex = 0; _colIndex < pars.getNumFiles(); _colIndex++){
			it->second.write(outF);
		}
	}
	
	pars.outputFormat = outputFormat;
	
	return true;
}
