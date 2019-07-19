//
//  dtafilter.cpp
//  DTarray_pro
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

base::StringMap* Protein::_locDB = nullptr;
base::StringMap* Protein::_fxnDB = nullptr;
utils::Residues* Protein::_mwdb = nullptr;
utils::FastaFile* Protein::_seqDB = nullptr;
saint::BaitFile* Protein::_baitFile = nullptr;
locReport::LocDB* Protein::_locTable = nullptr;

utils::Residues* Peptide::_mwdb = nullptr;
utils::FastaFile* Peptide::_seqDB = nullptr;

//diffmod symbols to search for
const char* DIFFMODS = "*";

void Protein::consolidate(const Protein& rhs)
{
	_col[*rhs._colIndex]._colname = rhs._col[*rhs._colIndex]._colname;
	_col[*rhs._colIndex]._count = rhs._col[*rhs._colIndex]._count;
	_col[*rhs._colIndex]._uniquePeptides = rhs._col[*rhs._colIndex]._uniquePeptides;
	_col[*rhs._colIndex]._coverage = rhs._col[*rhs._colIndex]._coverage;
	_col[*rhs._colIndex].sequenceCount = rhs._col[*rhs._colIndex].sequenceCount;
	_col[*rhs._colIndex]._modPeptides = rhs._col[*rhs._colIndex]._modPeptides;
	_col[*rhs._colIndex]._modPeptidesSC = rhs._col[*rhs._colIndex]._modPeptidesSC;
}

void Peptide::consolidate(const Peptide& rhs)
{
	_col[*rhs._colIndex]._colname = rhs._col[*rhs._colIndex]._colname;
	_col[*rhs._colIndex]._count += rhs._col[*rhs._colIndex]._count;
	_col[*rhs._colIndex]._scan = rhs._col[*rhs._colIndex]._scan;
	_col[*rhs._colIndex]._parentFile = rhs._col[*rhs._colIndex]._parentFile;
	_col[*rhs._colIndex]._obsMH = rhs._col[*rhs._colIndex]._obsMH;
	_col[*rhs._colIndex]._modPeptidesSC += rhs._col[*rhs._colIndex]._modPeptidesSC;
	_col[*rhs._colIndex]._xCorr = rhs._col[*rhs._colIndex]._xCorr;
	_mods.insert(rhs._mods.begin(), rhs._mods.end());
}

template <class _Tp>
void base::ProteinDataTemplate<_Tp>::initialize(const std::vector<_Tp>& tempColNames,
												size_t colNamesLen, size_t* colIndex)
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
	auto it = _locDB->find(_ID);
	if(it == _locDB->end()){
		loc = DAT_NOT_FOUND;
	}
	else loc = it->second;
}

void Protein::addFxn()
{
	assert(_fxnDB != nullptr);
	auto it = _fxnDB->find(_ID);
	if(it == _fxnDB->end()){
		fxn = DAT_NOT_FOUND;
	}
	else fxn = it->second;
}

void Protein::operator = (const Protein& rhs)
{
	_ID = rhs._ID;
	_protein = rhs._protein;
	_description = rhs._description;
	MW = rhs.MW;
	pI = rhs.pI;
	_length = rhs._length;
	loc = rhs.loc;
	fxn = rhs.fxn;
	fullDescription = rhs.fullDescription;
	_matchDirrection = rhs._matchDirrection;
	_col.insert(_col.begin(), rhs._col.begin(), rhs._col.end());
	_avgMass = rhs._avgMass;
	_monoMass = rhs._monoMass;
	_sequence = rhs._sequence;
	_supDataAdded = rhs._supDataAdded;
}

void Peptide::operator = (const Peptide& rhs)
{
	_key = rhs._key;
	_calcSequence = rhs._calcSequence;
	_length = rhs._length;
	_proteinID = rhs._proteinID;
	_matchDirrection = rhs._matchDirrection;
	_calcMH = rhs._calcMH;
	_fileName = rhs._fileName;
	_protein = rhs._protein;
	_description = rhs._description;
	_charge = rhs._charge;
	unique = rhs.unique;
	_compareSequence = rhs._compareSequence;
	_baseSequence = rhs._baseSequence;
	_mods = rhs._mods;
	_begin = rhs._begin;
	_end = rhs._end;
	_col.insert(_col.begin(), rhs._col.begin(), rhs._col.end());
	_avgMass = rhs._avgMass;
	_monoMass = rhs._monoMass;
	_sequence = rhs._sequence;
	_supDataAdded = rhs._supDataAdded;
}

//parse proten header line and extract desired data
bool Protein::getProteinData(std::string line)
{
	//split line by tabs
	std::vector<std::string> elems;
	utils::split(line, IN_DELIM, elems);
	utils::removeEmptyStrings(elems);
	assert(elems.size() == 9);
	
	//parse elem[0]
	if(!parse_matchDir_ID_Protein(elems[0])) return false;
	
	//keep fullDescription but separate by spaces instead of tabs
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
	
	//check _par->filter is true
	if(_par->getFilter())
	{
		//transform subject str to lowercase
		std::string temp;
		if(_par->getToLower()) temp = utils::toLower(_description);
		else temp = _description;
		
		//if(utils::strContains("keratin", temp))
		//	std::cout << "Found!\n";
		
		//if exclude
		if(!_par->getExcludeStr().empty())
		{
			if(_par->getMatchRegex())
			{
				//bool match = std::regex_search(temp, std::regex(_par->getExcludeStr()));
				if(std::regex_search(temp, std::regex(_par->getExcludeStr())))
				{
					std::cout << "Skipping: " << _description << std::endl;
					return false;
				}
			}
			else {
				if(utils::strContains(_par->getExcludeStr(), temp))
				{
					std::cout << "Skipping: " << temp << std::endl;
					return false;
				}
			} //end else
		}// end if exclude
		
		//if add
		if(!_par->getAddStr().empty())
		{
			if(_par->getMatchRegex())
			{
				//bool match = std::regex_search(temp, std::regex(_par->getAddStr()));
				if(!std::regex_search(temp, std::regex(_par->getAddStr()))) return false;
			}
			else {
				if(!utils::strContains(_par->getExcludeStr(), temp)) return false;
			} //end else
		} //end if add
	}//end if filter
	
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
			_matchDirrection = utils::toLower(elems.at(0));
			_ID = elems.at(1);
			
			size_t underScoreI = elems.at(2).find_last_of("_");
			_protein = elems.at(2).substr(0, underScoreI);
		}
		else{
			_matchDirrection = utils::toLower(elems.at(0));
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
	std::string tmp_sequence = _seqDB->getSequence(_ID);
	
	if(tmp_sequence == utils::PROT_SEQ_NOT_FOUND)
	{
		_sequence = utils::PROT_SEQ_NOT_FOUND;
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
	_avgMass = _mwdb->calcMW(_calcSequence);
	_monoMass = _mwdb->calcMono(_calcSequence);
	_formula = _mwdb->calcFormula(_calcSequence, _par->unicode);
}

/**
 \brief Calculate the _begin, _end, and _mods Peptide members. <br>
 
 Peptide::_proteinID is used to lookup parent protein sequence in
 Peptide::_seqDB. The index at which the peptide begins in the protein
 is then used to calculate the _mods and _end members.
 */
void Peptide::calcMod()
{
	std::string _protSeq = _seqDB->getSequence(_proteinID);
	if(_protSeq == utils::PROT_SEQ_NOT_FOUND){
		_begin = 0; _end = 0;
		_mods.insert(utils::PROT_SEQ_NOT_FOUND);
		return;
	}
	
	if(_baseSequence == "CLALGMSRDAVK")
		std::cout << "Found!" << NEW_LINE;
	
	//get start index of peptide sequence in protein sequence
	size_t length = _calcSequence.length();
	_begin = _protSeq.find(_baseSequence) + 1; // 1 added to _begin here
	if(_begin == std::string::npos){
		_begin = 0; _end = 0;
		_mods.insert(utils::PEP_SEQ_NOT_FOUND);
		return;
	}
	
	size_t posCount = 0;
	bool modFound;
	for(size_t i = 0; i < length; i++)
	{
		modFound = false;
		for(const char* c = DIFFMODS; *c; c++)
		{
			if(_calcSequence[i] == *c)
			{
				modFound = true;
				_mods.insert(std::string(1, _calcSequence[i - 1]) +
							 std::to_string(_begin + posCount - 1));
				break;
			}
		}
		
		if(!modFound) posCount++;
	}
	_end = _begin + posCount - 1;
}

/**
 Loop through protein header and peptide lines in DTA filter file and add data to Proteins
 */
bool Proteins::readIn(params::Params* const pars,
					  const std::vector <base::SampleData_protein>& colNamesTemp,
					  const std::vector<base::SampleData_peptide>& pColNamesTemp,
					  Peptides& peptides)
{
	std::string fname = pars->getwd() + pars->getFilePath(_colIndex);
	
	std::ifstream inF(fname);
	if(!inF) return false;
	
	DataType::iterator proteinIndex;
	bool inProtein = false;
	bool foundHeader = false;
	std::string line;
	std::map<std::string, Peptide>::iterator peptidesIndex;
	size_t colNamesLen = pars->getNumFiles();
	
	while(utils::safeGetline(inF, line))
	{
		if(utils::trim(line).empty())
			continue;
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
				//initialize Protein to hold data for current line
				Protein newProtein;
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
					std::streampos pos;
					while(utils::safeGetline(inF, line, pos))
					{
						if(utils::trim(line).empty())
							continue;
						
						//break if starting new protein or end of file
						if(utils::strContains('%', line) || line == "\tProteins\tPeptide IDs\tSpectra")
							break;
						
						if(pars->includePeptides)
						{
							Peptide newPeptide;
							newPeptide.initialize(pColNamesTemp, colNamesLen, &peptides._colIndex);
							newPeptide._proteinID = newProtein._ID;
							newPeptide._protein = newProtein._protein;
							newPeptide._description = newProtein._description;
							newPeptide._matchDirrection = newProtein._matchDirrection;
							newPeptide._parsePeptide(line);
							
							if(pars->peptideGroupMethod == params::Params::byScan){
								peptides.data[newPeptide._key] = newPeptide;
							}
							else{
								peptidesIndex = peptides.data.find(newPeptide._key);
								if(peptidesIndex == peptides.data.end()){
									peptides.data[newPeptide._key] = newPeptide;
								}
								else{
									peptides.data[newPeptide._key].consolidate(newPeptide);
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
					inF.seekg(pos);
					inProtein = false;
				}//end of if
			}//end of else
		}//end of if
	}//end of while
	return true;
}

void Peptide::_parsePeptide(const std::string& line)
{
	std::vector<std::string> elems;
	utils::split(line, IN_DELIM, elems);
	assert(elems.size() == 13);
	
	unique = elems[0] == "*";
	_parseSequence(elems[12]); //get calcSequence and compare sequence
	_length = std::to_string(_calcSequence.length());
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
	_key = makeKey();
	_col[*_colIndex]._scan = elems[2];
	_col[*_colIndex]._parentFile = elems[0];
}

std::string Peptide::makeKey() const {
	switch(_par->peptideGroupMethod){
		case params::Params::byScan : return _fileName;
			break;
		case params::Params::byProtein : return _proteinID + "_" + _compareSequence + "_" + _charge;
			break;
		case params::Params::byCharge : return _proteinID + "_" + _compareSequence;
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
	_key.clear();
	_fileName.clear();
	_avgMass = 0;
	_monoMass = 0;
	_proteinID.clear();
}

/**
 Public Proteins::readIn function which adds all files in filterFileParams to Proteins
 and summarizes progress for user.
 
 \param par Initialized Params object.
 \param peptides Peptides object to add peptide data to.
 
 \return true if I/O was suecessful.
 */
bool Proteins::readIn(params::Params* const par, Peptides& peptides)
{
	//initialize static Protein variables
	Protein::_par = par;
	Protein::_locDB = &_locDB;
	Protein::_fxnDB = &_fxnDB;
	Protein::_mwdb = &_mwdb;
	Protein::_seqDB = &_seqDB;
	Protein::_baitFile = &_baitFile;
	Protein::_locTable = &_locTable;
	
	//initialize statid Peptide variables
	Peptide::_par = par;
	Peptide::_mwdb = &peptides._mwdb;
	Peptide::_seqDB = peptides._seqDB;
	
	size_t colNamesLen = par->getNumFiles();
	
	std::vector<base::SampleData_protein> colNamesTemp;
	std::vector<base::SampleData_peptide> pColNamesTemp;
	
	for(int i = 0; i < colNamesLen; i++)
	{
		base::SampleData_protein temp (_colNames[i]);
		colNamesTemp.push_back(temp);
		base::SampleData_peptide pTemp (_colNames[i]);
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
		std::cout << "Adding " << par->getFileColname(Proteins::_colIndex) << "..." << std::endl;
	}
	return true;
}

bool Proteins::readInMWdb(const params::Params& par)
{
	return _mwdb.initialize(par.atomCountTableFname,
							par.atomMassTableFname);
}

bool Proteins::readInSeqDB(std::string fname)
{
	return _seqDB.read(fname);
}

/**
 Read data from tsv file into base::StringMap.
 
 \param fname File name of .tsv file.
 \param dat Empty map to populate.
 \param columns Array where first element key column name and second is value name.
 \return true if successful, false if \p fname can't be read, or \p columns to not
 exist in \p fname.
 */
bool Proteins::_readDB(std::string fname,
					   base::StringMap& dat,
					   const std::string columns [2])
{
	dat.clear();
	//read tsv file
	utils::TsvFile tsvFile;
	if(!tsvFile.read(fname)) return false;
	
	//make sure required columns exist
	for(int i = 0; i < 2; i++){
		if(!tsvFile.colExists(columns[i])){
			std::cerr << "Required column: " << columns[i]
			<< " does not exist in " << fname << NEW_LINE;
			return false;
		}
	}
	
	size_t nRow = tsvFile.getNrow();
	for(size_t i = 0; i < nRow; i++){
		dat[tsvFile.getValStr(i, columns[0])] = tsvFile.getValStr(i, columns[1]);
	}
	
	return true;
}

bool Proteins::readInFxnDB(std::string fname)
{
	std::string cols [] = {"id", "panther_category"};
	
	return _readDB(fname, _fxnDB, cols);
}

/**
 Read data from csv file containing protein subcelluar localizations.
 
 \param fname Path to file.
 \param col Name of column to get data from.
 
 \return true if sucessful.
 */
bool Proteins::readInLocDB(std::string fname, std::string col)
{
	std::string cols [] = {"id", col};
	
	return _readDB(fname, _locDB, cols);
}

bool Proteins::readBaitFile(std::string fname)
{
	return _baitFile.read(fname);
}

bool Peptides::readInMWdb(const params::Params& par)
{
	return _mwdb.initialize(par.atomCountTableFname,
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

void Protein::writePrey(std::ostream& outF) const
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

void Protein::writeProtein(std::ostream& outF)
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
			outF << OUT_DELIM << parseSample(_col[*_colIndex]._colname, _par->sampleNamePrefix,
											 _par->parseSampleName, 1, _par->getMatchRegex())
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

void Protein::writeCount(std::ostream& outF) const
{
	assert(outF);
	for(int i = 0; i < _colSize; i++)
		outF << OUT_DELIM << _col[i]._count;
}
void Protein::writeUnique(std::ostream& outF) const
{
	assert(outF);
	for(int i = 0; i < _colSize; i++)
		outF << OUT_DELIM << _col[i]._uniquePeptides;
}
void Protein::writeCoverage(std::ostream& outF) const
{
	assert(outF);
	for(int i = 0; i < _colSize; i++)
		outF << OUT_DELIM << _col[i]._coverage;
}
void Protein::writeSequenceCount(std::ostream& outF) const
{
	assert(outF);
	for(int i = 0; i < _colSize; i++)
		outF << OUT_DELIM << _col[i].sequenceCount;
}

void Protein::writeModStat(std::ostream& outF) const
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
	
	params::Params::OutputFormat outputFormat = par.outputFormat;
	par.outputFormat = params::Params::wideFormat;
	
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
	if(par.getSubCelluarLoc){
		it = std::find(headers.begin(), headers.end(), "Mass(Da)");
		headers.insert(it + 1, "subcelluar_loc");
	}
	if(par.getFxn){
		it = std::find(headers.begin(), headers.end(), "Mass(Da)");
		headers.insert(it + 1, "Function");
	}
	if(par.calcMW){
		int mwHeadersLen = MWCALC_HEADERS_LENGTH;
		if(par.getSeq)
			mwHeadersLen++;
		
		it = std::find(headers.begin(), headers.end(), "Mass(Da)");
		headers.insert(it + 1, MWCALC_HEADERS, MWCALC_HEADERS + mwHeadersLen);
		
	}
	else if(par.getSeq && !par.calcMW){
		it = std::find(headers.begin(), headers.end(), "Mass(Da)");
		headers.insert(it + 1, MWCALC_HEADERS[3]);
	}
	
	colNamesLength = int(headers.size());
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
				outF << parseSample(_colNames[i], par.sampleNamePrefix, false, 0, par.getMatchRegex());
			else outF << tabs << parseSample(_colNames[i], par.sampleNamePrefix, false, 0, par.getMatchRegex());
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
				ofColNames.push_back(parseSample(_colNames[i], par.sampleNamePrefix, false, 0, par.getMatchRegex()));
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
	
	params::Params::OutputFormat outputFormat = par.outputFormat;
	par.outputFormat = params::Params::longFormat;
	
	//print column headers
	std::vector<std::string> headers;
	std::vector<std::string>::iterator it = headers.begin();
	
	for (int i = 0; i < DEFAULT_COL_NAMES_DB_LENGTH ; i++)
		headers.push_back(DEFAULT_COL_NAMES_DB[i]);
	
	if(par.parseSampleName)
	{
		it = std::find(headers.begin(), headers.end(), "Long_sample_name");
		headers.insert(it + 1, PARSE_SAMPLE_NAME_HEADERS, PARSE_SAMPLE_NAME_HEADERS + PARSE_SAMPLE_NAME_HEADERS_LEN);
	}
	if(par.getSubCelluarLoc)
	{
		it = std::find(headers.begin(), headers.end(), "Mass(Da)");
		headers.insert(it + 1, "Subcelluar_loc");
	}
	if(par.getFxn)
	{
		it = std::find(headers.begin(), headers.end(), "Mass(Da)");
		headers.insert(it + 1, "Function");
	}
	if(par.calcMW)
	{
		int mwHeadersLen = MWCALC_HEADERS_LENGTH;
		if(par.getSeq)
			mwHeadersLen++;
		
		it = std::find(headers.begin(), headers.end(), "Mass(Da)");
		headers.insert(it + 1, MWCALC_HEADERS, MWCALC_HEADERS + mwHeadersLen);
	}
	else if(par.getSeq && !par.calcMW)
	{
		it = headers.begin();
		it = std::find(headers.begin(), headers.end(), "Mass(Da)");
		headers.insert(it + 1, MWCALC_HEADERS[3]);
	}
	if(par.includeModStat)
	{
		it = std::find(headers.begin(), headers.end(), "Spectral_counts");
		headers.insert(it + 1, SUP_INFO_HEADERS[5]);
		headers.insert(it + 1, SUP_INFO_HEADERS[4]);
	}
	if(par.includeSequenceCount)
	{
		it = std::find(headers.begin(), headers.end(), "Spectral_counts");
		headers.insert(it + 1, SUP_INFO_HEADERS[3]);
	}
	if(par.includeCoverage)
	{
		it = std::find(headers.begin(), headers.end(), "Spectral_counts");
		headers.insert(it + 1, SUP_INFO_HEADERS[2]);
	}
	if(par.includeUnique)
	{
		it = std::find(headers.begin(), headers.end(), "Spectral_counts");
		headers.insert(it + 1, SUP_INFO_HEADERS[1]);
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

/**
 Populate Protein::_locTable with localization data.
 
 \param summary Should sub-organell locations be omitted from loc report?
 */
void Proteins::buildLocTable(bool summary)
{
	//init _locTable if necissary
	if(summary) _locTable.set_summaryLocs();
	
	//apply loc build fxn across protein data
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
	
	headers.insert(headers.end(), std::begin(LOC_REPORT_HEADERS), std::end(LOC_REPORT_HEADERS));
	
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
	
	repeatedHeadersV.insert(repeatedHeadersV.end(), std::begin(LOC_REPORT_HEADERS), std::end(LOC_REPORT_HEADERS));
	
	if(pars.includeUnique){
		it = std::find(repeatedHeadersV.begin(), repeatedHeadersV.end(), "Sum_SC");
		repeatedHeadersV.insert(it + 1, "Sum_uniq_SC");
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
			outF << parseSample(_colNames[i], pars.sampleNamePrefix, false, 0, pars.getMatchRegex());
			else outF << tabs << parseSample(_colNames[i], pars.sampleNamePrefix, false, 0, pars.getMatchRegex());
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
				ofColNames.push_back(parseSample(_colNames[i], pars.sampleNamePrefix, false, 0, pars.getMatchRegex()));
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

/**
 optional fxn to parse long sample name
 */
std::string parseSample(std::string sampleName, std::string prefix,
						bool parseSampleName, bool outputFormat, bool re)
{
	 std::regex pattern = std::regex(prefix);
	
	 //return unparsed sampleName if prefix is empty string or is not found in sampleName
	 if(!parseSampleName && (!(re ? std::regex_search(sampleName, pattern) :
							   utils::strContains(prefix, sampleName)) ||
							 prefix.length() == 0)){
		 return sampleName;
	 }
	 else if(parseSampleName && prefix.empty()){
		 return sampleName.substr(0, sampleName.find_last_of("_"));
	 }
	 else {
		 std::string sample = re ? std::regex_replace(sampleName, pattern, std::string("")) :
								   utils::removeSubstr(prefix, sampleName); //remove prefix from sampleName
		 return outputFormat ? sample.substr(0, sample.find_last_of("_")) : sample;
	 }
}

//!get replicate number from sample name
std::string parseReplicate(std::string str)
{
	return str.substr(str.find_last_of("_")+1);
}

/**
 Get number of spectral counts for a peptide from a line in DTA filter file
 */
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

/**
 Generate _calcSequence, _baseSequence, and _compareSequence by parsing \p str.
 \param str full peptide sequence to parse
 */
void Peptide::_parseSequence(const std::string& str)
{
	size_t firstP = str.find(".");
	size_t secP = str.find_last_of(".");
	
	//calc calcSequence
	if(firstP == std::string::npos || secP == std::string::npos)
		_calcSequence = str;
	else _calcSequence = str.substr(firstP + 1, secP - (str.length() - secP));
	
	//calc compareSequence
	_baseSequence = _calcSequence;
	for(const char* p = DIFFMODS; *p; p++)
		_baseSequence = utils::removeChars(*p, _baseSequence);
	
	//get _compareSequence
	_compareSequence = _par->modGroupMethod == 0 ? _calcSequence : _baseSequence;
}

/**
 Concat Peptide::_mods into a nicely formatted string.
 
 \param delim deliminator between multiple mods
 \return mods as a single string.
 */
std::string Peptide::getMods(std::string delim) const
{
	std::string ret = "";
	for(auto it = _mods.begin(); it != _mods.end(); ++it)
	{
		if(it == _mods.begin())
			ret += *it;
		else ret += delim + *it;
	}
	return ret;
}

void Peptide::write(std::ostream& outF)
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
		if(_par->getSeq)
			calcMod();
		_supDataAdded = true;
	}
	
	if(_par->outputFormat == 2)
		if(!_par->includeNullPeptides && _col[*_colIndex].isNull())
			return;
	
	outF << _proteinID
	<< OUT_DELIM <<	_protein
	<< OUT_DELIM <<	_description
	<< OUT_DELIM << _sequence
	<< OUT_DELIM <<	_compareSequence
	<< OUT_DELIM <<	_length;
	
	if(_par->peptideGroupMethod != params::Params::byCharge)
		outF << OUT_DELIM << _charge;
	
	outF << OUT_DELIM << unique;
	
	if(_par->getSeq)
		outF << OUT_DELIM << _begin
		<< OUT_DELIM << _end
		<< OUT_DELIM << getMods();
	
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
		if(_par->peptideGroupMethod != params::Params::byCharge)
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
			outF << OUT_DELIM << parseSample(_col[*_colIndex]._colname, _par->sampleNamePrefix,
											 _par->parseSampleName, 1, _par->getMatchRegex())
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
	
	params::Params::OutputFormat outputFormat = pars.outputFormat;
	pars.outputFormat = params::Params::wideFormat;
	
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
	headers.insert(headers.begin(),
				   DEFALUT_PEPTIDE_COLNAMES,
				   std::end(DEFALUT_PEPTIDE_COLNAMES));
	
	if(pars.peptideGroupMethod != params::Params::byCharge){
		it = std::find(headers.begin(), headers.end(), "Length(aa)");
		headers.insert(it + 1, "Charge");
	}
	if(pars.calcMW){
		it = std::find(headers.begin(), headers.end(), "Unique");
		headers.insert(it + 1, MWCALC_HEADERS, MWCALC_HEADERS + 3);
	}
	if(pars.getSeq){
		it = std::find(headers.begin(), headers.end(), "Unique");
		headers.insert(it + 1, PEP_SEQ_HEADERS, PEP_SEQ_HEADERS + PEP_SEQ_HEADERS_LEN);
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
				outF << parseSample(_colNames[i], pars.sampleNamePrefix, false, 0, pars.getMatchRegex());
			else outF << tabs << parseSample(_colNames[i], pars.sampleNamePrefix, false, 0, pars.getMatchRegex());
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
				ofColNames.push_back(parseSample(_colNames[i], pars.sampleNamePrefix, false, 0, pars.getMatchRegex()));
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
	
	params::Params::OutputFormat outputFormat = pars.outputFormat;
	pars.outputFormat = params::Params::longFormat;
	
	//generate headers based off params
	std::vector<std::string> headers;
	std::vector<std::string>::iterator it;
	int defaultColNamesLen = DEFALUT_PEPTIDE_DB_COLNAMES_LEN;
	
	if(pars.parseSampleName)
		defaultColNamesLen += 2;
	
	for(int i = 0; i < defaultColNamesLen; i++)
		headers.push_back(DEFALUT_PEPTIDE_DB_COLNAMES[i]);
	if(pars.peptideGroupMethod != params::Params::byCharge)
	{
		//headers to add
		std::string add [] = {"ObsMH", "xCorr", "Scan", "Parent_file"};
		if(pars.peptideGroupMethod != params::Params::byCharge){
			it = std::find(headers.begin(), headers.end(), "Length(aa)");
			headers.insert(it + 1, "Charge");
		}
		it = std::find(headers.begin(), headers.end(), "CalcMH");
		headers.insert(it + 1, std::begin(add), std::end(add));
	}
	if(pars.calcMW){
		it = std::find(headers.begin(), headers.end(), "Unique");
		headers.insert(it + 1, MWCALC_HEADERS, MWCALC_HEADERS + 3);
	}
	if(pars.getSeq){
		it = std::find(headers.begin(), headers.end(), "Unique");
		headers.insert(it + 1, PEP_SEQ_HEADERS, PEP_SEQ_HEADERS + PEP_SEQ_HEADERS_LEN);
	}
	if(pars.includeModStat){
		it = std::find(headers.begin(), headers.end(), "Spectral_counts");
		headers.insert(it+1, "Mod_pep_SC");
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
