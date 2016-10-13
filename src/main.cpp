//
//  main.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 3/25/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "DTarray_AJM.hpp"
#include "dtafilter.cpp"
#include "subCelluarLoc.cpp"
#include "utils.cpp"
#include "calcMW.cpp"

using namespace std;

int main (int argc, char *argv[])
{
	//check arguements
	assert(argc == 4);
	string wd = string(argv[1]);
	string flistName = string(argv[2]);
	string paramsName = string(argv[3]);
	assert(util::dirExists(wd));
	
	//read in names of files to combine and output params
	FilterFileParams filterFileParams;
	if (!filterFileParams.readFlist(flistName, wd))
	{
		cout << "Failed to read file list! Exiting..." << endl;
		return 0;
	}
	if (!filterFileParams.readDTParams(paramsName, wd))
	{
		cout << "Failed to read params file! Exiting..." << endl;
		return 0;
	}

	//read in and combine files
	Proteins proteins(filterFileParams);
	if(!proteins.readIn(wd, filterFileParams))
		return 0;
	
	//read in subcellular locations database and add sub cell locations to proteins
	if(filterFileParams.getSubCelluarLoc)
	{
		if(!proteins.readInLocDB(filterFileParams.locDBfname))
		{
			cout <<"Failed to read protein location DB file! Exiting..." << endl;
			return 0;
		}
		cout << endl << "Searching for subcellular locations of proteins in dataset..." << endl;
		proteins.addSubcelluarLoc();
	}
	
	//calculate mass of peptides or proteins from sequence and amino acid mass databases
	if(filterFileParams.calcMW)
	{
		MWDB mwDB;
		if(!mwDB.readIn(wd, filterFileParams))
		{
			cout << "Failed to read mwDB files! Exiting..." << endl;
			return 0;
		}
		cout << "Calculating protein molecular weights from " << filterFileParams.mwDBFname << endl;
		proteins.calcMW(mwDB);
	}
	
	if((!filterFileParams.calcMW && filterFileParams.includeSeq ) ||
	   (filterFileParams.seqDBfname != filterFileParams.mwDBFname))
	{
		SeqDB seqDB;
		if(!seqDB.readIn(wd + filterFileParams.seqDBfname))
		{
			cout << "Failed to read seqDB file! Exiting..." << endl;
			return 0;
		}
		cout << "Getting protein sequences from " << filterFileParams.seqDBfname << endl;
		proteins.addSeq(seqDB);
	}
	
	if(!filterFileParams.sampleNamePrefix.empty())
		cout << "Parsing colnames by prefix: " << filterFileParams.sampleNamePrefix << endl;
	
	//write out combined data to OF_NAME
	if (filterFileParams.outputFormat == "standard")
	{
		if(!proteins.writeOut(wd + filterFileParams.ofname, filterFileParams))
		{
			cout << "Could not write out file! Exiting..." << endl;
			return 0;
		}
	}
	else if (filterFileParams.outputFormat == "db")
	{
		if(!proteins.writeOutDB(wd + filterFileParams.ofname, filterFileParams))
		{
			cout << "Could not write out file! Exiting..." << endl;
			return 0;
		}
		cout << endl << "Results written in database format." << endl;
	}
	
	//summarize results for user
	cout << endl << proteins.colNames.size() << " files combined." << endl;
	cout << "Results written to: " << filterFileParams.ofname << endl;
	
	return 0;
}

