//
//  main.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 3/25/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "DTarray_AJM.hpp"
#include "FilterFile.cpp"
#include "dtafilter.cpp"
#include "utils.cpp"
#include "calcMW.cpp"
#include "dbase.cpp"

using namespace std;

int main (int argc, char *argv[])
{
	//check arguements
	assert(argc == 4);
	string wd = string(argv[1]);
	string flistName = string(argv[2]);
	string paramsName = string(argv[3]);
	assert(utils::dirExists(wd));
	
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
	
	Proteins proteins(filterFileParams);
	Peptides peptides(filterFileParams);
	Peptides * peptidesPtr = new Peptides(peptides);

	//read in subcellular locations database and add sub cell locations to proteins
	if(filterFileParams.getSubCelluarLoc)
	{
		cout << endl << "Reading subcellular locations database..." << endl;
		if(!proteins.readInLocDB(filterFileParams.locDBfname))
		{
			cout <<"Failed to read protein location DB file! Exiting..." << endl;
			return 0;
		}
	}
	
	//read in subcellular locations database and add sub cell locations to proteins
	//filterFileParams.getFxn = true;
	if(filterFileParams.getFxn)
	{
		cout << endl << "Reading protein function database..." << endl;
		if(!proteins.readInFxnDB(filterFileParams.fxnDBfname))
		{
			cout <<"Failed to read protein function DB file! Exiting..." << endl;
			return 0;
		}
	}
	
	//calculate mass of peptides or proteins from sequence and amino acid mass databases
	if(filterFileParams.calcMW)
	{
		if(!proteins.readInMWdb(wd, filterFileParams))
		{
			cout << "Failed to read mwDB files! Exiting..." << endl;
			return 0;
		}
		cout << "Calculating protein molecular weights from: " << filterFileParams.mwDBFname << endl;
	}
	
	if((!filterFileParams.calcMW && filterFileParams.getSeq ) ||
	   (filterFileParams.seqDBfname != filterFileParams.mwDBFname))
	{
		if(!proteins.readInSeqDB(wd + filterFileParams.seqDBfname))
		{
			cout << "Failed to read seqDB file! Exiting..." << endl;
			return 0;
		}
		cout << "Getting protein sequences from " << filterFileParams.seqDBfname << endl;
	}
	
	//read in and combine files
	if(!proteins.readIn(wd, filterFileParams, peptidesPtr))
		return 0;
	
	cout << endl << filterFileParams.numFiles << " files combined." << endl;
	
	if(!filterFileParams.sampleNamePrefix.empty())
		cout << "Parsing colnames by prefix: " << filterFileParams.sampleNamePrefix << endl;
	
	//write out combined protein data
	if(filterFileParams.includeProteins)
	{
		assert(filterFileParams.includeProteins != 0);		
		if (filterFileParams.outputFormat == 1 || filterFileParams.outputFormat == 3)
		{
			if(!proteins.writeOut(wd + filterFileParams.ofname, filterFileParams))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return 0;
			}
		cout << endl << "Protein data written in wide format to: " << filterFileParams.ofname << endl;
		}
		
		if(filterFileParams.outputFormat == 2 || filterFileParams.outputFormat == 3)
		{
			if(!proteins.writeOutDB(wd + filterFileParams.dbOfname, filterFileParams))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return 0;
			}
			cout << endl << "Protein data written in long format to: " << filterFileParams.dbOfname << endl;
		}
	}
	
	//write out combined peptide data
	if(filterFileParams.includePeptides)
	{
		assert(filterFileParams.peptideOutput != 0);
		if(filterFileParams.peptideOutput == 1 || filterFileParams.peptideOutput == 3)
		{
			if(!peptides.writeOut(wd + filterFileParams.peptideOfFname, filterFileParams))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return 0;
			}
		}
		cout << endl << "peptide data written in long format to: " << filterFileParams.peptideOfFname << endl;
		
		if(filterFileParams.peptideOutput == 2 || filterFileParams.peptideOutput == 3)
		{
			if(!peptides.writeOutDB(wd + filterFileParams.dbPeptideOfFname, filterFileParams))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return 0;
			}
			cout << endl << "Peptide data written in long format to: " << filterFileParams.dbPeptideOfFname << endl;
		}
	}
	
	return 0;
}

