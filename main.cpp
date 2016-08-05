/*
 DTarray_AJM reads in a specified number of dtaselect-filter files and writes protein, their molecular
 weights and spectral count data to OF_NAME in the working dirrectory.
 
 Written by Aaron Maurais
 25 May 2016
 */

#include <iostream>
#include <vector>
#include "dtafilter.h"

using namespace std;

int main (int argc, char *argv[])
{
	//check paramaters
	assert(argc == 3);
	string wd = string(argv[1]);
	string paramsName = string(argv[2]);
	assert(dirExists(wd));
	
	//read in names of files to combine from params file
	FilterFileParams filterFileParams;
	if (!filterFileParams.readDTParams(paramsName, wd))
	{
		cout <<"Failed to read params file! Exiting..." << endl;
		return 0;
	}
	
	//combine files
	cout << endl;
	Proteins proteins(filterFileParams);
	if(!proteins.readIn(wd, filterFileParams))
		return 0;
	
	//read in subcellular locations database and add sub cell locations to proteins
	if(filterFileParams.getSubCelluarLoc)
	{
		Btree proteinDBtree;
		if(!proteinDBtree.readInProteins(filterFileParams.locDBfname))
		{
			cout <<"Failed to read protein location DB file! Exiting..." << endl;
			return 0;
		}
		cout << endl << "Searching for subcellular locations of proteins in dataset..." << endl;
		proteins.addSubcelluarLoc(proteinDBtree);
	}
	
	if(!filterFileParams.sampleNamePrefix.empty())
		cout << "Parsing colnames by prefix: " << filterFileParams.sampleNamePrefix << endl;
	
	//write out combined data to OF_NAME
	if (filterFileParams.outputFormat == "standard")
	{
		if(!proteins.writeOut(wd + OF_NAME, filterFileParams))
		{
			cout << "Could not write out file! Exiting..." << endl;
			return 0;
		}
	}
	else if (filterFileParams.outputFormat == "DB")
	{
		if(!proteins.writeOutDB(wd + OF_NAME, filterFileParams))
		{
			cout << "Could not write out file! Exiting..." << endl;
			return 0;
		}
		cout << endl << "Results written in database format." << endl;
	}
	
	//summarize results for user
	cout << endl << proteins.colNames.size() << " files combined." << endl;
	cout << "Results written to: " << OF_NAME << endl;
	
	return 0;
}

