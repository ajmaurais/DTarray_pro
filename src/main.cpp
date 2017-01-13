//
//  main.cpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 3/25/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "DTarray_pro.hpp"

using namespace std;

int main (int argc, char* argv[])
{
	//check arguements
	assert(argc == 4);
	string wd = string(argv[1]);
	string flistName = string(argv[2]);
	string paramsName = string(argv[3]);
	assert(utils::dirExists(wd));
	
	//read in names of files to combine and output params
	filterFile::FilterFileParams par;
	if (!par.readFlist(flistName, wd))
	{
		cout << "\nFailed to read file list! Exiting..." << endl;
		return 0;
	}
	if (!par.readDTParams(paramsName, wd))
	{
		cout << "\nFailed to read params file! Exiting..." << endl;
		return 0;
	}
	if(!par.optionsCompatable())
		return 0;
	
	if(par.peptideOutput != filterFile::none)
		cout << "Grouping peptides " <<
		filterFile::groupFormatString(par.peptideGroupMethod) << endl;
	
	Proteins proteins(par);
	Peptides peptides(par);
	
	//read saint bait file
	if(par.includeSaint)
	{
		if(!proteins.readBaitFile(par.saintBaitFile))
		{
			cout << "Could not read bait file!" << endl;
			return 0;
		}
	}
	//read in subcellular locations database and add sub cell locations to proteins
	if(par.getSubCelluarLoc)
	{
		cout << endl << "Reading subcellular locations database...";
		if(!proteins.readInLocDB(par.locDBfname))
		{
			cout <<"Failed to read protein location DB file! Exiting..." << endl;
			return 0;
		}
		cout << " done!" << endl;
	}
	
	//read in subcellular locations database and add sub cell locations to proteins
	//filterFileParams.getFxn = true;
	if(par.getFxn)
	{
		cout << endl << "Reading protein function database...";
		if(!proteins.readInFxnDB(par.fxnDBfname))
		{
			cout <<"Failed to read protein function DB file! Exiting..." << endl;
			return 0;
		}
		cout << " done!" << endl;
	}
	
	//calculate mass of peptides or proteins from sequence and amino acid mass databases
	if(par.calcMW)
	{
		cout << endl << "Calculating protein molecular weights from: " << par.mwDBFname;
		if(!proteins.readInMWdb(wd, par))
		{
			cout << "Failed to read mwDB files! Exiting..." << endl;
			return 0;
		}
		cout << " done!" << endl;
	}
	
	if((!par.calcMW && par.getSeq ) || (par.seqDBfname != par.mwDBFname))
	{
		cout << endl << "Getting protein sequences from " << par.seqDBfname;
		if(!proteins.readInSeqDB(par.seqDBfname))
		{
			cout << "Failed to read seqDB file! Exiting..." << endl;
			return 0;
		}
		cout << " done!" << endl;
	}
	
	if(!par.includeReverse)
		cout << endl << "Skipping reverse matches..." << endl;
	
	//read in and combine files
	cout << endl;
	if(!proteins.readIn(wd, par, &peptides))
		return 0;
	
	cout << endl << par.numFiles << " files combined." << endl;
	
	if(!par.sampleNamePrefix.empty())
		cout << endl << "Parsing colnames by prefix: " << par.sampleNamePrefix << endl;
	
	//write out combined protein data
	if(par.includeProteins)
	{
		assert(par.outputFormat != 0);
		if (par.outputFormat == 1 || par.outputFormat == 3)
		{
			cout << endl << "Writing protein data...";
			if(!proteins.writeOut(wd + par.ofname, par))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return 0;
			}
			cout << " done!" << endl << "Protein data written in wide format to: "
			<< par.ofname << endl << endl;
		}
		if(par.outputFormat == 2 || par.outputFormat == 3)
		{
			cout << endl << "Writing protein data...";
			if(!proteins.writeOutDB(wd + par.dbOfname, par))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return 0;
			}
			cout << " done!" << endl << "Protein data written in long format to: "
			<< par.dbOfname << endl << endl;
		}
	}
	//write out combined peptide data
	if(par.includePeptides)
	{
		assert(par.peptideOutput != 0);
		if(par.peptideOutput == 1 || par.peptideOutput == 3)
		{
			cout << endl << "Writing peptide data...";
			if(!peptides.writeOut(wd + par.peptideOfFname, par))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return 0;
			}
			cout << " done!" << endl << "Peptide data written in wide format to: "
			<< par.peptideOfFname << endl << endl;
		}
		
		if(par.peptideOutput == 2 || par.peptideOutput == 3)
		{
			cout << endl << "Writing peptide data...";
			if(!peptides.writeOutDB(wd + par.dbPeptideOfFname, par))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return 0;
			}
			cout << " done!" << endl << "Peptide data written in long format to: "
			<< par.dbPeptideOfFname << endl << endl;
		}
	}
	//write saint input files
	if(par.includeSaint)
	{
		cout << endl << "Writing saint output files...";
		if(!proteins.writeSaint(wd + par.saintPreyFname, 1))
		{
			cout << "Could not write prey file! Exiting..." << endl;
			return 0;
		}
		if(!proteins.writeSaint(wd + par.saintInteractionFname, 2))
		{
			cout << "Could not write interaction file! Exiting..." << endl;
			return 0;
		}
		cout << " done!" << endl << "prey file written to: " << par.saintPreyFname << endl
		<< "interaction data written to: " << par.saintPreyFname << endl << endl;
 	}
	
	return 0;
}

