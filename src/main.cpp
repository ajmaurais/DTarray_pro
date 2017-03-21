//
//  main.cpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 3/25/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "DTarray_pro.hpp"

using namespace std;

int main (int argc, char* argv[])
{
	//read in names of files to combine and output params
	params::Params par;
	if(!par.getOpts(argc, argv))
		return -1;
	
	cout << endl << "DTarray_pro v" << BIN_VERSION_NUM << endl;
	
	if(!utils::fileExists(par.wd + par.flistName) || par.rewriteFlist)
	{
		assert(par.writeFlist());
		cout << endl << "Generating " << par.flistName << " using " << par.inputFormat << " input format." << endl;
	}
	
	if(!par.readFlist(par.flistName, par.wd))
	{
		cout << "\nFailed to read file list! Exiting..." << endl;
		return -1;
	}
	if(par.numFiles <= 0)
	{
		assert(par.numFiles == 0);
		cout << endl << "No filter files found." << endl << endl;
		return -1;
	}
	if(!par.optionsCompatable())
		return -1;
	
	if(par.peptideOutput != params::none)
		cout << endl << "Grouping peptides " <<
		params::groupFormatString(par.peptideGroupMethod) << "." << endl;
	
	Proteins proteins(par);
	Peptides peptides(par);
	
	//read saint bait file
	if(par.includeSaint)
	{
		if(!proteins.readBaitFile(par.saintBaitFile))
		{
			cout << endl << "Could not read bait file!" << endl;
			return -1;
		}
	}
	//read in subcellular locations database and add sub cell locations to proteins
	if(par.getSubCelluarLoc)
	{
		cout << endl << "Reading subcellular locations database...";
		if(!proteins.readInLocDB(par.locDBfname))
		{
			cout <<"Failed to read protein location DB file! Exiting..." << endl;
			return -1;
		}
		cout << " done!" << endl;
	}
	
	//read in subcellular locations database and add sub cell locations to proteins
	if(par.getFxn)
	{
		cout << endl << "Reading protein function database...";
		if(!proteins.readInFxnDB(par.fxnDBfname))
		{
			cout <<"Failed to read protein function DB file! Exiting..." << endl;
			return -1;
		}
		cout << " done!" << endl;
	}
	
	//calculate mass of peptides or proteins from sequence and amino acid mass databases
	if(par.calcMW)
	{
		if(!utils::fileExists(par.wd + par.staticModsFname))
			if(!par.writeSmod(par.wd))
				cout << "Failed to write new smod file" << endl;
		
		cout << endl << "Calculating protein molecular weights from: " << par.mwDBFname;
		if(!proteins.readInMWdb(par.wd, par))
		{
			cout << "Failed to read mwDB files! Exiting..." << endl;
			return -1;
		}
		cout << " done!" << endl;
	}
	
	if((!par.calcMW && par.getSeq ) || (par.calcMW  && (par.seqDBfname != par.mwDBFname)))
	{
		cout << endl << "Getting protein sequences from " << par.seqDBfname;
		if(!proteins.readInSeqDB(par.seqDBfname))
		{
			cout << "Failed to read seqDB file! Exiting..." << endl;
			return -1;
		}
		cout << " done!" << endl;
	}
	
	if(!par.includeReverse)
		cout << endl << "Skipping reverse matches..." << endl;
	
	//read in and combine files
	cout << endl;
	if(!proteins.readIn(&par, &peptides))
		return -1;
	
	cout << endl << par.numFiles << " files combined." << endl;
	
	if(!par.sampleNamePrefix.empty())
		cout << endl << "Parsing colnames by prefix: " << par.sampleNamePrefix << endl;
	
	//write out combined protein data
	if(par.includeProteins)
	{
		assert(par.outputFormat != params::none);
		if (par.outputFormat == params::wideFormat || par.outputFormat == params::both)
		{
			cout << endl << "Writing protein data...";
			if(!proteins.writeOut(par.wd + par.ofname, par))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return -1;
			}
			cout << " done!" << endl << "Protein data written in wide format to: "
			<< par.ofname << endl << endl;
		}
		if(par.outputFormat == params::longFormat || par.outputFormat == params::both)
		{
			cout << endl << "Writing protein data...";
			if(!proteins.writeOutDB(par.wd + par.dbOfname, par))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return -1;
			}
			cout << " done!" << endl << "Protein data written in long format to: "
			<< par.dbOfname << endl << endl;
		}
	}
	//write out combined peptide data
	if(par.includePeptides)
	{
		assert(par.peptideOutput != params::none);
		if(par.peptideOutput == params::wideFormat || par.peptideOutput == params::both)
		{
			cout << endl << "Writing peptide data...";
			if(!peptides.writeOut(par.wd + par.peptideOfFname, par))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return -1;
			}
			cout << " done!" << endl << "Peptide data written in wide format to: "
			<< par.peptideOfFname << endl << endl;
		}
		
		if(par.peptideOutput == params::longFormat || par.peptideOutput == params::both)
		{
			cout << endl << "Writing peptide data...";
			if(!peptides.writeOutDB(par.wd + par.dbPeptideOfFname, par))
			{
				cout << "Could not write out file! Exiting..." << endl;
				return -1;
			}
			cout << " done!" << endl << "Peptide data written in long format to: "
			<< par.dbPeptideOfFname << endl << endl;
		}
	}
	//write saintExpress input files
	if(par.includeSaint)
	{
		cout << endl << "Writing saint output files...";
		if(!proteins.writeSaint(par.wd + par.saintPreyFname, 1))
		{
			cout << "Could not write prey file! Exiting..." << endl;
			return -1;
		}
		if(!proteins.writeSaint(par.wd + par.saintInteractionFname, 2))
		{
			cout << "Could not write interaction file! Exiting..." << endl;
			return -1;
		}
		cout << " done!" << endl << "prey file written to: " << par.saintPreyFname << endl
		<< "interaction data written to: " << par.saintPreyFname << endl << endl;
 	}
	
	return 0;
}

