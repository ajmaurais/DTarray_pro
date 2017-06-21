//
//  main.cpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 3/25/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#include "DTarray_pro.hpp"

using namespace std;

int main(int argc, char* argv[])
{	
	//read in names of files to combine and output params
	params::Params par;
	if(!par.getOpts(argc, argv))
		return -1;
	
	cerr << endl << "DTarray_pro v" << BIN_VERSION_NUM << endl;
	
	if(!utils::fileExists(par.getwd() + par.flistName) || par.rewriteFlist)
	{
		assert(par.writeFlist());
		cerr << endl << "Generating " << par.flistName << " using " << par.inputFormat << " input format." << endl;
	}
	if(!par.readFlist(par.flistName, par.getwd()))
	{
		cerr << "\nFailed to read file list! Exiting..." << endl;
		return -1;
	}
	if(par.getNumFiles() <= 0)
	{
		cerr << endl << "No filter files found. Exiting..." << endl << endl;
		return -1;
	}
	if(!par.optionsCompatable())
		return -1;
	
	if(par.peptideOutput != params::none)
		cerr << endl << "Grouping peptides " <<
		params::groupFormatString(par.peptideGroupMethod) << "." << endl;
	
	if(par.modGroupMethod != 0)
		cerr << endl << "Grouping modified peptides." << endl;
	
	Proteins proteins(par);
	Peptides peptides(par);
	
	//read saint bait file
	if(par.includeSaint)
	{
		if(!proteins.readBaitFile(par.saintBaitFile))
		{
			cerr << endl << "Could not read bait file! Exiting..." << endl;
			return -1;
		}
	}
	//read in subcellular locations database
	if(par.getSubCelluarLoc)
	{
		cerr << endl << "Reading subcellular locations database...";
		if(!proteins.readInLocDB(par.locDBfname))
		{
			cerr <<"Failed to read protein location DB file! Exiting..." << endl;
			return -1;
		}
		cerr << " done!" << endl;
	}
	//read in fxn database
	if(par.getFxn)
	{
		cerr << endl << "Reading protein function database...";
		if(!proteins.readInFxnDB(par.fxnDBfname))
		{
			cerr <<"Failed to read protein function DB file! Exiting..." << endl;
			return -1;
		}
		cerr << " done!" << endl;
	}
	//calculate mass of peptides or proteins from sequence and amino acid mass databases
	if(par.calcMW)
	{
		if(!utils::fileExists(par.getwd() + par.staticModsFname))
			if(!par.writeSmod(par.getwd()))
				cerr << "Failed to write new smod file" << endl;
		
		cerr << endl << "Getting protein sequences from " << par.mwDBFname << "...";
		if(!proteins.readInMWdb(par.getwd(), par))
		{
			cerr << "Failed to read mwDB files! Exiting..." << endl;
			return -1;
		}
		cerr << " done!" << endl;
	}
	if((!par.calcMW && par.getSeq ) || (par.calcMW  && (par.seqDBfname != par.mwDBFname)))
	{
		cerr << endl << "Getting protein sequences from " << par.seqDBfname << "...";
		if(!proteins.readInSeqDB(par.seqDBfname))
		{
			cerr << "Failed to read seqDB file! Exiting..." << endl;
			return -1;
		}
		cerr << " done!" << endl;
	}
	if(!par.includeReverse)
		cerr << endl << "Skipping reverse matches..." << endl;
	
	//read in and combine files
	cerr << endl;
	if(!proteins.readIn(&par, &peptides))
		return -1;
	cerr << endl << par.getNumFiles() << " files combined." << endl;
	
	if(!par.sampleNamePrefix.empty())
		cerr << endl << "Parsing colnames by prefix: " << par.sampleNamePrefix << endl;
	
	//write out combined protein data
	if(par.includeProteins)
	{
		assert(par.outputFormat != params::none);
		if (par.outputFormat == params::wideFormat || par.outputFormat == params::both)
		{
			cerr << endl << "Writing protein data...";
			if(!proteins.writeOut(par.getwd() + par.ofname, par))
			{
				cerr << "Could not write out file! Exiting..." << endl;
				return -1;
			}
			cerr << " done!" << endl << "Protein data written in wide format to: "
			<< par.ofname << endl << endl;
		}
		if(par.outputFormat == params::longFormat || par.outputFormat == params::both)
		{
			cerr << endl << "Writing protein data...";
			if(!proteins.writeOutDB(par.getwd() + par.dbOfname, par))
			{
				cerr << "Could not write out file! Exiting..." << endl;
				return -1;
			}
			cerr << " done!" << endl << "Protein data written in long format to: "
			<< par.dbOfname << endl << endl;
		}
	}
	
	//write out combined peptide data
	if(par.includePeptides)
	{
		assert(par.peptideOutput != params::none);
		if(par.peptideOutput == params::wideFormat || par.peptideOutput == params::both)
		{
			cerr << endl << "Writing peptide data...";
			if(!peptides.writeOut(par.getwd() + par.peptideOfFname, par))
			{
				cerr << "Could not write out file! Exiting..." << endl;
				return -1;
			}
			cerr << " done!" << endl << "Peptide data written in wide format to: "
			<< par.peptideOfFname << endl << endl;
		}
		
		if(par.peptideOutput == params::longFormat || par.peptideOutput == params::both)
		{
			cerr << endl << "Writing peptide data...";
			if(!peptides.writeOutDB(par.getwd() + par.dbPeptideOfFname, par))
			{
				cerr << "Could not write out file! Exiting..." << endl;
				return -1;
			}
			cerr << " done!" << endl << "Peptide data written in long format to: "
			<< par.dbPeptideOfFname << endl << endl;
		}
	}
	
	//write saintExpress input files
	if(par.includeSaint)
	{
		cerr << endl << "Writing saint output files...";
		if(!proteins.writeSaint(par.getwd() + par.saintPreyFname, 1))
		{
			cerr << "Could not write prey file! Exiting..." << endl;
			return -1;
		}
		if(!proteins.writeSaint(par.getwd() + par.saintInteractionFname, 2))
		{
			cerr << "Could not write interaction file! Exiting..." << endl;
			return -1;
		}
		cerr << " done!" << endl << "prey file written to: " << par.saintPreyFname << endl
		<< "interaction data written to: " << par.saintPreyFname << endl << endl;
 	}
	
	//write sub celluar localization report
	if(par.locOutput != params::none)
	{
		proteins.buildLocTable();
		
		if(par.locOutput == params::wideFormat || par.locOutput == params::both)
		{
			cerr << endl << "Writing sub celluar loc summary table...";
			if(!proteins.writeWideLocTable(par.getwd() + par.locTableFname, par))
			{
				cerr << "Could not write loc summary table! Exiting..." << endl;
				return -1;
			}
			cerr << " done!" << endl << "Subcellular loc summary table written in wide format to: "
			<< par.locTableFname << endl << endl;
		}
		if(par.locOutput == params::longFormat || par.locOutput == params::both)
		{
			cerr << endl << "Writing sub celluar loc summary table...";
			if(!proteins.writeLongLocTable(par.getwd() + par.locTableLongFname, par))
			{
				cerr << "Could not write loc summary table! Exiting..." << endl;
				return -1;
			}
			cerr << " done!" << endl << "Subcellular loc summary table written in long format to: "
			<< par.locTableLongFname << endl << endl;
		}
	}
	
	return 0;
}

