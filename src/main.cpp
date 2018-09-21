//
//  main.cpp
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

#include <main.hpp>

int main(int argc, char* argv[])
{	
	//read in names of files to combine and output params
	params::Params par;
	if(!par.getOpts(argc, argv))
		return -1;
	
	std::cerr << std::endl << "DTarray_pro v" << BIN_VERSION_NUM << std::endl;
	
	if(!utils::fileExists(par.getwd() + par.flistName) || par.rewriteFlist)
	{
		assert(par.writeFlist());
		std::cerr << std::endl << "Generating " << par.flistName << " using " << par.inputFormat << " input format." << std::endl;
	}
	if(!par.readFlist(par.flistName, par.getwd()))
	{
		std::cerr << "\nFailed to read file list! Exiting..." << std::endl;
		return -1;
	}
	if(par.getNumFiles() <= 0)
	{
		std::cerr << std::endl << "No filter files found. Exiting..." << std::endl << std::endl;
		return -1;
	}
	if(!par.optionsCompatable())
		return -1;
	
	if(par.peptideOutput != params::none)
		std::cerr << std::endl << "Grouping peptides " <<
		params::groupFormatString(par.peptideGroupMethod) << "." << std::endl;
	
	if(par.modGroupMethod != 0)
		std::cerr << std::endl << "Grouping modified peptides." << std::endl;
	
	Proteins proteins(par);
	Peptides peptides(par);
	
	//read saint bait file
	if(par.includeSaint)
	{
		if(!proteins.readBaitFile(par.getwd() + par.saintBaitFile))
		{
			std::cerr << std::endl << "Could not read bait file! Exiting..." << std::endl;
			return -1;
		}
	}
	
	//only do this stuff if peptide output file is to be generated
	if(par.includeProteins)
	{
		//read in subcellular locations database
		if(par.getSubCelluarLoc)
		{
			std::cerr << std::endl << "Reading subcellular locations database...";
			if(!proteins.readInLocDB(par.locDBfname))
			{
				std::cerr <<"Failed to read protein location DB file! Exiting..." << std::endl;
				return -1;
			}
			std::cerr << " done!" << std::endl;
		}
		//read in fxn database
		if(par.getFxn)
		{
			std::cerr << std::endl << "Reading protein function database...";
			if(!proteins.readInFxnDB(par.fxnDBfname))
			{
				std::cerr <<"Failed to read protein function DB file! Exiting..." << std::endl;
				return -1;
			}
			std::cerr << " done!" << std::endl;
		}
		//calculate mass of peptides or proteins from sequence and amino acid mass databases
		if(par.calcMW)
		{
			std::cerr << std::endl << "Getting protein sequences from " << par.getmwDBFname() << "...";
			if(!proteins.readInMWdb(par))
			{
				std::cerr << "Failed to read mwDB files! Exiting..." << std::endl;
				return -1;
			}
			std::cerr << " done!" << std::endl;
			peptides.setMWdb(proteins.get_mwdb());
		}
		if((!par.calcMW && par.getSeq ) || (par.calcMW  && (par.getSeqDBfname() != par.getmwDBFname())))
		{
			std::cerr << std::endl << "Getting protein sequences from " << par.getSeqDBfname() << "...";
			if(!proteins.readInSeqDB(par.getSeqDBfname()))
			{
				std::cerr << "Failed to read seqDB file! Exiting..." << std::endl;
				return -1;
			}
			std::cerr << " done!" << std::endl;
		}
	}
	else if(par.calcMW) //read mwdb for peptides if skipped for proteins
	{
		std::cerr << std::endl << "Getting residue formulas from " << par.atomCountTableFname << "...";
		if(!peptides.readInMWdb(par))
		{
			std::cerr << "Failed to read mwDB files! Exiting..." << std::endl;
			return -1;
		}
		std::cerr << " done!" << std::endl;
	}//end if par.includeProteins
	
	if(!par.includeReverse)
		std::cerr << std::endl << "Skipping reverse matches..." << std::endl;
	
	//read in and combine files
	std::cerr << std::endl;
	if(!proteins.readIn(&par, &peptides))
		return -1;
	std::cerr << std::endl << par.getNumFiles() << " files combined." << std::endl;
	
	if(!par.sampleNamePrefix.empty())
		std::cerr << std::endl << "Parsing colnames by prefix: " << par.sampleNamePrefix << std::endl;
	
	//write out combined protein data
	if(par.includeProteins)
	{
		assert(par.outputFormat != params::none);
		if (par.outputFormat == params::wideFormat || par.outputFormat == params::both)
		{
			std::cerr << std::endl << "Writing protein data...";
			if(!proteins.writeOut(par.getwd() + par.ofname, par))
			{
				std::cerr << "Could not write out file! Exiting..." << std::endl;
				return -1;
			}
			std::cerr << " done!" << std::endl << "Protein data written in wide format to: "
			<< par.ofname << std::endl << std::endl;
		}
		if(par.outputFormat == params::longFormat || par.outputFormat == params::both)
		{
			std::cerr << std::endl << "Writing protein data...";
			if(!proteins.writeOutDB(par.getwd() + par.dbOfname, par))
			{
				std::cerr << "Could not write out file! Exiting..." << std::endl;
				return -1;
			}
			std::cerr << " done!" << std::endl << "Protein data written in long format to: "
			<< par.dbOfname << std::endl << std::endl;
		}
	}
	
	//write out combined peptide data
	if(par.includePeptides)
	{
		assert(par.peptideOutput != params::none);
		if(par.peptideOutput == params::wideFormat || par.peptideOutput == params::both)
		{
			std::cerr << std::endl << "Writing peptide data...";
			if(!peptides.writeOut(par.getwd() + par.peptideOfFname, par))
			{
				std::cerr << "Could not write out file! Exiting..." << std::endl;
				return -1;
			}
			std::cerr << " done!" << std::endl << "Peptide data written in wide format to: "
			<< par.peptideOfFname << std::endl << std::endl;
		}
		
		if(par.peptideOutput == params::longFormat || par.peptideOutput == params::both)
		{
			std::cerr << std::endl << "Writing peptide data...";
			if(!peptides.writeOutDB(par.getwd() + par.dbPeptideOfFname, par))
			{
				std::cerr << "Could not write out file! Exiting..." << std::endl;
				return -1;
			}
			std::cerr << " done!" << std::endl << "Peptide data written in long format to: "
			<< par.dbPeptideOfFname << std::endl << std::endl;
		}
	}
	
	//write saintExpress input files
	if(par.includeSaint)
	{
		std::cerr << std::endl << "Writing saint output files...";
		if(!proteins.writeSaint(par.getwd() + par.saintPreyFname, Proteins::preyFile))
		{
			std::cerr << "Could not write prey file! Exiting..." << std::endl;
			return -1;
		}
		if(!proteins.writeSaint(par.getwd() + par.saintInteractionFname, Proteins::interactionFile))
		{
			std::cerr << "Could not write interaction file! Exiting..." << std::endl;
			return -1;
		}
		std::cerr << " done!" << std::endl << "prey file written to: " << par.saintPreyFname << std::endl
		<< "interaction data written to: " << par.saintPreyFname << std::endl << std::endl;
 	}
	
	//write sub celluar localization report
	if(par.locOutput != params::none)
	{
		proteins.buildLocTable();
		
		if(par.locOutput == params::wideFormat || par.locOutput == params::both)
		{
			std::cerr << std::endl << "Writing sub celluar loc summary table...";
			if(!proteins.writeWideLocTable(par.getwd() + par.locTableFname, par))
			{
				std::cerr << "Could not write loc summary table! Exiting..." << std::endl;
				return -1;
			}
			std::cerr << " done!" << std::endl << "Subcellular loc summary table written in wide format to: "
			<< par.locTableFname << std::endl << std::endl;
		}
		if(par.locOutput == params::longFormat || par.locOutput == params::both)
		{
			std::cerr << std::endl << "Writing sub celluar loc summary table...";
			if(!proteins.writeLongLocTable(par.getwd() + par.locTableLongFname, par))
			{
				std::cerr << "Could not write loc summary table! Exiting..." << std::endl;
				return -1;
			}
			std::cerr << " done!" << std::endl << "Subcellular loc summary table written in long format to: "
			<< par.locTableLongFname << std::endl << std::endl;
		}
	}
	
	return 0;
}

