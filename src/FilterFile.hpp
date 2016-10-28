//
//  FilterFile.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef FilterFile_hpp
#define FilterFile_hpp

#include <fstream>
#include <cassert>

using namespace std;

/**********************/
/* class definitions */
/*********************/

class FilterFileParams;
class FilterFileParam;

//stores the data pertaining to a specific filter file (or MS run) for each protein
struct FilterFile{
	string colname, count;
	string uniquePeptides;
	
	//constructor
	FilterFile (string, string, string);
};

struct FilterFileParam{
	string path;
	string colname;
};

struct Param {
	string param;
	string value;
	
	//constructor
	Param (string);
};

//stores names and locations for DTA filter files and output paramaters found in
//params file.
class FilterFileParams{
	friend class Proteins;
	vector<FilterFileParam> file;
public:
	int numFiles;
	string outputFormat;
	string sampleNamePrefix;
	bool includeUnique;
	bool getSubCelluarLoc;
	string locDBfname;
	bool calcMW;
	string aaDBfanme, mwDBFname, staticModsFname;
	string ofname;
	bool includeSeq;
	string seqDBfname;
	
	//modifiers
	bool readDTParams(string, string);
	bool readFlist(string, string);
	
	//properties
	inline string getFileColname(int) const;
};


#endif /* FilterFile_hpp */
