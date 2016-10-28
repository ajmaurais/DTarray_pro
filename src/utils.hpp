//
//  utils.hpp
//  DTarray_AJM
//
//  Created by Aaron Maurais on 10/28/16.
//  Copyright © 2016 Aaron Maurais. All rights reserved.
//

#ifndef utils_hpp
#define utils_hpp

#include <cassert>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <cstring>
#include <stdlib.h>

namespace util {
	
	using namespace std;
	
	/******************************/
	/* globally scoped constants */
	/*****************************/
	
	string const WHITESPACE = " \f\n\r\t\v";
	string const COMMENT_SYMBOL = "#"; //if changed, paramsCommentSymbol must also be changed in DTarray_AJM.sh
	
	
	/*************/
	/* functions */
	/*************/
	
	bool dirExists (string);
	bool fileExists (string);
	string toString(int);
	int toInt(string);
	inline bool strContains(string, string);
	inline bool strContains(char, string);
	void split (const string, char, vector<string> &);
	inline string trimTraling(const string&);
	inline string trimLeading(const string&);
	inline string trim(const string&);
	bool isCommentLine(string);
	bool isInteger(string);
	inline void getLineTrim(ifstream&, string&);
	string removeSubstr(string, string);
	string toLower(string);
	template<class T> long binSearch(const vector<T>* const, const T&, long, long);
	template<class T> typename vector<T>::iterator insertSorted(vector<T>* const vec, const T& item);
}

#endif /* utils_hpp */
