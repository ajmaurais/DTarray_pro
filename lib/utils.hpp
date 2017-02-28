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
#include <stdexcept>

namespace utils {
	
	using namespace std;
	
	/******************************/
	/* namespace scoped constants */
	/*****************************/
	
	string const WHITESPACE = " \f\n\r\t\v";
	string const COMMENT_SYMBOL = "#"; //if changed, paramsCommentSymbol must also be changed in DTarray_AJM.sh
	enum newline_type {lf, crlf, cr, unknown};
	char const DEFAULT_LINE_DELIM = '\n';
	size_t const DEFAULT_BEGIN_LINE = 0;
	
	/******************************/
	/*     class definitions     */
	/*****************************/
	
	class File;
	
	class File{
	private:
		char* buffer;
		char delim;
		newline_type delimType;
		string fname;
		size_t beginLine;
		stringstream ss;
	public:
		File(){
			buffer = nullptr;
		}
		File(string str){
			buffer = nullptr;
			fname = str;
		}
		~File(){}
		
		//modifers
		bool read(string);
		bool read();
		inline string getLine();
		inline string getLine_skip();
		inline string getLine_trim();
		inline string getLine_skip_trim();
		inline string getLine_trim_skip();
		
		//properties
		inline bool end() const{
			return (ss.rdbuf()->in_avail() == 0);
		}
		inline string getFname() const{
			return fname;
		}
		inline char getDelim() const{
			return delim;
		}
		inline newline_type getNewLineType() const{
			return delimType;
		}
	};

	
	/*************/
	/* functions */
	/*************/
	
	inline char getDelim(newline_type);
	inline newline_type detectLineEnding_killStream(ifstream&);
	inline newline_type detectLineEnding(ifstream&);
	bool dirExists (string);
	bool fileExists (string);
	inline string toString(int);
	inline string toString(size_t);
	inline int toInt(string);
	inline double toDouble(string);
	inline bool strContains(string, string);
	inline bool strContains(char, string);
	inline void split (const string, const char, vector<string> &);
	inline string trimTraling(const string&);
	inline string trimLeading(const string&);
	inline string trim(const string&);
	bool isCommentLine(string);
	bool isInteger(string);
	inline void getLineTrim(istream& is, string& line, char delim = DEFAULT_LINE_DELIM, size_t beginLine = DEFAULT_BEGIN_LINE);
	inline void getLine(istream& is, string& line, char delim = DEFAULT_LINE_DELIM, size_t beginLine = DEFAULT_BEGIN_LINE);
	inline string removeSubstr(string, string);
	string toLower(string);
	string repeat(string, size_t);
}

#endif /* utils_hpp */
