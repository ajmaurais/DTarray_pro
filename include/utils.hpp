//
//  utils.hpp
//  costom general utils library
//
//  Created by Aaron Maurais on 8/31/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#pragma once

#include <cassert>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <stdexcept>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <ctime>
#include <vector>
#include <fstream>
#include <algorithm>

#ifndef PATH_MAX
	#define PATH_MAX 1024
#endif

namespace utils{
	
	/******************************/
	/* namespace scoped constants */
	/*****************************/
	
	std::string const WHITESPACE = " \f\n\r\t\v";
	std::string const COMMENT_SYMBOL = "#";
	enum newline_type {lf, crlf, cr, unknown};
	char const DEFAULT_LINE_DELIM = '\n';
	size_t const DEFAULT_BEGIN_LINE = 0;
	bool const IGNORE_HIDDEN_FILES = true; //ignore hidden files in utils::ls
	std::string const SUBSCRIPT_MAP [] = {"\u2080", "\u2081", "\u2082", "\u2083",
		"\u2084", "\u2085", "\u2086", "\u2087", "\u2088", "\u2089"};
	/*std::string const SUPERSCRIPT_MAP [] = {"\u00B0", "\u00B1", "\u00B2", "\u00B3",
		"\u00B4", "\u00B5", "\u00B6", "\u00B7", "\u00B8", "\u00B9"};*/
	
	/******************************/
	/*     class definitions     */
	/*****************************/
	
	class File;
	
	//file class for reading in text files line by line
	//automatically detects and handles line return characters from different operating systems
	class File{
	private:
		char* buffer;
		char delim;
		newline_type delimType;
		std::string fname;
		size_t beginLine;
		std::stringstream ss;
		unsigned long slen;
	public:
		File(){
			buffer = nullptr;
			slen = 0;
		}
		File(std::string str){
			buffer = nullptr;
			fname = str;
			slen = 0;
		}
		~File(){}
		
		//modifers
		bool read(std::string);
		bool read();
		std::string getLine();
		std::string getLine_skip();
		std::string getLine_trim();
		std::string getLine_skip_trim();
		std::string getLine_trim_skip();
		
		//properties
		bool end(){
			return (ss.tellg() >= slen);
		}
		std::string getFname() const{
			return fname;
		}
		char getDelim() const{
			return delim;
		}
		newline_type getNewLineType() const{
			return delimType;
		}
	};
	
	/*************/
	/* functions */
	/*************/
	
	//file utils
	char getDelim(newline_type);
	newline_type detectLineEnding_killStream(std::ifstream&);
	newline_type detectLineEnding(std::ifstream&);
	bool dirExists(const char*);
	bool dirExists(std::string);
	bool fileExists(const char*);
	bool fileExists(std::string);
	std::string pwd();
	std::string absPath(std::string);
	std::string absPath(const char*);
	bool isAbsPath(const std::string&);
	bool isAbsPath(const char*);
	bool ls(const char*, std::vector<std::string>&);
	bool ls(const char*, std::vector<std::string>&, std::string);
	bool mkdir(std::string);
	bool mkdir(const char*);
	void systemCommand(std::string command);
	std::string baseName(const std::string& path, const std::string& delims = "/\\");
	std::string removeExtension(const std::string&);
	std::string getExtension(const std::string&);
	
	//type conversions
	template <typename _Tp> std::string toString(_Tp num);
	int toInt(std::string);
	double toDouble(std::string);
	bool isInteger(std::string);
	bool isDouble(std::string);
	
	//std::string utils
	bool strContains(std::string, std::string);
	bool strContains(char, std::string);
	bool startsWith(std::string whithinStr, std::string findStr);
	bool endsWith(std::string whithinStr, std::string findStr);
	void split (const std::string&, const char, std::vector<std::string>&);
	std::string trimTraling(const std::string&);
	std::string trimLeading(const std::string&);
	std::string trim(const std::string&);
	void trimAll(std::vector<std::string>&);
	void removeBlanks(std::vector<std::string>&);
	bool isCommentLine(std::string);
	std::string removeSubstr(std::string, std::string);
	std::string removeChars(char, std::string);
	std::string toLower(std::string);
	std::string repeat(std::string, size_t);
	void getLineTrim(std::istream& is, std::string& line,
		char delim = DEFAULT_LINE_DELIM, size_t beg= DEFAULT_BEGIN_LINE);
	void getLine(std::istream& is, std::string& line,
		char delim = DEFAULT_LINE_DELIM, size_t beg= DEFAULT_BEGIN_LINE);
	std::string toSubscript(int);
	//std::string toSuperscript(int);
	
	//other
	bool isFlag(const char*);
	bool isArg(const char*);
	std::string ascTime();
	template <typename _Tp, size_t N> _Tp* begin(_Tp(&arr)[N]) {
		return &arr[0];
	}
	template <typename _Tp, size_t N> _Tp* end(_Tp(&arr)[N]) {
		return &arr[0]+N;
	}
}

/* utils_hpp */
