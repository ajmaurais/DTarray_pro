//
//  utils.hpp
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
	
	//std::string utils
	bool strContains(std::string findTxt, std::string whithinTxt);
	bool strContains(char findTxt, std::string whithinTxt);
	bool startsWith(std::string whithinStr, std::string findStr);
	bool endsWith(std::string whithinStr, std::string findStr);
	void split(const std::string&, const char, std::vector<std::string>&);
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
