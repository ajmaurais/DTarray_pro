//
//  utils.cpp
//  costom general utils library
//
//  Created by Aaron Maurais on 8/31/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#include <utils.hpp>

namespace utils{
	
	/*############# File ###############*/
	
	bool File::read(std::string _fname)
	{
		fname = _fname;
		return read();
	}
	
	bool File::read()
	{
		if(fname.empty())
			throw std::runtime_error("File must be specified!");
		
		std::ifstream inF(fname.c_str());
		
		if(!inF)
			return false;
		
		delimType = detectLineEnding(inF);
		if(delimType == unknown)
			throw std::runtime_error("Invalid delimiter in file: " + fname + "!");
		delim = utils::getDelim(delimType);
		
		inF.seekg(0, inF.end);
		const long size = inF.tellg();
		inF.seekg(0, inF.beg);
		buffer = (char*) calloc(size, sizeof(char));
		
		if(inF.read(buffer, size))
		{
			ss << buffer;
			if(delimType == crlf)
				beginLine = 1;
			else beginLine = 0;
			free(buffer);
			slen = ss.str().length();
			return true;
		} else {
			free(buffer);
			return false;
		}
	}
	
	inline std::string File::getLine()
	{
		std::string ret;
		utils::getLine(ss, ret, delim, beginLine);
		return ret;
	}
	
	inline std::string File::getLine_skip()
	{
		std::string ret;
		do{
			utils::getLine(ss, ret, delim, beginLine);
		} while ((isCommentLine(ret) || ret.empty()) && ss);
		return ret;
	}
	
	inline std::string File::getLine_trim()
	{
		std::string ret;
		utils::getLineTrim(ss, ret, delim, beginLine);
		return ret;
	}
	
	inline std::string File::getLine_skip_trim()
	{
		std::string ret;
		do{
			utils::getLineTrim(ss, ret, delim, beginLine);
		} while ((isCommentLine(ret) || ret.empty()) && ss);
		return ret;
	}
	
	inline std::string File::getLine_trim_skip()
	{
		return getLine_skip_trim();
	}
	
	/*############# functions ###############*/
	
	/*******************/
	/*  file utilities */
	/*******************/
	
	newline_type detectLineEnding_killStream(std::ifstream& inF) {
		char tmp;
		while(inF){
			inF.get(tmp);
			if(tmp == inF.widen('\r')) {	// old Mac or Windows
				inF.get(tmp);
				if(tmp == inF.widen('\n'))	// Windows
					return crlf;
				return cr;	// old Macs
			}
			if(tmp == inF.widen('\n'))	// Unix and modern Macs
				return lf;
		}
		return unknown;
	}
	
	newline_type detectLineEnding(std::ifstream& inF)
	{
		if(!inF)
			throw std::runtime_error("Bad file stream!");
		const std::streampos p = inF.tellg();
		const newline_type ret = detectLineEnding_killStream(inF);
		inF.seekg(p);
		return ret;
	}
	
	char getDelim(newline_type type)
	{
		switch(type){
			case lf : return '\n';
				break;
			case crlf : return '\r';
				break;
			case cr : return '\r';
				break;
			case unknown : return '\n';
				break;
			default : throw std::runtime_error("Invalid type!");
				break;
		}
	}
	
	//returns true if folder at end of path exists and false if it does not
	bool dirExists (const char* path)
	{
		struct stat buffer;
		return stat(path, &buffer) == 0 && S_ISDIR(buffer.st_mode);
	}
	
	//returns true if file at end of path exists and false if it does not
	bool fileExists(const char* path)
	{
		struct stat buffer;
		return stat(path, &buffer) == 0;
	}
	
	bool fileExists(std::string path)
	{
		return fileExists(path.c_str());
	}
	
	bool dirExists(std::string path)
	{
		return dirExists(path.c_str());
	}
	
	//returns dirrectory from which program is run
	std::string pwd()
	{
		char temp[PATH_MAX + 1];
		return (getcwd(temp, PATH_MAX) ? std::string(temp) : std::string(""));
	}
	
	//resolves relative and symbolic file references
	std::string absPath(const char* _fname)
	{
		char fbuff [PATH_MAX + 1];
		realpath(_fname, fbuff);
		return std::string(fbuff);
	}
	
	//resolves relative and symbolic file references
	std::string absPath(std::string _fname)
	{
		return(absPath(_fname.c_str()));
	}
	
	bool ls(const char* path, std::vector<std::string>& files)
	{
		files.clear();
		DIR* dirFile = opendir(path);
		if(!dirFile)
			return false;
		else{
			struct dirent* hFile;
			while((hFile = readdir(dirFile)))
			{
				//skip . and ..
				if(!strcmp(hFile->d_name, ".") || !strcmp(hFile->d_name, ".."))
					continue;
				
				//skip hidden files
				if(IGNORE_HIDDEN_FILES && (hFile->d_name[0] == '.'))
					continue;
				
				//add to files
				files.push_back(hFile->d_name);
			}
			closedir(dirFile);
		}
		return true;
	}
	
	bool ls(const char* path, std::vector<std::string>& files, std::string extension)
	{
		files.clear();
		std::vector<std::string> allFiles;
		if(!ls(path, allFiles))
			return false;
		
		for(std::vector<std::string>::iterator it = allFiles.begin(); it != allFiles.end(); ++it)
			if(endsWith(*it, extension))
				files.push_back(*it);
		return true;
	}
	
	//exicutes std::string arg as bash command
	void systemCommand(std::string command)
	{
		system(command.c_str());
	}
	
	//make dir.
	//returns true if sucessful
	bool mkdir(const char* path)
	{
		//get abs path and make sure that it dosen't already exist
		//return false if it does
		std::string rpath = absPath(path);
		if(dirExists(rpath))
			return false;
		
		//make sure that parent dir exists
		//return false if not
		size_t pos = rpath.find_last_of("/");
		assert(pos != std::string::npos);
		std::string parentDir = rpath.substr(0, pos);
		if(!dirExists(parentDir))
			return false;
		
		//make dir
		systemCommand("mkdir " + rpath);
		
		//test that new dir exists
		return dirExists(rpath);
	}
	
	std::string baseName(const std::string& path, const std::string& delims)
	{
		return path.substr(path.find_last_of(delims) + 1);
	}
	
	std::string removeExtension(const std::string& filename)
	{
		std::string::size_type const p(filename.find_last_of('.'));
		return p > 0 && p != std::string::npos ? filename.substr(0, p) : filename;
	}

	std::string getExtension(const std::string& filename)
	{
		std::string::size_type const p(filename.find_last_of('.'));
		return p > 0 && p != std::string::npos ? filename.substr(p) : filename;
	}
	
	/*********************/
	/*  type conversions */
	/*********************/
	
	//converts num to std::string because to_std::string does not work with stl 98
	template <typename _Tp>
	inline std::string toString(_Tp num)
	{
		std::string str;
		std::stringstream convert;
		
		convert << num;
		convert >> str;
		
		return str;
	}
	
	//converts std::string to int because atoi does not work with stl 98
	//Pre: str must be a std::string with a valid interger conversion
	inline int toInt(std::string str)
	{
		assert(isInteger(str));
		
		int num;
		std::stringstream convert;
		
		convert << str;
		convert >> num;
		
		return num;
	}
	
	//return true if str can be converted to an int
	bool isInteger(std::string str)
	{
		if(str.empty() || ((!isdigit(str[0])) && (str[0] != '-') && (str[0] != '+')))
			return false ;
		
		char * p ;
		strtol(str.c_str(), &p, 10) ;
		
		return (*p == 0) ;
	}
	
	//return true if str can be converted to a double
	bool isDouble(std::string str)
	{
		if(str.empty() || ((!isdigit(str[0])) && (str[0] != '-') && (str[0] != '+')))
			return false ;
		
		char * p ;
		strtod(str.c_str(), &p);
		
		return !(*p != '\0' || p == str);
	}
	
	//converts std::string to int because stod does not work with some c++ compilers
	//Precondition: str must be a std::string with a valid double conversion
	double toDouble(std::string str)
	{
		assert(isDouble(str));
		double num;
		std::stringstream convert;
		
		convert << str;
		convert >> num;
		
		return num;
	}
	
	/*****************/
	/*  std::string utils */
	/*****************/
	
	//returns true if findTxt is found in whithinTxt and false if it it not
	inline bool strContains(std::string findTxt, std::string whithinTxt)
	{
		return whithinTxt.find(findTxt) != std::string::npos;
	}
	
	//overloaded version of strContains, handels findTxt as char
	inline bool strContains(char findTxt, std::string whithinTxt)
	{
		return strContains(std::string(1, findTxt), whithinTxt);
	}
	
	inline bool startsWith(std::string whithinStr, std::string findStr)
	{
		return (whithinStr.find(findStr) == 0);
	}
	
	inline bool endsWith(std::string whithinStr, std::string findStr)
	{
		size_t pos = whithinStr.rfind(findStr);
		if(pos == std::string::npos)
			return false;
		
		return (whithinStr.substr(pos) == findStr);
	}
	
	//split str by delim and populate each split into elems
	inline void split (const std::string& str, const char delim, std::vector<std::string>& elems)
	{
		elems.clear();
		std::stringstream ss (str);
		std::string item;
		
		while(getline(ss, item, delim)) {
			elems.push_back(item);
		}
	}
	
	//remove trailing WHITESPACE
	std::string trimTraling(const std::string& str)
	{
		if(str.empty())
			return "";
		return str.substr(0, str.find_last_not_of(WHITESPACE) + 1);
	}
	
	//remove leading WHITESPACE
	std::string trimLeading(const std::string& str)
	{
		if(str.empty())
			return "";
		return str.substr(str.find_first_not_of(WHITESPACE));
	}
	
	//remove trailing and leading WHITESPACE
	std::string trim(const std::string& str)
	{
		if(str.empty())
			return "";
		return trimLeading(trimTraling(str));
	}
	
	void trimAll(std::vector<std::string>& elems)
	{
		for(std::vector<std::string>::iterator it = elems.begin(); it != elems.end(); ++it)
			*it = utils::trim(*it);
	}
	
	void removeBlanks(std::vector<std::string>& elems)
	{
		for(std::vector<std::string>::iterator it = elems.begin(); it != elems.end();)
		{
			if(it->empty())
				elems.erase(it);
			else ++it;
		}
	}
	
	//returns true if line begins with COMMENT_SYMBOL, ignoring leading whitespace
	bool isCommentLine(std::string line)
	{
		line = trimLeading(line);
		return line.substr(0, COMMENT_SYMBOL.length()) == COMMENT_SYMBOL;
	}
	
	//gets new line from is and removes trailing and leading whitespace
	inline void getLineTrim(std::istream& is, std::string& line, char delim, size_t beginLine)
	{
		utils::getLine(is, line, delim, beginLine);
		line = trim(line);
	}
	
	//gets new line from is and handels non zero start to line
	inline void getLine(std::istream& is, std::string& line, char delim, size_t beginLine)
	{
		getline(is, line, delim);
		if(beginLine > 0)
			line = line.substr(beginLine);
	}
	
	//removes findStr from whithinStr and returns whithinStr
	inline std::string removeSubstr(std::string findStr, std::string whithinStr)
	{
		std::string::size_type i = whithinStr.find(findStr);
		
		if(i != std::string::npos)
			whithinStr.erase(i, findStr.length());
		
		return whithinStr;
	}
	
	inline std::string removeChars(char findChar, std::string wStr)
	{
		wStr.erase(remove(wStr.begin(), wStr.end(), findChar), wStr.end());
		return wStr;
	}
	
	std::string toLower(std::string str)
	{
		transform(str.begin(), str.end(), str.begin(), ::tolower);
		return str;
	}
	
	std::string repeat(std::string str, size_t numTimes)
	{
		std::string ret = "";
		assert(!str.empty());
		for(int i = 0; i < numTimes; i++)
			ret += str;
		return ret;
	}
	
	std::string toSubscript(int _num)
	{
		std::string strNum = utils::toString(_num);
		std::string ret = "";
		
		size_t len = strNum.length();
		for(size_t i = 0; i < len; i++)
		{
			int tempInt = (int)strNum[i] - 48; //convert char to int
			assert(tempInt >= 0 && tempInt <= 9); //check that tempInt will not overrun buffer
			ret += SUBSCRIPT_MAP[tempInt];
		}
		return ret;
	}
	
	/*std::string toSuperscript(int _num)
	{
		std::string strNum = utils::toString(_num);
		std::string ret = "";
		
		size_t len = strNum.length();
		for(size_t i = 0; i < len; i++)
		{
			int tempInt = (int)strNum[i] - 48; //convert char to int
			assert(tempInt >= 0 && tempInt <= 9); //check that tempInt will not overrun buffer
			ret += SUPERSCRIPT_MAP[tempInt];
		}
		return ret;
	}*/
	
	/*********/
	/* other */
	/*********/
	
	bool isFlag(const char* tok)
	{
		if(tok == nullptr)
			return false;
		else return tok[0] == '-';
	}
	
	bool isArg(const char* tok)
	{
		if(tok == nullptr)
			return false;
		else return !isFlag(tok);
	}
	
	std::string ascTime()
	{
		//get current time
		const char* curTime;
		time_t rawTime;
		struct tm * timeInfo;
		time(&rawTime);
		timeInfo = localtime(&rawTime);
		curTime = asctime(timeInfo);
		return(std::string(curTime));
	}
}
