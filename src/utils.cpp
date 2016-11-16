
#include "utils.hpp"

namespace utils{
	
	/*############# File ###############*/
	
	bool File::read(string _fname)
	{
		fname = _fname;
		ifstream inF(fname.c_str());
		
		if(!inF)
			return false;
		
		delimType = detectLineEnding(inF);
		if(delimType == unknown)
			throw runtime_error("Invalid delimiter in file: " + fname + "!");
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
			return true;
		}
		else {
			free(buffer);
			return false;
		}
	}
	
	inline string File::getLine()
	{
		string ret;
		utils::getLine(ss, ret, delim, beginLine);
		return ret;
	}
	
	inline string File::getLine_skip()
	{
		string ret;
		do{
			utils::getLine(ss, ret, delim, beginLine);
		} while ((isCommentLine(ret) || ret.empty()) && ss);
		return ret;
	}
	
	inline string File::getLine_trim()
	{
		string ret;
		utils::getLineTrim(ss, ret, delim, beginLine);
		return ret;
	}
	
	inline string File::getLine_skip_trim()
	{
		string ret;
		do{
			utils::getLineTrim(ss, ret, delim, beginLine);
		} while ((isCommentLine(ret) || ret.empty()) && ss);
		return ret;
	}
	
	inline string File::getLine_trim_skip()
	{
		return getLine_skip_trim();
	}
	
	/*############# functions ###############*/
	
	inline newline_type detectLineEnding_killStream(ifstream& inF) {
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
	
	inline newline_type detectLineEnding(ifstream& inF)
	{
		if(!inF)
			throw runtime_error("Bad file stream!");
		const streampos p = inF.tellg();
		const newline_type ret = detectLineEnding_killStream(inF);
		inF.seekg(p);
		return ret;
	}
	
	inline char getDelim(newline_type type)
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
			default : throw runtime_error("Invalid type!");
				break;
		}
	}
	
	//returns true if folder at end of path exists and false if it does not
	bool dirExists (string path)
	{
		struct stat buffer;
		return stat(path.c_str(), &buffer) == 0 && S_ISDIR(buffer.st_mode);
	}
	
	//returns true if file at end of path exists and false if it does not
	bool fileExists (string path)
	{
		struct stat buffer;
		return stat(path.c_str(), &buffer) == 0;
	}
	
	//converts int to string because to_string does not work with some c++ compilers
	string toString(int num)
	{
		string str;
		stringstream convert;
		
		convert << num;
		convert >> str;
		
		return str;
	}
	
	//converts string to int because atoi does not work with some c++ compilers
	//Precondition: str must be a string with a valid interger conversion
	int toInt(string str)
	{
		assert(isInteger(str));
		
		int num;
		stringstream convert;
		
		convert << str;
		convert >> num;
		
		return num;
	}
	
	//returns true if findTxt is found in whithinTxt and false if it it not
	inline bool strContains(string findTxt, string whithinTxt)
	{
		return whithinTxt.find(findTxt) != string::npos;
	}
	
	//overloaded version of strContains handels findTxt as char
	inline bool strContains(char findTxt, string whithinTxt)
	{
		return strContains(string(1, findTxt), whithinTxt);
	}
	
	//split str by delim and populate each split into elems
	inline void split (const string str, char delim, vector<string> & elems)
	{
		elems.clear();
		stringstream ss (str);
		string item;
		
		while (getline(ss, item, delim)) {
			elems.push_back(item);
		}
	}
	
	//remove trailing WHITESPACE
	inline string trimTraling(const string& str)
	{
		if (str.empty())
			return "";
		return str.substr( 0, str.find_last_not_of(WHITESPACE) + 1 );
	}
	
	//remove leading WHITESPACE
	inline string trimLeading(const string& str)
	{
		if (str.empty())
			return "";
		return str.substr(str.find_first_not_of(WHITESPACE));
	}
	
	//remove trailing and leading WHITESPACE
	inline string trim(const string& str)
	{
		if (str.empty())
			return "";
		return trimLeading(trimTraling(str));
	}
	
	//returns true if line begins with COMMENT_SYMBOL, ignoring leading whitespace
	bool isCommentLine(string line)
	{
		line = trimLeading(line);
		return line.substr(0, COMMENT_SYMBOL.length()) == COMMENT_SYMBOL;
	}
	
	//return true if str can be converted to an int
	bool isInteger(string str)
	{
		if(str.empty() || ((!isdigit(str[0])) && (str[0] != '-') && (str[0] != '+')))
			return false ;
		
		char * p ;
		strtol(str.c_str(), &p, 10) ;
		
		return (*p == 0) ;
	}
	
	bool isDouble(string str)
	{
		if(str.empty() || ((!isdigit(str[0])) && (str[0] != '-') && (str[0] != '+')))
			return false ;
		
		char * p ;
		strtod(str.c_str(), &p);
		
		return (*p != '\0' || p == str) ? false : true;
	}
	
	//converts string to int because stod does not work with some c++ compilers
	//Precondition: str must be a string with a valid double conversion
	double toDouble(string str)
	{
		assert(isDouble(str));
		double num;
		stringstream convert;
		
		convert << str;
		convert >> num;
		
		return num;
	}
	
	//gets new line from is and removes trailing and leading whitespace
	inline void getLineTrim(istream& is, string& line, char delim, size_t beginLine)
	{
		utils::getLine(is, line, delim, beginLine);
		line = trim(line);
	}
	
	//gets new line from is and handels non zero start to line
	inline void getLine(istream& is, string& line, char delim, size_t beginLine)
	{
		getline(is, line, delim);
		if(beginLine > 0)
			line = line.substr(beginLine);
	}
	
	//removes findStr from whithinStr and returns whithinStr
	string removeSubstr(string findStr, string whithinStr)
	{
		string::size_type i = whithinStr.find(findStr);
		
		if(i != string::npos)
			whithinStr.erase(i, findStr.length());
		
		return whithinStr;
	}
	
	string toLower(string str)
	{
		transform(str.begin(), str.end(), str.begin(), ::tolower);
		return str;
	}
	
	template<class _Tp>
	typename vector<_Tp>::iterator insertSorted(const vector<_Tp>& vec, const _Tp& item)
	{
		return vec.insert(upper_bound(vec.begin(), vec.end(), item), item);
	}
	
	/**
	 Template binary search. 
	 Pre: vec must be sorted, _Tp must have == , < , and > operator members.
	 Post: returns iterator to positon at which findItem is found. If findItem is 
	 not found, returns vec.end().
	 */
	template <class _Tp>
	typename vector<_Tp>::iterator binSearch(const vector<_Tp>& vec, const _Tp& findItem, long begin, long end)
	{
		if (begin > end)
			return vec.end();
		
		long mid = (begin + end)/2;
		
		if (vec[mid] == findItem)
			return vector<_Tp>::iterator (vec.begin() + mid);
		
		if (vec[mid] < findItem)
			return binSearch(vec, findItem, mid + 1, end);
		else
			return binSearch(vec, findItem, begin, mid - 1);
	}
}

