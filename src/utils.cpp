
#include "utils.hpp"

namespace util{
	
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
	void split (const string str, char delim, vector<string> & elems)
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
	
	//gets new line from inF and removes trailing and leading whitespace
	inline void getLineTrim(ifstream& inF, string& line)
	{
		getline(inF, line);
		line = trim(line);
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
	
	template<class T>
	typename vector<T>::iterator insertSorted(vector<T> * const vec, T const& item)
	{
		return vec->insert(upper_bound(vec->begin(), vec->end(), item), item);
	}
	
	template <class T>
	long binSearch(const vector<T>* const vec, const T& findItem, long begin, long end)
	{
		if (begin > end)
			return -1;
		
		long mid = (begin + end)/2;
		
		if (vec->at(mid) == findItem)
			return mid;
		
		if (vec->at(mid) < findItem)
			return binSearch(vec, findItem, mid + 1, end);
		else
			return binSearch(vec, findItem, begin, mid - 1);
	}
}

