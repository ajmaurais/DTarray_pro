
#include <sys/stat.h>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <cassert>

using namespace std;

string const WHITESPACE = " \f\n\r\t\v";
string const COMMENT_SYMBOL = "#"; //if changed value must also be changed in DTarray_AJM.sh

bool dirExists (string);
bool fileExists (string);
string toString(int);
int toInt(string);
bool strContains(char, string);
void split (const string, char, vector<string> &);
inline string trimTraling(const string&);
inline string trimLeading(const string&);
inline string trim(const string&);
bool isCommentLine(string);
bool isInteger(string);

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
//Precondition: str must be an integer stored as a string
int toInt(string str)
{
	assert(isInteger(str));
	
	int num;
	stringstream convert;
	
	convert << str;
	convert >> num;
	
	return num;
}

//search through string for char and return true if char is found
bool strContains(char findTxt, string whithinTxt)
{
	int len = int(whithinTxt.length());
	
	for(int i = 0; i < len; i++)
		if(whithinTxt[i] == findTxt)
			return true;
	
	return false;
}

//split str by delim and populate each split into elems
void split (const string str, char delim, vector<string> & elems)
{
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

