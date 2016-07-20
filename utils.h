
#include <sys/stat.h>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

string const WHITESPACE = " \f\n\r\t\v";

//returns true if folder at end of path exists and false if it does not
bool dirExists (string path)
{
	struct stat buffer;
	if (stat(path.c_str(), &buffer) == 0 && S_ISDIR(buffer.st_mode))
		return true;
	return false;
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
int toInt(string str)
{
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
	if (str == "")
		return "";
	return str.substr( 0, str.find_last_not_of(WHITESPACE) + 1 );
}

//remove leading WHITESPACE
inline string trimLeading(const string& str)
{
	if (str == "")
		return "";
	return str.substr(str.find_first_not_of(WHITESPACE));
}

//remove trailing and leading WHITESPACE
inline string trim(const string& str)
{
	if (str == "")
		return "";
	return trimLeading( trimTraling(str));
}

bool isCommentLine(string line)
{
	line = trimLeading(line);
	if(line[0] == '#')
		return true;
	return false;
}

