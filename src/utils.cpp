
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

inline void getLineTrim(ifstream& inF, string& line)
{
	getline(inF, line);
	line = trim(line);
}

inline int strComp(string str1, string str2)
{
	return strcmp(str1.c_str(), str2.c_str());
}
