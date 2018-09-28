
#include <Rcpp.h>
#include <regex>
#include <iostream>
#include <exception>

namespace utils{
  //returns true if findTxt is found in whithinTxt and false if it it not
  bool strContains(std::string findTxt, std::string whithinTxt)
  {
    return whithinTxt.find(findTxt) != std::string::npos;
  }
  
  //removes findStr from whithinStr and returns whithinStr
  std::string removeSubstr(std::string findStr, std::string whithinStr)
  {
    std::string::size_type i = whithinStr.find(findStr);
    
    if(i != std::string::npos)
      whithinStr.erase(i, findStr.length());
    
    return whithinStr;
  }
}


std::string parseSample(std::string sampleName, std::string prefix, bool parseSampleName, bool outputFormat)
{
  //return unparsed sampleName if prefix is empty std::string or is not found in sampleName
  if(!parseSampleName && (!utils::strContains(prefix, sampleName) || prefix.length() == 0)){
    return sampleName;
  }
  else if(parseSampleName && prefix.empty()){
    return sampleName.substr(0, sampleName.find_last_of("_"));
  }
  else {
    std::string sample = utils::removeSubstr(prefix, sampleName); //remove prefix from sampleName
    return outputFormat ? sample.substr(0, sample.find_last_of("_")) : sample;
  }
}

std::string parseSample_new(std::string sampleName, std::string prefix,
                        bool parseSampleName, bool outputFormat, bool re)
{
  std::regex pattern = std::regex(prefix);
  
  //return unparsed sampleName if prefix is empty string or is not found in sampleName
  if(!parseSampleName && (!(re ? std::regex_search(sampleName, pattern) :
                              utils::strContains(prefix, sampleName)) ||
                                prefix.length() == 0)){
    return sampleName;
  }
  else if(parseSampleName && prefix.empty()){
    return sampleName.substr(0, sampleName.find_last_of("_"));
  }
  else {
    std::string sample = re ? std::regex_replace(sampleName, pattern, "") :
    utils::removeSubstr(prefix, sampleName); //remove prefix from sampleName
    return outputFormat ? sample.substr(0, sample.find_last_of("_")) : sample;
  }
}

// [[Rcpp::export]]
Rcpp::StringVector runTest(char func, Rcpp::CharacterVector sampleNames, std::string prefix, bool parseSampleName, bool outputFormat, bool re)
{
  Rcpp::StringVector ret;
  size_t len = sampleNames.size();
  
  for(int i = 0; i < len; i++)
  {
    switch (func) {
    case 'o':
      ret.push_back(parseSample(std::string(sampleNames[i]), prefix, parseSampleName, outputFormat).c_str());
      break;
    case 'n':
      ret.push_back(parseSample_new(std::string(sampleNames[i]), prefix, parseSampleName, outputFormat, re).c_str());
      break;
    default:
      throw std::runtime_error("Invalid arg for function!");
      break;
    }
  }
  
  return ret;
}

