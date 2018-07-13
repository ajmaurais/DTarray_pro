//
//  dtarray_pro.hpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 1/2/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#pragma once

#include <stdio.h>
#include <string>

//deal with older c++ compilers
#if (__cplusplus == 199711L || __cplusplus == 1)
	#define nullptr NULL
#endif

//make sure this value is defined
#ifndef PATH_MAX
	#define PATH_MAX 1024
#endif

#define BIN_VERSION_NUM "1.7"
#define MIN_BIN_VERSION_NUM 1.60

#ifndef OUT_DELIM
#define OUT_DELIM '\t'
#endif

#ifndef IN_DELIM
#define IN_DELIM '\t'
#endif

std::string const REVERSE_MATCH = "Reverse_";

std::string parseReplicate(std::string);
std::string parseSample(std::string, std::string, bool, bool);

/* dtarray_pro.hpp */
