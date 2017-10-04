//
//  DTarray_pro.hpp
//  DTarray_pro
//
//  Created by Aaron Maurais on 1/2/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#ifndef DTarray_pro_hpp
#define DTarray_pro_hpp

//deal with older c++ compilers
#if (__cplusplus == 199711L || __cplusplus == 1)
	#define nullptr NULL
#endif

//make sure this value is defined
#ifndef PATH_MAX
	#define PATH_MAX 1024
#endif

#define BIN_VERSION_NUM "1.63"
#define MIN_BIN_VERSION_NUM 1.60

#include "gitVersion.hpp"
#include "dtafilter.hpp"
#include "../lib/hashTable.cpp"
#include "../lib/utils.cpp"
#include "params.cpp"
#include "dtafilter.cpp"
#include "calcMW.cpp"
#include "dbase.cpp"
#include "saintOutput.cpp"
#include "locReport.cpp"

#endif /* DTarray_pro.hpp */
