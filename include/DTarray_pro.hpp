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

#include <gitVersion.hpp>
#include <dtafilter.hpp>
#include "../src/hashTable.cpp"
#include "../src/utils.cpp"
#include "../src/params.cpp"
#include "../src/dtafilter.cpp"
#include "../src/calcMW.cpp"
#include "../src/dbase.cpp"
#include "../src/saintOutput.cpp"
#include "../src/locReport.cpp"

#endif /* DTarray_pro.hpp */
