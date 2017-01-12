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

#define BIN_VERSION_NUM "1.50"
#define MIN_BIN_VERSION_NUM 1.3

#include "dtafilter.hpp"
#include "../lib/BinTree.cpp"
#include "../lib/hashTable.cpp"
#include "../lib/utils.cpp"
#include "FilterFile.cpp"
#include "dtafilter.cpp"
#include "calcMW.cpp"
#include "dbase.cpp"

#endif /* DTarray_pro.hpp */
