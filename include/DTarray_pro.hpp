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

#define BIN_VERSION_NUM "1.7"
#define MIN_BIN_VERSION_NUM 1.60

#include <gitVersion.hpp>
#include <utils.hpp>
#include <params.hpp>
#include <dtafilter.hpp>
#include <calcMW.hpp>
#include <dbase.hpp>
#include <saintOutput.hpp>
#include <locReport.hpp>
#include <molecularFormula.hpp>
#include <dtafilter.hpp>

#endif /* DTarray_pro.hpp */
