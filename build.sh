
#recompile DTarray_pro

GIT_VERSION_HPP="src/gitVersion.hpp"

#if dir is under git versioning, get updated git info
if [ -d .git ] ; then
	echo "Getting git version number..."
	gitDate=$(git log -1 --format=%cd --date=local)
	gitCommit=$(git rev-parse HEAD)
	echo -e "\n#ifndef gitVersion_hpp" > $GIT_VERSION_HPP
	echo -e "#define gitVersion_hpp\n" >> $GIT_VERSION_HPP
	echo "const char* GIT_DATE = \"$gitDate\";" >> $GIT_VERSION_HPP
	echo "const char* GIT_COMMIT = \"$gitCommit\";" >> $GIT_VERSION_HPP
	echo -e "\n#endif /* gitVersion_hpp */" >> $GIT_VERSION_HPP
fi

echo "Recompiling DTarray_pro..."

g++ -o bin/DTarray_pro src/main.cpp

echo "Done!"
