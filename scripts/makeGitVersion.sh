
#update gitVersion.hpp if dir is under version control.

GIT_VERSION_HPP="include/gitVersion.hpp"

#if dir is under git versioning, get updated git info
if [ -d .git ] ; then
	echo "Updateing gitVersion.hpp..."
	gitDate=$(git log -1 --format=%cd --date=local)
	gitCommit=$(git rev-parse HEAD)
	echo -e "\n#ifndef gitVersion_hpp" > $GIT_VERSION_HPP
	echo -e "#define gitVersion_hpp\n" >> $GIT_VERSION_HPP
	echo "const char* GIT_DATE = \"$gitDate\";" >> $GIT_VERSION_HPP
	echo "const char* GIT_COMMIT = \"$gitCommit\";" >> $GIT_VERSION_HPP
	echo -e "\n#endif /* gitVersion_hpp */" >> $GIT_VERSION_HPP
fi


