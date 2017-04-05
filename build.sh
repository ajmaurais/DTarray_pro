
#recompile DTarray_pro

#if dir is under git versioning, get updated git info
if [ -d .git ] ; then
	echo "Getting git version number..."
	gitDate=$(git log -1 --format=%cd --date=local)
	gitCommit=$(git rev-parse HEAD)
	echo -e "\n#ifndef gitVersion_h" > src/gitVersion.hpp
	echo -e "#define gitVersion_h\n" >> src/gitVersion.hpp
	echo "const char* GIT_DATE = \"$gitDate\";" >> src/gitVersion.hpp
	echo "const char* GIT_COMMIT = \"$gitCommit\";" >> src/gitVersion.hpp
	echo -e "\n#endif /* gitVersion_h */" >> src/gitVersion.hpp
fi

echo "Recompiling DTarray_pro..."

g++ -o bin/DTarray_pro src/main.cpp

echo "Done!"
