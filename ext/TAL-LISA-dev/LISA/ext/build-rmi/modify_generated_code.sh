set -v
#Prefix
P=$1
m=$2
sed -i 's/#include <filesystem>//g' $P.cpp
sed -i 's/std::filesystem::path(dataPath) \/ \"/\"RMI\/rmi_data\//g' $P.cpp
sed -i 's/void cleanup/bool save(char const *filename);\nvoid cleanup/g' $P.h
sed -i "s/void cleanup/bool save(char const *filename)\n{\n    int64_t L1_SIZE = $m;\n    std::ofstream outstream(filename, std::ofstream::binary);\n    outstream.seekp(0);\n    outstream.write((char *)\&(L0_PARAMETER0), sizeof(double));\n    outstream.write((char *)\&(L0_PARAMETER1), sizeof(double));\n    outstream.write((char *)\&(L1_SIZE), sizeof(int64_t));\n    outstream.write(L1_PARAMETERS, L1_SIZE * 24);\n    outstream.close();\n}\n\nvoid cleanup/g" $P.cpp
