#ifndef CONTEXT_DB_H
#define CONTEXT_DB_H

#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "DBReader.h"

struct ContextFeatureRecord {
    uint32_t scaffoldId;
    uint32_t start;
    uint32_t end;
    DBKeyType targetKey;
    uint64_t featureNameOffset;
    int8_t strand;
};

class ContextDb {
public:
    static size_t featureRecordSize();
    static size_t contextRecordSize();
    static void appendContext(std::string &payload, uint32_t firstFeatureId, uint8_t leftCount, uint8_t rightCount);
    static void readContext(const char *&ptr, const char *end, uint32_t &firstFeatureId,
                            uint8_t &leftCount, uint8_t &rightCount);

    template <typename T>
    static void appendValue(std::string &out, T value) {
        out.append(reinterpret_cast<const char *>(&value), sizeof(T));
    }

    template <typename T>
    static T readValue(const char *&ptr, const char *end) {
        if (ptr > end || static_cast<size_t>(end - ptr) < sizeof(T)) {
            failTruncatedPayload();
        }
        T value;
        memcpy(&value, ptr, sizeof(T));
        ptr += sizeof(T);
        return value;
    }

private:
    static void failTruncatedPayload();
};

class ContextDbWriter {
public:
    ContextDbWriter(const std::string &contextDb);

    uint32_t writeFeature(const std::string &featureId, const std::string &targetSequenceId,
                          DBKeyType targetKey, uint32_t scaffoldId, uint32_t start,
                          uint32_t end, int8_t strand);
    void writeScaffold(const std::string &scaffold);

private:
    std::ofstream featureFile;
    std::ofstream nameFile;
    std::ofstream scaffoldFile;
    uint32_t nextFeatureId;
};

class ContextDbReader {
public:
    ContextDbReader(const std::string &contextDb);

    ContextFeatureRecord readFeature(uint32_t featureId);
    std::pair<std::string, std::string> readNames(const ContextFeatureRecord &record);
    const std::string &scaffoldName(uint32_t scaffoldId) const;

private:
    static std::string readNullTerminatedName(std::ifstream &input);
    static std::vector<std::string> readScaffolds(const std::string &filename);

    std::ifstream featureFile;
    std::ifstream nameFile;
    std::vector<std::string> scaffolds;
};

#endif
