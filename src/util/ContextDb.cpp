#include "ContextDb.h"

#include <cstdlib>
#include <limits>

#include "Debug.h"
#include "Util.h"

namespace {
void validateContext(uint32_t firstFeatureId, uint8_t leftCount, uint8_t rightCount) {
    const uint32_t featureCount = static_cast<uint32_t>(leftCount) + rightCount + 1;
    if (firstFeatureId > std::numeric_limits<uint32_t>::max() - (featureCount - 1)) {
        Debug(Debug::ERROR) << "Invalid context feature range\n";
        EXIT(EXIT_FAILURE);
    }
}
}

size_t ContextDb::featureRecordSize() {
    return sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint32_t) + sizeof(DBKeyType) + sizeof(uint64_t) + sizeof(int8_t);
}

size_t ContextDb::contextRecordSize() {
    return sizeof(uint32_t) + sizeof(uint8_t) + sizeof(uint8_t);
}

void ContextDb::appendContext(std::string &payload, uint32_t firstFeatureId, uint8_t leftCount, uint8_t rightCount) {
    validateContext(firstFeatureId, leftCount, rightCount);
    appendValue(payload, firstFeatureId);
    appendValue(payload, leftCount);
    appendValue(payload, rightCount);
}

void ContextDb::readContext(const char *&ptr, const char *end, uint32_t &firstFeatureId,
                            uint8_t &leftCount, uint8_t &rightCount) {
    firstFeatureId = readValue<uint32_t>(ptr, end);
    leftCount = readValue<uint8_t>(ptr, end);
    rightCount = readValue<uint8_t>(ptr, end);
    validateContext(firstFeatureId, leftCount, rightCount);
}

void ContextDb::failTruncatedPayload() {
    Debug(Debug::ERROR) << "Truncated context payload\n";
    EXIT(EXIT_FAILURE);
}

ContextDbWriter::ContextDbWriter(const std::string &contextDb) :
    featureFile((contextDb + "_features").c_str(), std::ios::binary),
    nameFile((contextDb + "_feature_names").c_str(), std::ios::binary),
    scaffoldFile((contextDb + "_scaffolds").c_str(), std::ios::binary),
    nextFeatureId(0) {
    if (featureFile.fail()) {
        Debug(Debug::ERROR) << "Could not open " << contextDb << "_features for writing\n";
        EXIT(EXIT_FAILURE);
    }
    if (nameFile.fail()) {
        Debug(Debug::ERROR) << "Could not open " << contextDb << "_feature_names for writing\n";
        EXIT(EXIT_FAILURE);
    }
    if (scaffoldFile.fail()) {
        Debug(Debug::ERROR) << "Could not open " << contextDb << "_scaffolds for writing\n";
        EXIT(EXIT_FAILURE);
    }
}

uint32_t ContextDbWriter::writeFeature(const std::string &featureId, const std::string &targetSequenceId,
                                       DBKeyType targetKey, uint32_t scaffoldId, uint32_t start,
                                       uint32_t end, int8_t strand) {
    if (nextFeatureId == std::numeric_limits<uint32_t>::max()) {
        Debug(Debug::ERROR) << "Too many context features for uint32 identifiers\n";
        EXIT(EXIT_FAILURE);
    }

    uint64_t nameOffset = static_cast<uint64_t>(nameFile.tellp());
    nameFile.write(featureId.c_str(), featureId.size() + 1);
    if (targetKey == DB_KEY_INVALID) {
        nameFile.write(targetSequenceId.c_str(), targetSequenceId.size() + 1);
    }
    if (nameFile.fail()) {
        Debug(Debug::ERROR) << "Could not write context feature name\n";
        EXIT(EXIT_FAILURE);
    }

    std::string record;
    record.reserve(ContextDb::featureRecordSize());
    ContextDb::appendValue(record, scaffoldId);
    ContextDb::appendValue(record, start);
    ContextDb::appendValue(record, end);
    ContextDb::appendValue(record, targetKey);
    ContextDb::appendValue(record, nameOffset);
    ContextDb::appendValue(record, strand);

    featureFile.write(record.data(), record.size());
    if (featureFile.fail()) {
        Debug(Debug::ERROR) << "Could not write context feature record\n";
        EXIT(EXIT_FAILURE);
    }

    return nextFeatureId++;
}

void ContextDbWriter::writeScaffold(const std::string &scaffold) {
    scaffoldFile.write(scaffold.c_str(), scaffold.size() + 1);
    if (scaffoldFile.fail()) {
        Debug(Debug::ERROR) << "Could not write context scaffold\n";
        EXIT(EXIT_FAILURE);
    }
}

ContextDbReader::ContextDbReader(const std::string &contextDb) :
    featureFile((contextDb + "_features").c_str(), std::ios::binary),
    nameFile((contextDb + "_feature_names").c_str(), std::ios::binary),
    scaffolds(readScaffolds(contextDb + "_scaffolds")) {
    if (featureFile.fail()) {
        Debug(Debug::ERROR) << "File " << contextDb << "_features not found!\n";
        EXIT(EXIT_FAILURE);
    }
    if (nameFile.fail()) {
        Debug(Debug::ERROR) << "File " << contextDb << "_feature_names not found!\n";
        EXIT(EXIT_FAILURE);
    }
}

ContextFeatureRecord ContextDbReader::readFeature(uint32_t featureId) {
    std::string data(ContextDb::featureRecordSize(), '\0');
    featureFile.seekg(static_cast<std::streamoff>(static_cast<uint64_t>(featureId) *
                                                  ContextDb::featureRecordSize()));
    featureFile.read(&data[0], data.size());
    if (featureFile.fail()) {
        Debug(Debug::ERROR) << "Could not read context feature " << featureId << "\n";
        EXIT(EXIT_FAILURE);
    }

    const char *ptr = data.data();
    const char *end = ptr + data.size();
    ContextFeatureRecord record;
    record.scaffoldId = ContextDb::readValue<uint32_t>(ptr, end);
    record.start = ContextDb::readValue<uint32_t>(ptr, end);
    record.end = ContextDb::readValue<uint32_t>(ptr, end);
    record.targetKey = ContextDb::readValue<DBKeyType>(ptr, end);
    record.featureNameOffset = ContextDb::readValue<uint64_t>(ptr, end);
    record.strand = ContextDb::readValue<int8_t>(ptr, end);
    return record;
}

std::pair<std::string, std::string> ContextDbReader::readNames(const ContextFeatureRecord &record) {
    nameFile.clear();
    nameFile.seekg(static_cast<std::streamoff>(record.featureNameOffset));
    std::string featureId = readNullTerminatedName(nameFile);
    std::string targetId = record.targetKey == DB_KEY_INVALID ? readNullTerminatedName(nameFile) : "";
    return std::make_pair(featureId, targetId);
}

const std::string &ContextDbReader::scaffoldName(uint32_t scaffoldId) const {
    if (scaffoldId >= scaffolds.size()) {
        Debug(Debug::ERROR) << "Missing context scaffold " << scaffoldId << "\n";
        EXIT(EXIT_FAILURE);
    }
    return scaffolds[scaffoldId];
}

std::string ContextDbReader::readNullTerminatedName(std::ifstream &input) {
    std::string name;
    std::getline(input, name, '\0');
    if (input.good() == false) {
        Debug(Debug::ERROR) << "Could not read context name\n";
        EXIT(EXIT_FAILURE);
    }
    return name;
}

std::vector<std::string> ContextDbReader::readScaffolds(const std::string &filename) {
    std::ifstream input(filename.c_str(), std::ios::binary);
    if (input.fail()) {
        Debug(Debug::ERROR) << "File " << filename << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::vector<std::string> scaffolds;
    while (input.peek() != std::ifstream::traits_type::eof()) {
        scaffolds.push_back(readNullTerminatedName(input));
    }
    if (input.bad() || input.eof() == false) {
        Debug(Debug::ERROR) << "Could not read context scaffolds\n";
        EXIT(EXIT_FAILURE);
    }
    return scaffolds;
}
