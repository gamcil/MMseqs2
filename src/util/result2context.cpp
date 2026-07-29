#include <cstdint>
#include <limits>
#include <set>
#include <string>

#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"

struct ContextField {
    const char *data;
    size_t length;
};

struct ContextFeatureView {
    ContextField targetKey;
    ContextField targetSequenceId;
    ContextField featureId;
    ContextField start;
    ContextField end;
    ContextField strand;
};

// Each context DB payload contains one or more newline-delimited rows:
// scaffold, left count, right count, then six fields per feature:
// target key, target sequence ID, feature ID, start, end, strand.

static std::string lookupName(DBReader<DBKeyType> &reader, DBKeyType key) {
    size_t id = reader.getLookupIdByKey(key);
    return id == SIZE_MAX ? SSTR(key) : reader.getLookupEntryName(id);
}

static bool readField(const char *&ptr, const char *end, ContextField &field) {
    if (ptr >= end) {
        return false;
    }
    field.data = ptr;
    while (ptr < end && *ptr != '\t') {
        ptr++;
    }
    field.length = static_cast<size_t>(ptr - field.data);
    if (ptr < end) {
        ptr++;
    }
    return field.length != 0;
}

static bool readFeature(const char *&ptr, const char *end, ContextFeatureView &feature) {
    return readField(ptr, end, feature.targetKey) &&
           readField(ptr, end, feature.targetSequenceId) &&
           readField(ptr, end, feature.featureId) &&
           readField(ptr, end, feature.start) &&
           readField(ptr, end, feature.end) &&
           readField(ptr, end, feature.strand);
}

static uint16_t parseCount(const ContextField &field) {
    uint32_t value = 0;
    for (size_t i = 0; i < field.length; ++i) {
        if (field.data[i] < '0' || field.data[i] > '9') {
            Debug(Debug::ERROR) << "Invalid context neighbor count\n";
            EXIT(EXIT_FAILURE);
        }
        value = value * 10 + static_cast<uint32_t>(field.data[i] - '0');
        if (value > 255) {
            Debug(Debug::ERROR) << "Invalid context neighbor count\n";
            EXIT(EXIT_FAILURE);
        }
    }
    return static_cast<uint16_t>(value);
}

static DBKeyType parseTargetKey(const ContextField &field) {
    DBKeyType value = 0;
    for (size_t i = 0; i < field.length; ++i) {
        if (field.data[i] < '0' || field.data[i] > '9') {
            Debug(Debug::ERROR) << "Invalid mapped context target key\n";
            EXIT(EXIT_FAILURE);
        }
        const DBKeyType digit = static_cast<DBKeyType>(field.data[i] - '0');
        if (value > (std::numeric_limits<DBKeyType>::max() - digit) / 10) {
            Debug(Debug::ERROR) << "Invalid mapped context target key\n";
            EXIT(EXIT_FAILURE);
        }
        value = value * 10 + digit;
    }
    return value;
}

static void appendField(std::string &out, const ContextField &field) {
    out.append(field.data, field.length);
}

static void appendFeature(std::string &out, DBKeyType hitTargetKey, const ContextField &anchorFeatureId,
                          int relPos, const ContextField &scaffold, const ContextFeatureView &feature,
                          DBReader<DBKeyType> &targetLookup) {
    const bool mapped = !(feature.targetKey.length == 1 && feature.targetKey.data[0] == '-');
    out.append(SSTR(hitTargetKey)).append("\t");
    appendField(out, anchorFeatureId);
    out.append("\t").append(SSTR(relPos)).append("\t");
    appendField(out, feature.targetKey);
    out.append("\t");
    if (mapped) {
        out.append(lookupName(targetLookup, parseTargetKey(feature.targetKey)));
    } else {
        appendField(out, feature.targetSequenceId);
    }
    out.append("\t");
    appendField(out, feature.featureId);
    out.append("\t");
    appendField(out, scaffold);
    out.append("\t");
    appendField(out, feature.start);
    out.append("\t");
    appendField(out, feature.end);
    out.append("\t");
    appendField(out, feature.strand);
    out.push_back('\n');
}

static void appendContextLine(std::string &out, DBKeyType hitTargetKey, const char *line,
                              const char *lineEnd, DBReader<DBKeyType> &targetLookup) {
    ContextField scaffold;
    ContextField leftField;
    ContextField rightField;
    const char *ptr = line;
    if (!readField(ptr, lineEnd, scaffold) ||
        !readField(ptr, lineEnd, leftField) ||
        !readField(ptr, lineEnd, rightField)) {
        Debug(Debug::ERROR) << "Invalid context row header\n";
        EXIT(EXIT_FAILURE);
    }

    const uint16_t leftCount = parseCount(leftField);
    const uint16_t rightCount = parseCount(rightField);
    const uint16_t featureCount = static_cast<uint16_t>(leftCount + rightCount + 1);

    const char *featuresStart = ptr;
    ContextField anchorFeatureId;
    for (uint16_t i = 0; i <= leftCount; ++i) {
        ContextFeatureView feature;
        if (!readFeature(ptr, lineEnd, feature)) {
            Debug(Debug::ERROR) << "Invalid context feature row\n";
            EXIT(EXIT_FAILURE);
        }
        if (i == leftCount) {
            anchorFeatureId = feature.featureId;
        }
    }

    ptr = featuresStart;
    for (uint16_t i = 0; i < featureCount; ++i) {
        ContextFeatureView feature;
        if (!readFeature(ptr, lineEnd, feature)) {
            Debug(Debug::ERROR) << "Invalid context feature row\n";
            EXIT(EXIT_FAILURE);
        }
        appendFeature(out, hitTargetKey, anchorFeatureId,
                      static_cast<int>(i) - static_cast<int>(leftCount),
                      scaffold, feature, targetLookup);
    }
    if (ptr != lineEnd) {
        Debug(Debug::ERROR) << "Unexpected fields in context row\n";
        EXIT(EXIT_FAILURE);
    }
}

static void appendContexts(std::string &out, DBKeyType hitTargetKey, char *payload, size_t payloadLength,
                           DBReader<DBKeyType> &targetLookup) {
    if (payloadLength == 0 || payload[payloadLength - 1] != '\0') {
        Debug(Debug::ERROR) << "Invalid context payload\n";
        EXIT(EXIT_FAILURE);
    }

    const char *ptr = payload;
    const char *end = payload + payloadLength - 1;
    while (ptr < end) {
        const char *lineEnd = ptr;
        while (lineEnd < end && *lineEnd != '\n') {
            lineEnd++;
        }
        if (lineEnd == end) {
            Debug(Debug::ERROR) << "Unterminated context row\n";
            EXIT(EXIT_FAILURE);
        }
        appendContextLine(out, hitTargetKey, ptr, lineEnd, targetLookup);
        ptr = lineEnd + 1;
    }
}

int result2context(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<DBKeyType> targetLookup(par.db2.c_str(), par.db2Index.c_str(), 1, DBReader<DBKeyType>::USE_LOOKUP);
    targetLookup.open(DBReader<DBKeyType>::NOSORT);

    DBReader<DBKeyType> resultReader(par.db3.c_str(), par.db3Index.c_str(), 1,
                                     DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    resultReader.open(DBReader<DBKeyType>::LINEAR_ACCCESS);

    DBReader<DBKeyType> contextReader(par.db4.c_str(), par.db4Index.c_str(), 1,
                                      DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    contextReader.open(DBReader<DBKeyType>::NOSORT);

    DBWriter writer(par.db5.c_str(), par.db5Index.c_str(), 1, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writer.open();

    std::string output;
    output.reserve(4096);
    Debug::Progress progress(resultReader.getSize());
    for (size_t i = 0; i < resultReader.getSize(); ++i) {
        progress.updateProgress();
        std::set<DBKeyType> seenTargets;
        char *data = resultReader.getData(i, 0);
        while (*data != '\0') {
            const DBKeyType targetKey = Util::fast_atoi<DBKeyType>(data);
            if (seenTargets.insert(targetKey).second) {
                const size_t contextEntryId = contextReader.getId(targetKey);
                if (contextEntryId != SIZE_MAX) {
                    appendContexts(output, targetKey, contextReader.getData(contextEntryId, 0),
                                   contextReader.getEntryLen(contextEntryId), targetLookup);
                }
            }
            data = Util::skipLine(data);
        }
        writer.writeData(output.data(), output.size(), resultReader.getDbKey(i), 0);
        output.clear();
    }

    writer.close();
    contextReader.close();
    resultReader.close();
    targetLookup.close();
    return EXIT_SUCCESS;
}
