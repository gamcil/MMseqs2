#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "CommandCaller.h"
#include "ContextDb.h"
#include "DBWriter.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Parameters.h"
#include "Util.h"
#include "createcontextdb.sh.h"

struct ContextFeature {
    std::string targetSequenceId;
    std::string featureId;
    std::string scaffold;
    uint32_t start;
    uint32_t end;
    int8_t strand;
    DBKeyType dbKey;
    uint32_t scaffoldId;
    uint32_t compactId;
    std::vector<DBKeyType> anchorKeys;
};

static const uint32_t FEATURE_ID_INVALID = std::numeric_limits<uint32_t>::max();

static std::string fieldString(const char *field) {
    return std::string(field, Util::skipNonTab(field));
}

static void failInvalidFeatureLine(size_t lineNumber, const std::string &message) {
    Debug(Debug::ERROR) << "Invalid feature table line " << lineNumber << ": " << message << "\n";
    EXIT(EXIT_FAILURE);
}

static void requireFeatureField(const std::string &value, size_t lineNumber, const char *fieldName) {
    if (value.empty()) {
        failInvalidFeatureLine(lineNumber, std::string("empty ") + fieldName);
    }
}

template <typename T>
static T parseUnsignedFeatureField(const std::string &value, size_t lineNumber, const char *fieldName) {
    requireFeatureField(value, lineNumber, fieldName);
    if (value[0] < '0' || value[0] > '9') {
        failInvalidFeatureLine(lineNumber, std::string("invalid numeric value for ") + fieldName);
    }
    char *end = NULL;
    errno = 0;
    const unsigned long long parsed = strtoull(value.c_str(), &end, 10);
    if (errno == ERANGE || *end != '\0' ||
        parsed > static_cast<unsigned long long>(std::numeric_limits<T>::max())) {
        failInvalidFeatureLine(lineNumber, std::string("invalid numeric value for ") + fieldName);
    }
    return static_cast<T>(parsed);
}

static int8_t parseStrandFeatureField(const std::string &value, size_t lineNumber) {
    if (value == "1" || value == "+") {
        return 1;
    }
    if (value == "-1" || value == "-") {
        return -1;
    }
    failInvalidFeatureLine(lineNumber, "strand must be one of 1, -1, +, -");
    return 1;
}

static bool isValidSortMemory(const std::string &value) {
    const size_t suffix = value.find_first_not_of("0123456789");
    if (suffix == 0) {
        return false;
    }
    if (suffix == std::string::npos) {
        return true;
    }
    static const std::string units = "%bBkKmMgGtTpPeEzZyY";
    return suffix + 1 == value.size() && units.find(value[suffix]) != std::string::npos;
}

static void validateContextWindow(int window) {
    if (window > std::numeric_limits<uint8_t>::max()) {
        Debug(Debug::ERROR) << "--context-window must not exceed 255\n";
        EXIT(EXIT_FAILURE);
    }
}

static uint64_t featureLength(const ContextFeature &feature) {
    return static_cast<uint64_t>(feature.end) - feature.start + 1;
}

static bool isRedundantFeature(const ContextFeature &left, const ContextFeature &right, float threshold) {
    if (threshold <= 0.0f || left.scaffold != right.scaffold || left.strand != right.strand ||
        left.end < right.start || right.end < left.start) {
        return false;
    }
    const uint64_t overlap = static_cast<uint64_t>(std::min(left.end, right.end)) -
                             std::max(left.start, right.start) + 1;
    return static_cast<float>(overlap) /
           static_cast<float>(std::min(featureLength(left), featureLength(right))) >= threshold;
}

static void appendAnchorKeys(std::vector<DBKeyType> &target, const std::vector<DBKeyType> &source) {
    for (size_t i = 0; i < source.size(); ++i) {
        if (std::find(target.begin(), target.end(), source[i]) == target.end()) {
            target.push_back(source[i]);
        }
    }
}

static void collapseFeature(ContextFeature &kept, const ContextFeature &candidate) {
    std::vector<DBKeyType> anchorKeys = kept.anchorKeys;
    appendAnchorKeys(anchorKeys, candidate.anchorKeys);
    const bool keptMapped = kept.anchorKeys.empty() == false;
    const bool candidateMapped = candidate.anchorKeys.empty() == false;
    const bool preferCandidate = (keptMapped == false && candidateMapped) ||
                                 (keptMapped == candidateMapped &&
                                  featureLength(candidate) > featureLength(kept));
    if (preferCandidate) {
        kept = candidate;
    }
    kept.anchorKeys.swap(anchorKeys);
}

static uint32_t ensureFeatureId(ContextDbWriter &writer, ContextFeature &feature) {
    if (feature.compactId == FEATURE_ID_INVALID) {
        feature.compactId = writer.writeFeature(
            feature.featureId, feature.targetSequenceId, feature.dbKey,
            feature.scaffoldId, feature.start, feature.end, feature.strand
        );
    }
    return feature.compactId;
}

static void emitContext(std::ofstream &out, ContextDbWriter &writer, std::deque<ContextFeature> &features,
                        size_t anchorIndex, int window) {
    ContextFeature &anchor = features[anchorIndex];
    if (anchor.anchorKeys.empty()) return;  // skip non-anchor features
    size_t windowStart = anchorIndex > static_cast<size_t>(window) ? anchorIndex - window : 0;
    size_t windowEnd = std::min(features.size() - 1, anchorIndex + static_cast<size_t>(window));
    uint32_t firstFeatureId = ensureFeatureId(writer, features[windowStart]);
    for (size_t idx = windowStart + 1; idx <= windowEnd; ++idx) {
        ensureFeatureId(writer, features[idx]);
    }
    for (size_t i = 0; i < anchor.anchorKeys.size(); ++i) {
        // anchor db key \t first feature id in context \t #features left of anchor \t #features right of anchor
        out << anchor.anchorKeys[i] << '\t'
            << firstFeatureId << '\t'
            << anchorIndex - windowStart << '\t'
            << windowEnd - anchorIndex << '\n';
    }
    if (out.fail()) {
        Debug(Debug::ERROR) << "Could not write context row\n";
        EXIT(EXIT_FAILURE);
    }
}

static void queueFeature(std::ofstream &out, ContextDbWriter &writer, std::deque<ContextFeature> &features,
                         size_t &nextAnchor, const ContextFeature &feature, int window) {
    features.push_back(feature);
    while (nextAnchor + static_cast<size_t>(window) < features.size()) {
        emitContext(out, writer, features, nextAnchor, window);
        nextAnchor++;
        if (nextAnchor > static_cast<size_t>(window)) {
            features.pop_front();
            nextAnchor--;
        }
    }
}

static void flushScaffold(std::ofstream &out, ContextDbWriter &writer, std::deque<ContextFeature> &features,
                          size_t &nextAnchor, int window) {
    while (nextAnchor < features.size()) {
        emitContext(out, writer, features, nextAnchor, window);
        nextAnchor++;
    }
    features.clear();
    nextAnchor = 0;
}

static ContextFeature parseResolvedFeature(const std::string &line, size_t lineNumber) {
    const char *fields[7];
    if (Util::getFieldsOfLine(line.c_str(), fields, 7) < 7) {
        Debug(Debug::ERROR) << "Invalid resolved feature table line " << lineNumber << ": expected 7 tab-separated columns\n";
        EXIT(EXIT_FAILURE);
    }

    ContextFeature feature;
    const std::string dbKey = fieldString(fields[0]);
    feature.targetSequenceId = fieldString(fields[1]);
    feature.featureId = fieldString(fields[2]);
    feature.scaffold = fieldString(fields[3]);
    const std::string start = fieldString(fields[4]);
    const std::string end = fieldString(fields[5]);
    const std::string strand = fieldString(fields[6]);
    requireFeatureField(feature.targetSequenceId, lineNumber, "target_sequence_id");
    requireFeatureField(feature.featureId, lineNumber, "feature_id");
    requireFeatureField(feature.scaffold, lineNumber, "scaffold");
    feature.start = parseUnsignedFeatureField<uint32_t>(start, lineNumber, "start");
    feature.end = parseUnsignedFeatureField<uint32_t>(end, lineNumber, "end");
    if (feature.start > feature.end) {
        failInvalidFeatureLine(lineNumber, "start is greater than end");
    }
    feature.strand = parseStrandFeatureField(strand, lineNumber);
    feature.dbKey = dbKey == "-" ? DB_KEY_INVALID : parseUnsignedFeatureField<DBKeyType>(dbKey, lineNumber, "db key");
    feature.scaffoldId = 0;
    feature.compactId = FEATURE_ID_INVALID;
    if (feature.dbKey != DB_KEY_INVALID) {
        feature.anchorKeys.push_back(feature.dbKey);
    }
    return feature;
}

static bool readNextLookupEntry(std::ifstream &input, size_t &lineNumber, const int idMode,
                                std::string &lookupId, std::string &lookupKey) {
    std::string line;
    while (std::getline(input, line)) {
        lineNumber++;
        if (line.empty()) continue;
        const char *fields[3];
        if (Util::getFieldsOfLine(line.c_str(), fields, 3) < 3) {
            Debug(Debug::ERROR) << "Invalid lookup line " << lineNumber << ": expected at least 3 tab-separated columns\n";
            EXIT(EXIT_FAILURE);
        }
        lookupKey = fieldString(fields[0]);
        lookupId = idMode == 0 ? lookupKey : fieldString(fields[1]);
        return true;
    }
    return false;
}

int createcontextresolve(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::ifstream lookupInput(par.db1);
    if (lookupInput.fail()) {
        Debug(Debug::ERROR) << "File " << par.db1 << " not found!\n";
        EXIT(EXIT_FAILURE);
    }
    std::ifstream featureInput(par.db2);
    if (featureInput.fail()) {
        Debug(Debug::ERROR) << "File " << par.db2 << " not found!\n";
        EXIT(EXIT_FAILURE);
    }
    std::ofstream out(par.db3);
    if (out.fail()) {
        Debug(Debug::ERROR) << "Could not open " << par.db3 << " for writing!\n";
        EXIT(EXIT_FAILURE);
    }

    size_t lookupLineNumber = 0;
    std::string lookupId;
    std::string lookupKey;
    bool hasLookupEntry = readNextLookupEntry(lookupInput, lookupLineNumber, par.contextIdMode,
                                              lookupId, lookupKey);

    std::string featureLine;
    std::string previousTargetId;
    std::string resolvedKey;
    size_t featureLineNumber = 0;
    size_t resolvedAnchors = 0;
    size_t featureCount = 0;
    bool resolved = false;

    while (std::getline(featureInput, featureLine)) {
        featureLineNumber++;
        if (featureLine.empty() || featureLine[0] == '#') {
            continue;
        }
        const char *fields[6];
        if (Util::getFieldsOfLine(featureLine.c_str(), fields, 6) < 6) {
            Debug(Debug::ERROR) << "Invalid feature table line " << featureLineNumber << ": expected 6 tab-separated columns\n";
            EXIT(EXIT_FAILURE);
        }
        std::string targetId = fieldString(fields[0]);
        if (targetId == "target_sequence_id" && fieldString(fields[1]) == "feature_id") {
            continue;
        }
        requireFeatureField(targetId, featureLineNumber, "target_sequence_id");
        if (targetId != previousTargetId) {
            while (hasLookupEntry && lookupId < targetId) {
                hasLookupEntry = readNextLookupEntry(lookupInput, lookupLineNumber, par.contextIdMode, lookupId, lookupKey);
            }
            resolved = hasLookupEntry && lookupId == targetId;
            resolvedKey = resolved ? lookupKey : "-";
            if (resolved) {
                hasLookupEntry = readNextLookupEntry(lookupInput, lookupLineNumber, par.contextIdMode, lookupId, lookupKey);
                if (par.contextIdMode == 1 && hasLookupEntry && lookupId == targetId) {
                    Debug(Debug::ERROR) << "Duplicate sequence identifier \"" << lookupId
                                        << "\" in lookup at line " << lookupLineNumber
                                        << "; use --context-id-mode 0 with database keys\n";
                    EXIT(EXIT_FAILURE);
                }
            }
            previousTargetId = targetId;
        }
        out << resolvedKey << '\t' << featureLine << '\n';
        if (resolved) {
            resolvedAnchors++;
        }
        if (out.fail()) {
            Debug(Debug::ERROR) << "Could not write resolved feature table\n";
            EXIT(EXIT_FAILURE);
        }
        featureCount++;
    }
    Debug(Debug::INFO) << "Resolved " << resolvedAnchors << " of " << featureCount << " features as context anchors\n";
    return EXIT_SUCCESS;
}

int createcontextcontexts(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    validateContextWindow(par.contextWindow);

    std::ifstream input(par.db1);
    if (input.fail()) {
        Debug(Debug::ERROR) << "File " << par.db1 << " not found!\n";
        EXIT(EXIT_FAILURE);
    }
    std::ofstream out(par.db2);
    if (out.fail()) {
        Debug(Debug::ERROR) << "Could not open " << par.db2 << " for writing!\n";
        EXIT(EXIT_FAILURE);
    }

    const std::string contextDb = par.db3;
    ContextDbWriter contextWriter(contextDb);

    std::string line;
    size_t lineNumber = 0;
    std::deque<ContextFeature> features;
    std::string currentScaffold;
    uint32_t currentScaffoldId = 0;
    size_t nextAnchor = 0;
    bool haveScaffold = false;
    bool scaffoldWritten = false;
    bool havePending = false;
    ContextFeature pending;

    while (std::getline(input, line)) {
        lineNumber++;
        if (line.empty()) {
            continue;
        }
        ContextFeature feature = parseResolvedFeature(line, lineNumber);
        if (haveScaffold == false || feature.scaffold != currentScaffold) {
            if (haveScaffold) {
                queueFeature(out, contextWriter, features, nextAnchor, pending, par.contextWindow);
                havePending = false;
                flushScaffold(out, contextWriter, features, nextAnchor, par.contextWindow);
                if (scaffoldWritten) {
                    currentScaffoldId++;
                }
            }
            currentScaffold = feature.scaffold;
            haveScaffold = true;
            scaffoldWritten = false;
        }
        feature.scaffoldId = currentScaffoldId;
        if (scaffoldWritten == false && feature.anchorKeys.empty() == false) {
            contextWriter.writeScaffold(currentScaffold);
            scaffoldWritten = true;
        }

        if (havePending && isRedundantFeature(pending, feature, par.contextCollapseOverlap)) {
            collapseFeature(pending, feature);
        } else {
            if (havePending) {
                queueFeature(out, contextWriter, features, nextAnchor, pending, par.contextWindow);
            }
            pending = feature;
            havePending = true;
        }
    }

    if (havePending) {
        queueFeature(out, contextWriter, features, nextAnchor, pending, par.contextWindow);
    }
    flushScaffold(out, contextWriter, features, nextAnchor, par.contextWindow);
    return EXIT_SUCCESS;
}

int createcontextdbcore(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::ifstream input(par.db1);
    if (input.fail()) {
        Debug(Debug::ERROR) << "File " << par.db1 << " not found!\n";
        EXIT(EXIT_FAILURE);
    }
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), 1, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writer.open();

    std::string line;
    bool haveKey = false;
    DBKeyType currentKey = 0;
    std::string payload;

    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        const char *fields[4];
        if (Util::getFieldsOfLine(line.c_str(), fields, 4) < 4) {
            Debug(Debug::ERROR) << "Invalid context row\n";
            EXIT(EXIT_FAILURE);
        }

        DBKeyType key = Util::fast_atoi<DBKeyType>(fields[0]);
        if (haveKey && key != currentKey) {
            writer.writeData(payload.data(), payload.size(), currentKey, 0);
            payload.clear();
        }
        haveKey = true;
        currentKey = key;

        const uint16_t leftCount = Util::fast_atoi<uint16_t>(fields[2]);
        const uint16_t rightCount = Util::fast_atoi<uint16_t>(fields[3]);
        if (leftCount > std::numeric_limits<uint8_t>::max() ||
            rightCount > std::numeric_limits<uint8_t>::max()) {
            Debug(Debug::ERROR) << "Invalid context row: neighbor count exceeds 255\n";
            EXIT(EXIT_FAILURE);
        }
        ContextDb::appendContext(payload, Util::fast_atoi<uint32_t>(fields[1]),
                                 static_cast<uint8_t>(leftCount), static_cast<uint8_t>(rightCount));
    }

    if (haveKey) {
        writer.writeData(payload.data(), payload.size(), currentKey, 0);
    }

    writer.close();
    return EXIT_SUCCESS;
}

int createcontextdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    validateContextWindow(par.contextWindow);

    if (par.sortMemory.empty() == false && isValidSortMemory(par.sortMemory) == false) {
        Debug(Debug::ERROR) << "Invalid --sort-memory value. Use a single sort -S value such as 50% or 16G.\n";
        EXIT(EXIT_FAILURE);
    }

    std::string tmp = par.filenames.back();
    std::string contextDb = par.db1 + "_context";
    if (FileUtil::directoryExists(tmp.c_str()) == false && FileUtil::makeDir(tmp.c_str()) == false) {
        Debug(Debug::ERROR) << "Can not create tmp folder " << tmp << ".\n";
        EXIT(EXIT_FAILURE);
    }

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("CONTEXT_ID_MODE", SSTR(par.contextIdMode).c_str());
    cmd.addVariable("SORT_MEMORY", par.sortMemory.empty() ? NULL : par.sortMemory.c_str());
    cmd.addVariable("CREATECONTEXTRESOLVE_PAR", par.createParameterString(par.createcontextresolve, true).c_str());
    cmd.addVariable("CREATECONTEXTCONTEXTS_PAR", par.createParameterString(par.createcontextcontexts, true).c_str());
    cmd.addVariable("CREATECONTEXTDBCORE_PAR", par.createParameterString(par.createcontextdbcore, true).c_str());

    std::string program = tmp + "/createcontextdb.sh";
    FileUtil::writeFile(program, createcontextdb_sh, createcontextdb_sh_len);
    std::vector<std::string> workflowArgs;
    workflowArgs.push_back(par.db1);
    workflowArgs.push_back(par.db2);
    workflowArgs.push_back(contextDb);
    workflowArgs.push_back(tmp);
    cmd.execProgram(program.c_str(), workflowArgs);

    return EXIT_SUCCESS;
}
