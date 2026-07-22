#include <set>
#include <string>
#include <utility>

#include "ContextDb.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"

static std::string lookupName(DBReader<DBKeyType> &reader, DBKeyType key) {
    size_t id = reader.getLookupIdByKey(key);
    return id == SIZE_MAX ? SSTR(key) : reader.getLookupEntryName(id);
}

static void appendFeature(std::string &out, DBKeyType hitTargetKey, const std::string &anchorFeatureId,
                          int relPos, const ContextFeatureRecord &feature, ContextDbReader &contextDb,
                          DBReader<DBKeyType> &targetLookup) {
    const std::pair<std::string, std::string> names = contextDb.readNames(feature);
    const bool mapped = feature.targetKey != DB_KEY_INVALID;
    out.append(SSTR(hitTargetKey)).append("\t")
       .append(anchorFeatureId).append("\t")
       .append(SSTR(relPos)).append("\t")
       .append(mapped ? SSTR(feature.targetKey) : "-").append("\t")
       .append(mapped ? lookupName(targetLookup, feature.targetKey) : names.second).append("\t")
       .append(names.first).append("\t")
       .append(contextDb.scaffoldName(feature.scaffoldId)).append("\t")
       .append(SSTR(feature.start)).append("\t")
       .append(SSTR(feature.end)).append("\t")
       .append(feature.strand < 0 ? "-1\n" : "1\n");
}

static void appendContexts(std::string &out, DBKeyType hitTargetKey, char *payload, size_t payloadLength,
                           ContextDbReader &contextDb, DBReader<DBKeyType> &targetLookup,
                           std::set<std::pair<DBKeyType, uint32_t> > &seenContexts) {
    if (payloadLength == 0 || payload[payloadLength - 1] != '\0' ||
        (payloadLength - 1) % ContextDb::contextRecordSize() != 0) {
        Debug(Debug::ERROR) << "Invalid context payload length\n";
        EXIT(EXIT_FAILURE);
    }
    const char *ptr = payload;
    const char *end = payload + payloadLength - 1;
    while (ptr < end) {
        uint32_t firstFeatureId;
        uint8_t leftCount;
        uint8_t rightCount;
        ContextDb::readContext(ptr, end, firstFeatureId, leftCount, rightCount);
        const uint32_t anchorId = firstFeatureId + leftCount;
        if (seenContexts.insert(std::make_pair(hitTargetKey, anchorId)).second == false) {
            continue;
        }

        const std::string anchorFeatureId = contextDb.readNames(contextDb.readFeature(anchorId)).first;
        const uint16_t featureCount = static_cast<uint16_t>(leftCount) + rightCount + 1;
        for (uint16_t featurePos = 0; featurePos < featureCount; ++featurePos) {
            appendFeature(out, hitTargetKey, anchorFeatureId,
                          static_cast<int>(featurePos) - static_cast<int>(leftCount),
                          contextDb.readFeature(firstFeatureId + featurePos), contextDb, targetLookup);
        }
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

    const std::string contextDb = par.db2 + "_context";
    DBReader<DBKeyType> contextReader(contextDb.c_str(), (contextDb + ".index").c_str(), 1,
                                      DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    contextReader.open(DBReader<DBKeyType>::NOSORT);
    ContextDbReader contextDbReader(contextDb);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), 1, false, Parameters::DBTYPE_GENERIC_DB);
    writer.open();

    std::string output;
    output.reserve(4096);
    Debug::Progress progress(resultReader.getSize());
    for (size_t i = 0; i < resultReader.getSize(); ++i) {
        progress.updateProgress();
        std::set<std::pair<DBKeyType, uint32_t> > seenContexts;
        char *data = resultReader.getData(i, 0);
        while (*data != '\0') {
            const DBKeyType targetKey = Util::fast_atoi<DBKeyType>(data);
            const size_t contextEntryId = contextReader.getId(targetKey);
            if (contextEntryId != SIZE_MAX) {
                appendContexts(output, targetKey, contextReader.getData(contextEntryId, 0),
                               contextReader.getEntryLen(contextEntryId), contextDbReader, targetLookup, seenContexts);
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
