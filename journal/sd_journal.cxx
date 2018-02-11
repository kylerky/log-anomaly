#include "sd_journal.hxx"

namespace LogAnomaly {
SDJournal &operator>>(SDJournal &object, std::string &str) {
    using std::generic_category;
    using std::string;
    using std::system_error;

    if (!object)
        return object;
    int result = sd_journal_next(object.m_context);
    if (result < 0) {
        object.m_states[object.s_bad] = 1;
        // throw system_error(-result, generic_category(),
        //                    "Failed to iterate to the next entry");
        return object;
    } else {
        if (result == 0) {
            if (object.m_timeout > 0) {
                result = sd_journal_wait(object.m_context, object.m_timeout);
                switch (result) {
                case SD_JOURNAL_INVALIDATE:
                    object.m_states[object.s_invalidated] = 1;
                    return object;
                case SD_JOURNAL_NOP:
                    object.m_states[object.s_eof] = 1;
                    return object;
                case SD_JOURNAL_APPEND:
                    result = sd_journal_next(object.m_context);
                    if (result < 0) {
                        object.m_states[object.s_bad] = 1;
                        return object;
                    }
                    break;
                default:
                    object.m_states[object.s_bad] = 1;
                    return object;
                }
            } else {
                object.m_states[object.s_eof] = 1;
                return object;
            }
        }
        const void *data;
        size_t length;

        result =
            sd_journal_get_data(object.m_context, "MESSAGE", &data, &length);
        if (result < 0) {
            object.m_states[object.s_bad] = 1;
            // throw system_error(-result, generic_category(),
            //                    "Failed to read message field");
            return object;
        }

        str = string(static_cast<const char *>(data), length);
    }
    return object;
}

} // namespace LogAnomaly
