#ifndef LOG_ANOMALY_SD_JOURNAL_HXX
#define LOG_ANOMALY_SD_JOURNAL_HXX

extern "C" {
#include <systemd/sd-journal.h>
}
#include <bitset>
#include <string>
#include <system_error>
#include <utility>

namespace LogAnomaly {
class SDJournal {
  private:
    enum State {
        s_fail = 0,
        s_bad,
        s_eof,
        s_invalidated,
    };

  public:
    // constructors
    explicit SDJournal(int flags = 0) {
        using std::generic_category;
        using std::system_error;

        int result = sd_journal_open(&m_context, flags);
        if (result < 0) {
            m_states[s_bad] = 1;
            throw system_error(-result, generic_category(),
                               "Failed to open journal");
        }
    }
    SDJournal(const SDJournal &other) = delete;
    SDJournal(SDJournal &&other) {
        m_context = other.m_context;
        m_states = other.m_states;

        other.m_context = nullptr;
    }
    // swap
    friend void swap(SDJournal &left, SDJournal &right) {
        using std::swap;

        swap(left.m_context, right.m_context);
        swap(left.m_states, right.m_states);
    }
    // destructor
    ~SDJournal() { sd_journal_close(m_context); }

    // operators
    friend SDJournal &operator>>(SDJournal &object, std::string &str);
    SDJournal &operator=(SDJournal other) {
        std::swap(*this, other);
        return *this;
    }

    uint64_t timeout() const { return m_timeout; }
    void set_timeout(uint64_t timeout_usec = 0) { m_timeout = timeout_usec; }

    // states
    typedef unsigned char journal_state;

    static constexpr journal_state goodbit = 0;
    static constexpr journal_state failbit = 1 << s_fail;
    static constexpr journal_state badbit = 1 << s_bad;
    static constexpr journal_state eofbit = 1 << s_eof;
    static constexpr journal_state invalidatedbit = 1 << s_invalidated;

    bool good() const { return m_states.none(); }
    bool invalidated() const { return m_states[s_invalidated]; }
    bool eof() const { return m_states[s_eof]; }
    bool bad() const { return m_states[s_bad]; }
    bool fail() const { return m_states[s_fail] | m_states[s_bad]; }
    void clear(journal_state state = goodbit) { m_states = state; }
    void setstate(journal_state state) { m_states |= state; }
    journal_state rdstate() const { return m_states.to_ulong(); }

    explicit operator bool() const noexcept { return good(); }

    typedef std::string journalpos;
    journalpos tellg() {
        using std::generic_category;
        using std::string;
        using std::system_error;

        char *c_str;
        int result = sd_journal_get_cursor(m_context, &c_str);
        if (result < 0) {
            throw system_error(-result, generic_category(),
                               "Failed to get cursor");
        }

        string str(c_str);
        free(c_str);
        return str;
    }

    SDJournal &seekg(journalpos pos) {
        using std::generic_category;
        using std::system_error;

        int result = sd_journal_seek_cursor(m_context, pos.data());
        if (result < 0) {
            throw system_error(-result, generic_category(),
                               "Failed to seek to the position");
        }
        return *this;
    }

  private:
    sd_journal *m_context = nullptr;

    // state bits
    // |  3    |  2  |  1  |  0   |
    // |invalid| eof | bad | fail |
    std::bitset<4> m_states;
    uint64_t m_timeout;
};
} // namespace LogAnomaly

#endif // LOG_ANOMALY_SD_JOURNAL_HXX
