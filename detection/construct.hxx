#ifndef LOG_ANOMALY_CONSTRUCT_HXX
#define LOG_ANOMALY_CONSTRUCT_HXX

#include "construct.hxx"
#include "gpsoinn/gpsoinn.hxx"
#include "graph/graph.hxx"
#include "infomap/Infomap.h"
#include "journal/sd_journal.hxx"
#include "stemmer/cstemmer"

#include "cereal/archives/binary.hpp"

#include "cereal/types/deque.hpp"
#include "cereal/types/unordered_map.hpp"
#include "cereal/types/utility.hpp"

#include "tclap/CmdLine.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
namespace {

class Extractor {
  private:
    using string = std::string;
    enum State {
        s_space = 0,
        s_word,
        s_zero,
        s_zerox,
        s_b16,
        s_b10,
        s_b16_c,
    };

  public:
    enum Type {
        eof = 0,
        newline,
        base16,
        base10,
        c_base16,
        word,
    };
    using result_type = std::pair<Type, string>;

    Extractor(string &&str, size_t beg) : m_text(str), m_word_beg(beg) {}
    ~Extractor() {}

    result_type extract();

  private:
    inline bool is_base16(char ch) {
        return (ch <= 'f' && ch >= 'a') || (ch >= '0' && ch <= '9');
    };
    inline bool is_base10(char ch) { return ch >= '0' && ch <= '9'; };

    inline result_type construct_result(size_t pos) {
        Type t;
        if (m_state > 3)
            t = static_cast<Type>(m_state - 2);
        else if (m_state == s_zero)
            t = base10;
        else
            t = word;
        return {t, m_text.substr(m_word_beg, pos - m_word_beg)};
    }

    string m_text;
    size_t m_word_beg;
    State m_state = s_space;
};
struct Message {
    using string = std::string;
    using WordMap = std::unordered_map<string, size_t>;
    string cursor;
    WordMap map;
    size_t count;
    float average;
    float variance;
};
struct ClusterInfo {
    typedef std::vector<std::string> vector;
    typedef std::unique_ptr<LogAnomaly::GPNet<0>> P_GPNet;
    vector tokens;
    P_GPNet gpnet;
};

std::pair<std::vector<double>, unsigned> construct_gpnet_input(
    std::unordered_map<std::string, size_t> const &map,
    std::unordered_map<std::string, size_t> const &word_map,
    LogAnomaly::UndirectedGraph<std::pair<unsigned, unsigned>, size_t> const
        &word_graph,
    std::unordered_map<unsigned, ClusterInfo> &clusters);

constexpr unsigned block_number = 6;
constexpr unsigned block_size = 2000;
constexpr float vertex_fade_coeff = 2;
constexpr float edge_fade_coeff = 2.5;
constexpr unsigned vertex_fade_threshold = 1;
constexpr unsigned edge_fade_threshold = 1;

void process_log(std::string log_text,
                 std::unordered_map<std::string, size_t> &word_map,
                 LogAnomaly::UndirectedGraph<std::pair<unsigned, unsigned>,
                                             size_t> &word_graph,
                 std::deque<Message> &messages,
                 std::unordered_set<std::string> const &stop_words,
                 LogAnomaly::SDJournal &journal);
void fade_out(std::unordered_map<std::string, size_t> &word_map,
              LogAnomaly::UndirectedGraph<std::pair<unsigned, unsigned>, size_t>
                  &word_graph);

} // namespace
CEREAL_REGISTER_TYPE(
    LogAnomaly::UndirectedGraph<
        LogAnomaly::GPNet<0>::Node, unsigned,
        std::less<LogAnomaly::GPNet<0>::NodeVector>,
        Eigen::aligned_allocator<LogAnomaly::GPNet<0>::NodeVector>>);
CEREAL_REGISTER_TYPE(
    LogAnomaly::UndirectedGraph<std::pair<unsigned, unsigned>, size_t>);
#endif // LOG_ANOMALY_CONSTRUCT_HXX
