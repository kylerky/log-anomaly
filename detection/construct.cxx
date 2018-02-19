#include "construct.hxx"
int main() {
    using LogAnomaly::GPNet;
    using LogAnomaly::PorterStemmer::stem;
    using LogAnomaly::SDJournal;
    using LogAnomaly::UndirectedGraph;
    using infomap::Infomap;
    using std::accumulate;
    using std::cerr;
    using std::cout;
    using std::deque;
    using std::endl;
    using std::isalnum;
    using std::isinf;
    using std::isnan;
    using std::istringstream;
    using std::make_unique;
    using std::map;
    using std::multimap;
    using std::pair;
    using std::pow;
    using std::replace_if;
    using std::stoull;
    using std::string;
    using std::unordered_map;
    using std::unordered_set;
    using std::vector;

    const unordered_set<string> stop_words = {
        "a",      "able",    "about", "across",  "after", "all",    "almost",
        "also",   "am",      "among", "an",      "and",   "any",    "are",
        "as",     "at",      "be",    "because", "been",  "but",    "by",
        "can",    "cannot",  "could", "dear",    "did",   "do",     "does",
        "either", "else",    "ever",  "every",   "for",   "from",   "get",
        "got",    "had",     "has",   "have",    "he",    "her",    "hers",
        "him",    "his",     "how",   "however", "i",     "if",     "in",
        "into",   "is",      "it",    "its",     "just",  "least",  "let",
        "like",   "likely",  "may",   "me",      "might", "most",   "must",
        "my",     "neither", "no",    "nor",     "not",   "of",     "off",
        "often",  "on",      "only",  "or",      "other", "our",    "own",
        "rather", "said",    "say",   "says",    "she",   "should", "since",
        "so",     "some",    "than",  "that",    "the",   "their",  "them",
        "then",   "there",   "these", "they",    "this",  "tis",    "to",
        "too",    "twas",    "us",    "wants",   "was",   "we",     "were",
        "what",   "when",    "where", "which",   "while", "who",    "whom",
        "why",    "will",    "with",  "would",   "yet",   "you",    "your"};

    SDJournal journal;

    // extract tokens & numberic features
    string log_text;
    unordered_map<string, size_t> word_map;
    // first count, second cluster index
    UndirectedGraph<pair<unsigned, unsigned>, size_t> word_graph;
    deque<Message> messages;
    // start of log text

    bool to_wait = false;
    unsigned cycle_count = 0;
    unsigned retrieve_number = block_number * block_size;
    bool recluster = true;
    bool is_online = true;
    unordered_map<unsigned, double> prob_max;
    unordered_map<unsigned, ClusterInfo> clusters;

    while (true) {
        for (unsigned i = 0; i != retrieve_number; ++i) {
            if (journal >> log_text) {
                process_log(std::move(log_text), word_map, word_graph, messages,
                            stop_words, journal);

                // fade out phase
                ++cycle_count;
                if (cycle_count == block_size) {
                    cycle_count = 0;
                    if (is_online) {
                        recluster = true;
                        prob_max.clear();
                        clusters.clear();
                    }
                    fade_out(word_map, word_graph);
                }
                // end of log text
            } else {
                if (!journal.bad() && journal.eof())
                    i -= to_wait;
                else
                    throw; // todo
            }
        }

        map<unsigned, multimap<double, string>> probabilities;
        if (recluster) {
            Infomap infomapWrapper("--two-level -N2 --silent");
            for (auto ver_iter = word_graph.cbegin();
                 ver_iter != word_graph.cend(); ++ver_iter) {
                for (auto &edge : *ver_iter) {
                    infomapWrapper.addLink(ver_iter.index(), edge.head,
                                           edge.weight);
                }
            }
            infomapWrapper.run();
            for (infomap::LeafIterator leafIt(
                     &infomapWrapper.tree.getRootNode());
                 !leafIt.isEnd(); ++leafIt) {
                auto leaf_index = leafIt->originalLeafIndex;
                if (word_graph.valid(leaf_index)) {
                    auto &index = word_graph[leaf_index].value().second;
                    index = leafIt.moduleIndex();
                }
            }

            // kernel density estimation

            /* ************ */
            /* get clusters */

            for (auto word : word_map) {
                auto cluster_id = word_graph[word.second].value().second;
                auto iter = word_map.find(word.first);
                clusters[cluster_id].tokens.push_back(iter);
            }
            // for (auto &cluster : clusters) {
            //     cout << "cluster" << cluster.first << "\n";
            //     for (auto word_pair : cluster.second.tokens) {
            //         cout << word_pair->first << "\n";
            //     }
            //     cout << endl;
            // }

            for (auto &cluster : clusters) {
                cluster.second.gpnet = make_unique<GPNet<0>>(
                    cluster.second.tokens.size() + 3, 250, 50, 1, 1e-7);
            }

            typedef unordered_map<unsigned, unsigned> cluster_map;
            size_t progress_count = 0;
            for (auto &msg : messages) {
                // count cluster (via tokens)
                auto[input, msg_type] = construct_gpnet_input(
                    msg.map, word_map, word_graph, clusters);
                input.push_back(msg.count);
                input.push_back(msg.average);
                input.push_back(msg.variance);
                auto &gpnet = *clusters[msg_type].gpnet;
                gpnet.train(input);
                // ++progress_count;
                // cout << "\r" << progress_count << std::flush;
            }

            cout << endl;
            progress_count = 0;
            for (auto &msg : messages) {
                auto[input, msg_type] = construct_gpnet_input(
                    msg.map, word_map, word_graph, clusters);
                input.push_back(msg.count);
                input.push_back(msg.average);
                input.push_back(msg.variance);
                auto &gpnet = *clusters[msg_type].gpnet;
                double prob = gpnet.predict(input);
                if (is_online) {
                    if (!isnan(prob) && !isinf(prob) &&
                        prob_max[msg_type] < prob)
                        prob_max[msg_type] = prob;
                } else
                    probabilities[msg_type].insert({prob, msg.cursor});
                // ++progress_count;
                // cout << "\r" << progress_count << std::flush;
                // cout << prob << "\n";
            }
            cout << endl;
            if (is_online) {
                recluster = false;
                messages.erase(messages.cbegin(),
                               messages.cbegin() + block_size);
            }
        }

        // get anomalies
        if (is_online) {
            retrieve_number = 1;
            auto &msg = messages.back();
            auto[input, msg_type] =
                construct_gpnet_input(msg.map, word_map, word_graph, clusters);
            if (prob_max[msg_type]) {
                input.push_back(msg.count);
                input.push_back(msg.average);
                input.push_back(msg.variance);
                auto &gpnet = *clusters[msg_type].gpnet;
                double prob = gpnet.predict(input);
                if (prob < pow(prob_max[msg_type], -2)) {
                    cout << msg.cursor << endl;
                }
            }
        } else {
            vector<string> anomalies;
            for (auto &cluster : probabilities) {
                // for (auto &msg : cluster.second)
                //     cout << msg.first << "\n" << msg.second << "\n";
                // cout << endl;
                auto &cluster_logs = cluster.second;

                double max_prob = 0;
                for (auto riter = cluster_logs.crbegin();
                     riter != cluster_logs.crend(); ++riter) {
                    double prob = riter->first;
                    if (!isnan(prob) && !isinf(prob)) {
                        max_prob = prob;
                        break;
                    }
                }
                if (max_prob) {
                    for (auto iter = cluster_logs.cbegin();
                         iter != cluster_logs.cend(); ++iter) {
                        double prob = iter->first;
                        if (isnan(prob))
                            continue;
                        if (prob >= pow(max_prob, -2))
                            break;

                        // cout << prob << "\n";
                        // cout << max_prob << "\n\n";
                        anomalies.push_back(iter->second);
                    }
                }
            }
            cout << endl;
            cout << anomalies.size() << "\n";
            for (auto &cursor : anomalies)
                cout << cursor << "\n";

            cout << endl;
            break;
        }
    }

    return 0;
}

namespace {

Extractor::result_type Extractor::extract() {
    for (size_t pos = m_word_beg; pos != m_text.size(); ++pos) {
        char ch = m_text[pos];
        switch (ch) {
        case '\n':
            if (m_state) {
                auto result = construct_result(pos);
                m_word_beg = pos;
                m_state = s_space;
                return result;
            }
            m_word_beg = pos + 1;
            return {newline, string()};
            break;
        case ' ':
            if (m_state) {
                auto result = construct_result(pos);
                m_word_beg = pos + 1;
                m_state = s_space;
                return result;
            }
            break;
        default:
            switch (m_state) {
            case s_space:
                m_word_beg = pos;
                if (ch == '0')
                    m_state = s_zero;
                else if (is_base10(ch))
                    m_state = s_b10;
                else if (is_base16(ch))
                    m_state = s_b16;
                else
                    m_state = s_word;
                break;
            case s_zero:
                if (ch == 'x')
                    m_state = s_zerox;
                else if (is_base10(ch))
                    m_state = s_b10;
                else if (is_base16(ch))
                    m_state = s_b16;
                else {
                    auto result = construct_result(pos);
                    m_word_beg = pos;
                    m_state = s_space;
                    return result;
                }
                break;
            case s_zerox:
                if (is_base16(ch))
                    m_state = s_b16_c;
                else
                    m_state = s_word;
                break;
            case s_b16_c:
            case s_b16:
                if (!is_base16(ch))
                    m_state = s_word;
                break;
            case s_b10:
                if (!is_base10(ch)) {
                    if (is_base16(ch))
                        m_state = s_b16;
                    // else
                    //     m_state = s_word;
                    else {
                        auto result = construct_result(pos);
                        m_word_beg = pos;
                        m_state = s_space;
                        return result;
                    }
                }
                break;
            case s_word:
                if (is_base10(ch)) {
                    auto result = construct_result(pos);
                    m_word_beg = pos;
                    m_state = s_space;
                    return result;
                }
                break;
            default:
                break;
            }
            break;
        }
    }
    if (m_state) {
        auto sz = m_text.size();
        auto result = construct_result(sz);
        m_word_beg = sz;
        m_state = s_space;
        return result;
    }
    return {eof, string()};
}
void process_log(std::string log_text,
                 std::unordered_map<std::string, size_t> &word_map,
                 LogAnomaly::UndirectedGraph<std::pair<unsigned, unsigned>,
                                             size_t> &word_graph,
                 std::deque<Message> &messages,
                 std::unordered_set<std::string> const &stop_words,
                 LogAnomaly::SDJournal &journal) {
    using LogAnomaly::PorterStemmer::stem;
    using std::accumulate;
    using std::make_pair;
    using std::vector;
    // store message related
    messages.push_back(Message{
        .cursor = journal.tellg(), .count = 0, .average = 0, .variance = 0});
    Message &record = messages.back();
    vector<unsigned long long> numbers;

    // only words and numbers
    // in lower case
    transform(log_text.begin(), log_text.end(), log_text.begin(),
              [](char ch) -> char {
                  //   if (!std::isalnum(ch) && ch != '-' && ch != '\n')
                  if (!std::isalnum(ch) && ch != '\n')
                      return ' ';
                  return std::tolower(ch);
              });

    Extractor extraction(std::move(log_text), 8);
    bool predecessor = false;
    size_t pre_index;
    // word extraction start
    while (true) {
        typedef Extractor::Type Type;

        auto[t, word] = extraction.extract();
        switch (t) {
        // numeric features
        case Type::base16:
            if (word.size() < 6)
                break;
        case Type::c_base16:
            try {
                numbers.push_back(stoull(word, nullptr, 16));
            } catch (std::out_of_range &e) {
            }
            break;
        case Type::base10:
            numbers.push_back(stoull(word));
            break;

        // tokens
        case Type::word: {
            if (stop_words.count(word) == 1)
                break;

            // constexpr unsigned WORD_LENGTH_MAX = 15;
            // constexpr unsigned WORD_LENGTH_MIN = 2;
            word = stem(word);
            // if (word.size() > WORD_LENGTH_MAX || word.size() <
            // WORD_LENGTH_MIN)
            //     break;
            auto[iter, success] = word_map.insert({word, 0});
            if (success) {
                auto ver_iter = word_graph.insert_vertex(
                    make_pair<unsigned, unsigned>(0, 0));
                iter->second = ver_iter.index();
            }

            ++word_graph[iter->second].value().first;
            ++record.map[word];
            if (predecessor) {
                for (auto &edge : word_graph[iter->second]) {
                    if (edge.head == pre_index) {
                        ++edge.weight;
                        goto end_predecessor;
                    }
                }
                word_graph.insert_edge(iter->second, pre_index, 1);
            }
        end_predecessor:
            predecessor = true;
            pre_index = iter->second;
            break;
        }
        case Type::eof:
            goto extraction_done;
        case Type::newline:
            predecessor = false;
            break;
        }
    }
    // word extraction end
extraction_done:
    // calculating the numeric features
    if (numbers.size() > 0) {
        record.count = numbers.size();
        record.average = accumulate(numbers.cbegin(), numbers.cend(), 0ULL) /
                         static_cast<double>(record.count);
        if (numbers.size() > 1) {
            double sum = 0;
            for (auto &num : numbers) {
                sum += (num - record.average) * (num - record.average);
            }
            record.variance = sum / (record.count - 1);
        }
    }
}
void fade_out(std::unordered_map<std::string, size_t> &word_map,
              LogAnomaly::UndirectedGraph<std::pair<unsigned, unsigned>, size_t>
                  &word_graph) {
    for (auto ver_iter = word_graph.begin(); ver_iter != word_graph.end();) {
        // vertex fade
        ver_iter->value().first /= vertex_fade_coeff;
        if (ver_iter->value().first < vertex_fade_threshold)
            ver_iter = word_graph.erase_vertex(ver_iter);
        else {
            for (auto edge_iter = ver_iter->begin(),
                      pre_iter = ver_iter->before_begin();
                 edge_iter != ver_iter->end();) {
                edge_iter->weight /= edge_fade_coeff;
                if (edge_iter->weight < edge_fade_threshold) {
                    edge_iter = word_graph.erase_after_edge(ver_iter, pre_iter);
                } else {
                    ++edge_iter;
                    ++pre_iter;
                }
            }
            ++ver_iter;
        }
    }
    for (auto iter = word_map.cbegin(); iter != word_map.cend();) {
        if (!word_graph.valid(iter->second))
            iter = word_map.erase(iter);
        else
            ++iter;
    }
}
std::pair<std::vector<double>, unsigned> construct_gpnet_input(
    std::unordered_map<std::string, size_t> const &map,
    std::unordered_map<std::string, size_t> &word_map,
    LogAnomaly::UndirectedGraph<std::pair<unsigned, unsigned>, size_t> const
        &word_graph,
    std::unordered_map<unsigned, ClusterInfo> &clusters) {
    typedef std::unordered_map<unsigned, unsigned> cluster_map;
    using std::vector;
    cluster_map cluster_count;
    for (auto &word : map) {
        auto token_index = word_map[word.first];
        if (word_graph.valid(token_index)) {
            auto cluster_index = word_graph[token_index].value().second;
            ++cluster_count[cluster_index];
        }
    }

    // find max cluster
    size_t max = cluster_count[0];
    auto miter = cluster_count.cbegin();
    for (auto iter = cluster_count.cbegin(); iter != cluster_count.cend();
         ++iter) {
        if (iter->second > max) {
            max = iter->second;
            miter = iter;
        }
    }

    // construct input
    auto msg_type = miter->first;
    auto &cluster_tokens = clusters[msg_type].tokens;
    // GPSOINN input
    vector<double> input;
    for (auto token : cluster_tokens) {
        auto riter = map.find(token->first);

        if (riter == map.end()) {
            input.push_back(0);
        } else {
            input.push_back(riter->second);
        }
    }
    return {input, msg_type};
    // input.push_back(msg.count);
    // input.push_back(msg.average);
    // input.push_back(msg.variance);
}
} // namespace
