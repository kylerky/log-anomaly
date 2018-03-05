#include "check.hxx"

int main(int argc, char *argv[]) {
    using std::string;
    // arguments parsing
    string path;
    try {
        using namespace TCLAP;

        CmdLine cmd("Log viewer", ' ', "0.1");

        UnlabeledValueArg<string> filenameArg(
            "path", "The path to the cursor file", true, "", "path", cmd);
        cmd.parse(argc, argv);
        path = filenameArg.getValue();
    } catch (TCLAP::ArgException &e) {
        using std::cerr;
        using std::endl;
        cerr << "Failed to parse the arguments: " << e.error() << " for arg "
             << e.argId() << endl;
        return 1;
    }

    {
        using namespace LogAnomaly::Viewer;
        Ncurses curses_scr;
        Viewer viewer(path);
        while (true) {
            int ch = cscr.getChar();

            bool quit = false;
            switch (ch) {
            case 'j':
                viewer.moveCursor(0, 1);
                break;
            case 'k':
                viewer.moveCursor(0, -1);
                break;
            case 'h':
                viewer.moveCursor(-1);
                break;
            case 'l':
                viewer.moveCursor(1);
                break;
            case 'q':
                quit = true;
                break;
            case '0':
                viewer.scrollXbeg();
                break;
            case ']':
                viewer.switchInfo();
                break;
            case '[':
                viewer.switchLog();
                break;
            default:
                break;
            }
            if (quit)
                break;
        }
    }
    return 0;
}

namespace LogAnomaly {

namespace Viewer {

Viewer::Viewer(const string &path)
    : m_file(path), m_buffer(&m_lines), m_cluster(true) {
    using std::runtime_error;

    if (m_file.bad())
        throw runtime_error("Failed to open the file");
    auto[x, y] = cscr.getMaxXY();
    if (y < 2)
        throw runtime_error("Too small");
    m_buffer_view = NWindow(y - 1, x, 0, 0);
    m_status_bar = NWindow(1, x, y - 1, 0);
    m_buffer_view.scrollOK(false);
    m_status_bar.scrollOK(false);

    updateStatus();
    fetch(y);
    paint();
}

void Viewer::paint(int line) {
    auto[origX, origY] = m_buffer_view.getXY();
    auto[maxX, maxY] = m_buffer_view.getMaxXY();

    ch_string format_str;
    for (auto &ch : (*m_buffer)[line + m_beg_y])
        format_str.push_back(ch);

    m_buffer_view.mvAddCh(0, line, format_str);
    if (format_str.size() < maxX) {
        m_buffer_view.clearToEol(format_str.size(), line);
    }
    m_buffer_view.moveCursor(origX, origY) << flush;
}
void Viewer::moveCursor(int delta_x, int delta_y) {
    auto[orig_x, orig_y] = m_buffer_view.getXY();
    auto[max_x, max_y] = m_buffer_view.getMaxXY();
    int dest_x = orig_x + delta_x;
    int dest_y = orig_y + delta_y;
    bool y_update = false;

    if (dest_y >= max_y) {
        y_update = scrollY(dest_y - max_y + 1);
        if (!y_update) {
            fetch(max_y);
            y_update = scrollY(dest_y - max_y + 1);
        }
        dest_y = max_y - 1;
    } else if (dest_y < 0) {
        y_update = scrollY(dest_y);
        dest_y = 0;
    }
    if (dest_x >= max_x) {
        scrollX(dest_x - max_x + 1);
        dest_x = max_x - 1;
    } else if (dest_x < 0) {
        scrollX(dest_x);
        dest_x = 0;
    }

    auto line_length = (*m_buffer)[m_beg_y + dest_y].size();
    if (delta_y) {
        size_t index = m_indices[m_beg_y + dest_y];
        if (m_beg_x) {
            if (!m_cluster || index != m_log_index) {
                m_beg_x = 0;
                scrollX(0);
            } else {
                dest_x = 4;
                scrollX(line_length - m_beg_x - dest_x);
            }
        }
        m_log_index = index;
    }
    if (m_beg_x + dest_x >= line_length) {
        dest_x = line_length - m_beg_x - 1;
    }
    m_buffer_view.moveCursor(dest_x, dest_y) << flush;
    updateStatus();
}
void Viewer::paint() {
    using std::runtime_error;
    // m_buffer_view << "hi" << flush;
    auto[origX, origY] = m_buffer_view.getXY();
    auto[maxX, maxY] = m_buffer_view.getMaxXY();
    size_t i = 0;
    for (; i != maxY && i + m_beg_y != m_buffer->size(); ++i) {
        // m_buffer_view.mvAdd(0, i, m_lines[i + m_beg_y], maxX);
        ch_string format_str;
        for (auto &ch : (*m_buffer)[i + m_beg_y])
            format_str.push_back(ch);

        m_buffer_view.mvAddCh(0, i, format_str);
        if (format_str.size() < maxX) {
            m_buffer_view.clearToEol(format_str.size(), i);
        }
    }
    if (i != maxY) {
        string &last = (*m_buffer)[i - 1];
        if (last.size() > maxX) {
            if (i == maxY - 1)
                return;
            m_buffer_view.clearToBot(0, i + 1);
        } else
            m_buffer_view.clearToBot(last.size(), maxY - 1);
    }
    m_buffer_view.moveCursor(origX, origY) << flush;
}

unsigned Viewer::fetch(unsigned log_cnt) {
    using LogAnomaly::SDJournal;
    using LogAnomaly::Viewer::flush;
    using std::asctime;
    using std::find;
    using std::localtime;
    using std::ostringstream;
    using std::stoull;
    using std::string;
    typedef std::vector<string::const_iterator> IterVector;

    int cnt = 0;
    for (; cnt != log_cnt; ++cnt) {
        string cursor;
        // double prob;
        if (m_file >> cursor) {
            SDJournal::LogMap log_map;

            // m_file >> cursor >> prob;
            m_cursors.push_back(cursor);

            m_journal.seekg(cursor);
            // m_journal >> text;
            m_journal >> log_map;
            string text = log_map["MESSAGE"].substr(8);

            ostringstream info;

            time_t time_stamp;
            if (log_map.count("_SOURCE_REALTIME_TIMESTAMP")) {
                auto &real_timestamp = log_map["_SOURCE_REALTIME_TIMESTAMP"];
                time_stamp =
                    static_cast<time_t>(stoull(real_timestamp.substr(27)) /
                                        1e6); // 1e6 for micro second to second
            } else {
                auto &real_timestamp = log_map["__REALTIME_TIMESTAMP"];
                time_stamp =
                    static_cast<time_t>(m_journal.get_realtime_usec() /
                                        1e6); // 1e6 for micro second to second
            }
            char time_buffer[16];
            strftime(time_buffer, 16, "%b %e %T", localtime(&time_stamp));
            info << time_buffer << " ";
            info << log_map["_HOSTNAME"].substr(10) << " "
                 << log_map["_COMM"].substr(6) << "["
                 << log_map["_PID"].substr(5) << "]: ";
            {
                string line = info.str();
                bool first = true;
                for (size_t i = 0; i != text.size();) {
                    auto pos = text.find('\n', i);
                    if (pos == string::npos) {
                        line.append(text.substr(i));
                        m_lines.push_back(line);
                        m_indices.push_back(m_cursors.size() - 1);
                        break;
                    } else {
                        line.append(text.substr(i, pos + 1));
                        m_lines.push_back(line);
                        m_indices.push_back(m_cursors.size() - 1);
                        i = pos + 1;
                    }
                    if (first) {
                        first = false;
                        line = string(info.str().size(), ' ');
                    }
                }
            }
        } else
            break;
    }
    return cnt;
}
void Viewer::scrollX(int n) {
    auto[curr_x, curr_y] = m_buffer_view.getXY();
    auto[max_x, max_y] = m_buffer_view.getMaxXY();

    size_t beg = m_beg_y + curr_y;
    size_t end = m_beg_y + curr_y + 1;
    size_t index = m_indices[beg];
    m_log_index = index;
    if (m_cluster) {
        while (beg != 0) {
            if (m_indices[beg - 1] == index)
                --beg;
            else
                break;
        }
        while (end != m_indices.size()) {
            if (m_indices[end] == index)
                ++end;
            else
                break;
        }
    }

    size_t max_length = 0;
    for (size_t i = beg; i != end; ++i) {
        size_t length = (*m_buffer)[i].size();
        if (length > max_length)
            max_length = length;
    }

    size_t dest_x = m_beg_x + n;
    if ((n > 0 && dest_x < m_beg_x) || (n < 0 && dest_x > m_beg_x) ||
        (max_length > max_x && dest_x > max_length - max_x))
        return;

    size_t refresh_beg = beg > m_beg_y ? beg - m_beg_y : 0;
    size_t refresh_end =
        end < m_beg_y + max_y ? end - m_beg_y : m_beg_y + max_y - 1;
    for (size_t i = refresh_beg; i != refresh_end; ++i) {
        // m_buffer_view.mvAdd(0, refresh_beg,
        //                     m_lines[refresh_beg + m_beg_y].substr(dest_x),
        //                     max_x);
        size_t line_n = i + m_beg_y;
        string &str = (*m_buffer)[line_n];
        if (str.size() <= dest_x) {
            m_buffer_view.clearToEol(0, line_n);
            continue;
        }

        ch_string format_str;
        for (auto iter = str.cbegin() + dest_x;
             iter != str.cend() && format_str.size() < max_x; ++iter)
            format_str.push_back(*iter);
        // for (auto &ch : m_lines[refresh_beg + m_beg_y])
        //     format_str.push_back(ch);

        m_buffer_view.mvAddCh(0, line_n, format_str);

        if (format_str.size() < max_x)
            m_buffer_view.clearToEol(format_str.size(), line_n);
    }
    m_beg_x = dest_x;
    m_buffer_view << flush;
}
void Viewer::switchInfo() {
    auto[curr_x, curr_y] = m_buffer_view.getXY();

    m_pre_beg_x = m_beg_x;
    m_pre_beg_y = m_beg_y;

    m_pre_x = curr_x;
    m_pre_y = curr_y;

    m_beg_x = 0;
    m_beg_y = 0;

    auto &cursor = m_cursors[m_indices[m_pre_beg_y + curr_y]];
    m_journal.seekg(cursor);
    SDJournal::LogMap log_map;
    m_journal >> log_map;
    size_t field_max = 0;
    for (auto &elem : log_map) {
        auto length = elem.first.size();
        if (length > field_max)
            field_max = length;
    }
    for (auto &elem : log_map) {
        string text = elem.second.substr(elem.first.size() + 1);
        string field = elem.first;
        field.resize(field_max + 2, ' ');
        string line = field;
        bool first = true;
        for (size_t i = 0; i != text.size();) {
            auto pos = text.find('\n', i);
            if (pos == string::npos) {
                line.append(text.substr(i));
                m_log_lines.push_back(line);
                break;
            } else {
                line.append(text.substr(i, pos + 1));
                m_log_lines.push_back(line);
                i = pos + 1;
            }
            if (first) {
                first = false;
                line = string(field.size(), ' ');
            }
        }
    }
    m_cluster = false;
    m_buffer = &m_log_lines;
    paint();
    m_buffer_view.moveCursor(0, 0) << flush;
    updateStatus();
}

} // namespace Viewer

} // namespace LogAnomaly
