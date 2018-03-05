#include "journal/sd_journal.hxx"

#include "tclap/CmdLine.h"

#include <ncurses.h>

#include <algorithm>
#include <cerrno>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace LogAnomaly {

namespace Viewer {

class NWindow {
    friend class Ncurses;
    using string = std::string;
    using ch_string = std::basic_string<chtype>;

  public:
    typedef std::pair<int, int> dimension;
    NWindow() : m_win(nullptr){};
    NWindow(int nlines, int ncols, int begin_y, int begin_x) {
        using std::runtime_error;

        m_win = newwin(nlines, ncols, begin_y, begin_x);
        if (!m_win)
            throw runtime_error("Failed to create window");
    }
    NWindow(const NWindow &other) = delete;
    NWindow(NWindow &&other) {
        using std::swap;

        swap(m_win, other.m_win);
    }
    ~NWindow() { delwin(m_win); };

    friend void swap(NWindow &left, NWindow &right) {
        using std::swap;
        swap(left.m_win, right.m_win);
    }
    NWindow &operator=(NWindow other) {
        using std::swap;
        swap(*this, other);
        return *this;
    }
    friend NWindow &operator<<(NWindow &window, std::string str) {
        waddstr(window.m_win, str.c_str());
        return window;
    }
    friend NWindow &operator<<(NWindow &window, int integer) {
        wprintw(window.m_win, "%d", integer);
        return window;
    }
    friend NWindow &operator<<(NWindow &window, unsigned integer) {
        wprintw(window.m_win, "%u", integer);
        return window;
    }
    friend NWindow &operator<<(NWindow &window, unsigned long integer) {
        wprintw(window.m_win, "%lu", integer);
        return window;
    }
    friend NWindow &operator<<(NWindow &window, long integer) {
        wprintw(window.m_win, "%ld", integer);
        return window;
    }
    friend NWindow &operator<<(NWindow &window,
                               NWindow &(*manipulator)(NWindow &)) {
        return manipulator(window);
    }
    NWindow &refresh() {
        using std::runtime_error;

        int result = wrefresh(m_win);
        if (result == ERR)
            throw runtime_error("Failed to refresh window");
        return *this;
    }
    NWindow &mvAdd(int x, int y, string str) {
        using std::runtime_error;
        int result = mvwaddstr(m_win, y, x, str.c_str());
        if (result == ERR)
            throw runtime_error("Failed to add to window");

        return *this;
    }
    NWindow &mvAdd(int x, int y, string str, int n) {
        using std::runtime_error;

        int result = mvwaddnstr(m_win, y, x, str.c_str(), n);
        if (result == ERR)
            throw runtime_error("Failed to add to window");
        return *this;
    }
    NWindow &mvAddCh(int x, int y, ch_string str) {
        using std::runtime_error;
        int result = mvwaddchstr(m_win, y, x, str.c_str());
        if (result == ERR)
            throw runtime_error("Failed to add to window");

        return *this;
    }
    NWindow &mvAddCh(int x, int y, ch_string str, int n) {
        using std::runtime_error;
        int result = mvwaddchnstr(m_win, y, x, str.c_str(), n);
        if (result == ERR)
            throw runtime_error("Failed to add to window");

        return *this;
    }
    NWindow &clearToBot() {
        wclrtobot(m_win);
        return *this;
    }
    NWindow &clearToBot(int x, int y) {
        wmove(m_win, y, x);
        wclrtobot(m_win);
        return *this;
    }
    NWindow &clearToEol() {
        wclrtoeol(m_win);
        return *this;
    }
    NWindow &clearToEol(int x, int y) {
        wmove(m_win, y, x);
        wclrtoeol(m_win);
        return *this;
    }
    NWindow &moveCursor(int x, int y) {
        wmove(m_win, y, x);
        return *this;
    }
    int getChar() { return wgetch(m_win); }
    void scrollOK(bool ok) { scrollok(m_win, ok); }

    dimension getMaxXY() {
        int x, y;
        getmaxyx(m_win, y, x);
        return {x, y};
    }
    dimension getXY() {
        int x, y;
        getyx(m_win, y, x);
        return {x, y};
    }

  protected:
    NWindow(WINDOW *win) : m_win(win) {}
    WINDOW *m_win;
} cscr;

inline NWindow &flush(NWindow &win) { return win.refresh(); }

class Ncurses {
  public:
    Ncurses() {
        initscr();              // start curses mode
        cbreak();                  // disable line buffer
        noecho();               // hide echo characters
        keypad(stdscr, TRUE);   // get keys
        cscr = NWindow(stdscr); // construct cscr
        cscr << flush;
    }
    ~Ncurses() { endwin(); }
};

class Viewer {
    using SDJournal = LogAnomaly::SDJournal;
    using string = std::string;
    using ch_string = std::basic_string<chtype>;
    using StringVector = std::vector<string>;

  public:
    Viewer(const string &path);

    ~Viewer() {}

    bool scrollY() { return scrollY(1); }
    bool scrollY(int n) {
        size_t result = m_beg_y + n;
        if ((n < 0 && m_beg_y < result) || (n > 0 && m_beg_y > result))
            return false;
        int maxY = m_buffer_view.getMaxXY().second;
        size_t bot = result + maxY;
        if (result > bot || bot > m_buffer->size())
            return false;

        m_beg_y = result;
        paint();
        return true;
    }
    void scrollX(int n);
    void scrollXbeg() {
        scrollX(-m_beg_x);
        auto[x, y] = m_buffer_view.getXY();
        moveCursor(-x);
    }
    void moveCursor(int delta_x, int delta_y = 0);
    void paint();
    void switchInfo();
    void updateStatus() {
        auto[x, y] = m_buffer_view.getXY();
        m_status_bar.clearToEol(0, 0)
            << x + m_beg_x << ", " << y + m_beg_y << flush;
        m_buffer_view.moveCursor(x, y) << flush;
    }
    void switchLog() {
        m_beg_x = m_pre_beg_x;
        m_beg_y = m_pre_beg_y;

        m_log_lines.clear();
        m_buffer = &m_lines;
        m_cluster = true;
        paint();
        m_buffer_view.moveCursor(m_pre_x, m_pre_y) << flush;
        updateStatus();
    }
    void paint(int line);

  private:
    StringVector *m_buffer;
    StringVector m_lines;
    StringVector m_log_lines;
    StringVector m_cursors;
    std::vector<size_t> m_indices;

    NWindow m_buffer_view;
    NWindow m_status_bar;

    bool m_cluster;

    unsigned fetch(unsigned count);
    SDJournal m_journal;
    std::ifstream m_file;

    size_t m_beg_y = 0;
    size_t m_beg_x = 0;
    int m_pre_x;
    int m_pre_y;
    size_t m_pre_beg_y = 0;
    size_t m_pre_beg_x = 0;
    size_t m_log_index = 0;
};

} // namespace Viewer
} // namespace LogAnomaly
