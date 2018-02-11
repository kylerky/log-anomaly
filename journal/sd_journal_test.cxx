#include "sd_journal.hxx"
#include <iostream>

int main() {
    using namespace LogAnomaly;
    using std::string;

    SDJournal journal;
    string str;
    // for (unsigned i = 0; i != 100; ++i) {
    while (journal >> str) {
        // std::cout << str << std::endl;
    }
    journal.clear();
    journal.set_timeout(30e6);

    while (!journal.fail()) {
        while (journal >> str) {
            std::cout << str << std::endl;
        }
        if (journal.invalidated() || journal.eof())
            journal.clear();
    }
    std::cout << journal.good() << std::endl;
    std::cout << journal.bad() << std::endl;
    std::cout << journal.eof() << std::endl;

    return 0;
    // }
}
