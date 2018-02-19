#include "multiset.hxx"
#include <gtest/gtest.h>
#include <iostream>
#include <random>
#include <unordered_set>

using namespace LogAnomaly;

#define SEQ_ELEM_CNT 5000
#define RAND_ERASE_TRY_CNT 10
#define RAND_ERASE_CNT 1000
#define FREE_LIST_DELTA 0.2

class SequenceSet : public testing::Test {
  protected:
    virtual void SetUp() {
        for (int i = 0; i != SEQ_ELEM_CNT; ++i)
            set.insert(i);
    }
    virtual void TearDown() {}
    multiset<int> set;
};

TEST(iterator, value_init) {
    multiset<int>::iterator iter1;
    multiset<int>::iterator iter2;
    multiset<int>::const_iterator citer1;
    multiset<int>::const_iterator citer2;
    ASSERT_EQ(iter1, iter2);
    ASSERT_EQ(citer1, citer2);
}

TEST_F(SequenceSet, op) {
    for (int i = 0; i != SEQ_ELEM_CNT; ++i) {
        ASSERT_EQ(i, set[i]);
    }
}

TEST_F(SequenceSet, iter_increment) {
    {
        int i = 0;
        for (auto iter = set.begin(); iter != set.end(); ++iter) {
            ASSERT_EQ(*iter, i);
            ++i;
        }
    }
    {
        int i = 0;
        for (auto iter = set.cbegin(); iter != set.cend(); ++iter) {
            ASSERT_EQ(*iter, i);
            ++i;
        }
    }
}

TEST_F(SequenceSet, iter_copy) {
    auto iter = set.begin();
    for (int i = 0; i != SEQ_ELEM_CNT; ++i) {
        auto iter2 = iter;
        auto iter3 = iter;
        ASSERT_EQ(*iter2, *iter);
        ASSERT_EQ(*iter3, *iter);
        for (unsigned j = i; j != SEQ_ELEM_CNT; ++j) {
            ASSERT_EQ(*iter2, *iter3);
            ++iter2;
            ++iter3;
        }
        ++iter;
    }
}

TEST(erase, random) {
    std::random_device rand_dev;
    for (unsigned i = 0; i != RAND_ERASE_TRY_CNT; ++i) {
        multiset<int> set;
        for (int i = 0; i != SEQ_ELEM_CNT; ++i)
            set.insert(i);
        std::random_device::result_type seed = rand_dev();
        std::default_random_engine eng(seed);
        std::uniform_int_distribution<int> distribution(0, SEQ_ELEM_CNT - 1);

        std::unordered_set<int> drops;
        for (unsigned j = 0; j != RAND_ERASE_CNT; ++j) {
            int pos = distribution(eng);
            set.erase(pos);
            drops.insert(pos);
        }

        int index = 0;
        for (auto iter = set.cbegin(); iter != set.cend(); ++iter) {
            while (drops.count(index) == 1) {
                EXPECT_EQ(set.valid(index), false);
                ++index;
            }
            EXPECT_EQ(set.valid(index), true);
            EXPECT_EQ(*iter, index) << "seed: " << seed << "\nat try " << i;
            ++index;
        }
    }
}

TEST(freelist, part) {
    std::random_device rand_dev;
    for (unsigned i = 0; i != RAND_ERASE_TRY_CNT; ++i) {
        multiset<int> set;
        std::unordered_set<int> std_set;
        for (int i = 0; i != SEQ_ELEM_CNT; ++i) {
            set.insert(i);
            std_set.insert(i);
        }
        ASSERT_EQ(set.size(), std_set.size())
            << "unmatched size after first insertion";
        std::random_device::result_type seed = rand_dev();
        std::default_random_engine eng(seed);
        std::uniform_int_distribution<int> distribution(0, SEQ_ELEM_CNT - 1);

        for (unsigned j = 0; j != RAND_ERASE_CNT; ++j) {
            int pos = distribution(eng);
            set.erase(pos);
            std_set.erase(pos);
        }
        ASSERT_EQ(set.size(), std_set.size()) << "unmatched size after erase";
        for (unsigned j = SEQ_ELEM_CNT;
             j != SEQ_ELEM_CNT + FREE_LIST_DELTA * RAND_ERASE_CNT; ++j) {
            set.insert(j);
            std_set.insert(j);
        }
        ASSERT_EQ(set.size(), std_set.size())
            << "unmatched size after second insertion";
        for (auto elem : set)
            EXPECT_TRUE(std_set.count(elem));
    }
}

TEST(freelist, overwhelm) {
    std::random_device rand_dev;
    for (unsigned i = 0; i != RAND_ERASE_TRY_CNT; ++i) {
        multiset<int> set;
        std::unordered_set<int> std_set;
        for (int i = 0; i != SEQ_ELEM_CNT; ++i) {
            set.insert(i);
            std_set.insert(i);
        }
        ASSERT_EQ(set.size(), std_set.size())
            << "unmatched size after first insertion";
        std::random_device::result_type seed = rand_dev();
        std::default_random_engine eng(seed);
        std::uniform_int_distribution<int> distribution(0, SEQ_ELEM_CNT - 1);

        for (unsigned j = 0; j != RAND_ERASE_CNT; ++j) {
            int pos = distribution(eng);
            set.erase(pos);
            std_set.erase(pos);
        }
        ASSERT_EQ(set.size(), std_set.size()) << "unmatched size after erase";
        for (unsigned j = SEQ_ELEM_CNT;
             j != SEQ_ELEM_CNT + 2 * FREE_LIST_DELTA * RAND_ERASE_CNT; ++j) {
            set.insert(j);
            std_set.insert(j);
        }
        ASSERT_EQ(set.size(), std_set.size())
            << "unmatched size after second insertion";
        for (auto elem : set)
            EXPECT_TRUE(std_set.count(elem));
    }
}
