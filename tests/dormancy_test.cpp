

#include <gtest/gtest.h>

#include "../src/dormancy.h"


// function to test the tests
TEST(test_sleepy_create_mock_nodes, sleepy_create_mock_nodes) {
  tsk_node_table_t nodes;
  tsk_node_table_init(&nodes, 0);
  sleepy_create_mock_nodes(&nodes, 11);
  EXPECT_EQ(nodes.num_rows, 11);
  tsk_node_table_free(&nodes);
}

TEST(test_sleepy_init_tables, sleepy_init_tables) {
  tsk_table_collection_t tables;
  tsk_table_collection_init(&tables, 0);
  int N = 3, m = 1, gen = 1;
  sleepy_init_tables(&tables, N, m, gen);
  EXPECT_EQ(tables.nodes.num_rows, 2*N*m);
  tsk_table_collection_free(&tables);


  tsk_table_collection_t tables2;
  tsk_table_collection_init(&tables2, 0);
  N = 12, m = 5, gen = 45;
  sleepy_init_tables(&tables2, N, m, gen);
  EXPECT_EQ(tables2.nodes.num_rows, 2*N*m);
  tsk_table_collection_free(&tables2);
}

TEST(test_sleepy_reverse_time, sleepy_reverse_time) {
  tsk_table_collection_t tables;  
  tsk_table_collection_init(&tables, 0);
  sleepy_init_tables(&tables, 3, 3, 0);

  int begin = 2;
  int end = 10;
  int length = end - begin;
  std::vector<int> times;

  for (int i=begin; i<end; i++) {
      times.emplace_back(tables.nodes.time[i]);
  }

  sleepy_reverse_time(tables.nodes, begin, end);
  std::vector<int> times_after_reversing;
  for (int i=begin; i<end; i++) {
      times_after_reversing.emplace_back(tables.nodes.time[i]);
  }

  for (int i = 0; i < times.size(); ++i) {
    ASSERT_NE(times_after_reversing[i], times[i]) << "After first reversing equal at " << times[i];
  }

  sleepy_reverse_time(tables.nodes, begin, end);
  ASSERT_EQ(times.size(), length) << "Unequal size " << times.size() << " " << length;
  for (int i = 0; i < times.size(); ++i) {
    ASSERT_EQ(times[i], tables.nodes.time[begin+i]) << "After second reversing unequal at " << times[i] << " " << tables.nodes.time[i];
  }

}

TEST(test_sleepy_reverse_time2, sleepy_reverse_time) {
  tsk_table_collection_t tables;  
  tsk_table_collection_init(&tables, 0);
  sleepy_init_tables(&tables, 3, 3, 0);

  int begin = 0;
  int end = 10;
  int length = end - begin;
  std::vector<int> times;

  for (int i=begin; i<end; i++) {
      times.emplace_back(tables.nodes.time[i]);
  }

  sleepy_reverse_time(tables.nodes, begin, end);
  std::vector<int> times_after_reversing;
  for (int i=begin; i<end; i++) {
      times_after_reversing.emplace_back(tables.nodes.time[i]);
  }

  for (int i = 0; i < times.size(); ++i) {
    ASSERT_NE(times_after_reversing[i], times[i]) << "After first reversing equal at " << times[i];
  }

  sleepy_reverse_time(tables.nodes, begin, end);
  ASSERT_EQ(times.size(), length) << "Unequal size " << times.size() << " " << length;
  for (int i = 0; i < times.size(); ++i) {
    ASSERT_EQ(times[i], tables.nodes.time[begin+i]) << "After second reversing unequal at " << times[i] << " " << tables.nodes.time[i];
  }

}



TEST(test_sleepy_num_non_sample_node, sleepy_num_non_sample_node) {
  int gen = 0;
  tsk_table_collection_t tables;
  tsk_table_collection_init(&tables, 0);
  sleepy_init_tables(&tables, 3, 2, gen);
  tables.nodes.flags[0] = 0;
  tables.nodes.flags[4] = 0;
  tsk_id_t num_non_sample = sleepy_num_non_sample_node(tables.nodes);
  EXPECT_EQ(num_non_sample, 2);
}