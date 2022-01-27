

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

TEST(test_sleepy_dormancy_weights, sleepy_dormancy_weights) {

  int m = 5;
  double b = 1.0;
  std::vector<double> dorm_weights;
  sleepy_dormancy_weights(dorm_weights, b, m);

  ASSERT_EQ(dorm_weights.size(), m) << "dorm_weights: " << dorm_weights.size() << "m " << m;


  //dorm_weights.clear();
  //b = 1.2;
  //ASSERT_DEATH(sleepy_dormancy_weights(dorm_weights, b, m), "");

}

TEST(test_sleepy_dormancy_generation, sleepy_dormancy_generation) {

  int m = 5;
  double b = 1.0;
  int N = 2;
  std::vector<double> dorm_weights;
  sleepy_dormancy_weights(dorm_weights, b, m);
  std::vector<tsk_id_t> dormancy_generations; 
  sleepy_dormancy_generation(dormancy_generations, dorm_weights, 2*N);
  ASSERT_EQ(dormancy_generations.size(), 2*N);
  

}

TEST(test_sleepy_recombination_events, sleepy_recombination_events) {
  std::vector<recombination_event> recombination_events;
  double r = 1;
  std::pair<tsk_id_t, tsk_id_t> parent_idxs = {3, 4};
  tsk_id_t next_offspring_id = 12;
  double L = 10;

  sleepy_recombination_events(recombination_events, r, parent_idxs, next_offspring_id, L);
  EXPECT_EQ(recombination_events[0].left, 0);
  EXPECT_EQ(recombination_events[recombination_events.size()-1].right, L);

}
