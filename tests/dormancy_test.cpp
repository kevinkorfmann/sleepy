

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
