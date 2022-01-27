#pragma once

#include <tskit.h>
#include <tskit/tables.h>
#include <tskit/trees.h>

#include <numeric>
#include <random>

#include <fstream>

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <iomanip>  // std::setprecision()

#include <new>
#include <cassert>

#include<map>

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        fprintf(stderr, "line %d: %s", __LINE__, tsk_strerror(val));                    \
        exit(EXIT_FAILURE);                                                             \
    }

void sleepy_create_mock_nodes(tsk_node_table_t *self, int length=10) {
    
    for (int i=0; i<length; i++) {
        tsk_node_table_add_row(self, TSK_NODE_IS_SAMPLE, i, TSK_NULL, TSK_NULL, NULL, 0);
    }
    
}


void sleepy_init_tables(tsk_table_collection_t *self, int N=3, int m=1, int gen=1) {
    for (int i=0; i<m; i++) {
        for (int j=0; j<2*N; j++) {
            tsk_node_table_add_row(&self->nodes, TSK_NODE_IS_SAMPLE, (double) gen,
                                   TSK_NULL, TSK_NULL, NULL, 0);
        }
        gen++;
    }
}

tsk_id_t sleepy_num_non_sample_node(tsk_node_table_t &nodes) {
    tsk_id_t num_non_sample{0};
    for (int i=0; i<nodes.num_rows; i++) {
        if (nodes.flags[i] != 1) {
            num_non_sample++;
        }
    }
    return num_non_sample;
}


struct mutation_info {
    std::vector<int> origin_generation;
    double selection_coefficient;
    double dominance_coefficient;
    std::vector<int> abs_population_freq;
    std::vector<int> num_active_mutations;
    std::vector<int> current_generation;
};

struct recorder_t {

    // keys
    std::vector<double> positions;
    std::unordered_map<double, mutation_info> data;
    void insert(double position, int origin_generation, int abs_population_freq, double selection_coefficient,double dominance_coefficient, int num_active_mutations, int current_generation) 
        
    {
        if (data.find(position) == data.end()) {
            data[position] = mutation_info{ .origin_generation = std::vector<int>{origin_generation},
                                     .selection_coefficient = selection_coefficient,
                                     .dominance_coefficient = dominance_coefficient,
                                     .abs_population_freq = std::vector<int>{abs_population_freq},
                                     .num_active_mutations = std::vector<int>{num_active_mutations},
                                     .current_generation = std::vector<int>{current_generation},
                                   };
        } else {
            data[position].abs_population_freq.emplace_back(abs_population_freq);
            data[position].num_active_mutations.emplace_back(num_active_mutations);
            data[position].origin_generation.emplace_back(origin_generation); 
            data[position].current_generation.emplace_back(current_generation); 
        }
    }



    void save(std::string file) 
    {
        std::ofstream outfile{file};
        outfile << "position,origin_generation,abs_population_freq,selection_coefficient,dominance_coefficient,active_mutations,current_generation" << std::endl;
        for (auto &it : data) { 
            std::vector<int> abs_population_freq = it.second.abs_population_freq; 
            std::vector<int> num_active_mutations = it.second.num_active_mutations; 
            std::vector<int> origin_generation = it.second.origin_generation;
            std::vector<int> current_generation = it.second.current_generation;
            double selection_coefficient = it.second.selection_coefficient;
            double dominance_coefficient = it.second.dominance_coefficient;
            

            for (int j=0; j<abs_population_freq.size(); j++) {
                double position = it.first;
                outfile << std::setprecision(10) << position << "," << origin_generation[j] << ","
                    <<  abs_population_freq[j] << "," << selection_coefficient << "," << dominance_coefficient << ","
                    << num_active_mutations[j] <<  "," << current_generation[j] << std::endl;
            }
        }
        
    }
};

std::vector<int> sleepy_num_roots(tsk_table_collection_t &tables)
{
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    int ret, iter;
    std::vector<int> num_roots{};

    ret = tsk_table_collection_build_index(&tables, 0);
    check_tsk_error(ret);
    ret = tsk_treeseq_init(&ts, &tables, 0);
    check_tsk_error(ret);
    ret = tsk_tree_init(&tree, &ts, 0); 
    check_tsk_error(ret);

    for (iter = tsk_tree_first(&tree); iter == 1; iter = tsk_tree_next(&tree)) {
        num_roots.emplace_back(tsk_tree_get_num_roots(&tree));
    }

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
    return num_roots;
}

bool sleepy_mrca_found(tsk_table_collection_t &tables)
{
    std::vector<int> num_roots = sleepy_num_roots(tables);
    for (auto num_root : num_roots) {
        if (num_root != 1) { 
            return false; };
    }
    return true;
}

void sleepy_reverse_time(tsk_node_table_t &nodes, int begin, int end = 0)
{
    std::vector<double> time_before; 
    time_before.reserve(begin);
    std::vector<double> time_mid; 
    time_before.reserve(begin+nodes.num_rows-end);
    std::vector<double> time_after;
    
    if (end == 0) end = nodes.num_rows;
    for (int i=0; i<begin; i++) time_before.emplace_back(nodes.time[i]);
    
    double max_element{0};
    for (int i=begin; i<end; i++) {
        if (nodes.time[i] >= max_element) {
            max_element = nodes.time[i];
        }
    }
    
    for (int i=begin; i<end; i++) {
        time_mid.emplace_back((nodes.time[i] - max_element) * (-1));
    }
    
    for (int i=end; i<nodes.num_rows; i++) {
        time_after.emplace_back(nodes.time[i]);
    }
    
    time_before.insert(std::end(time_before), std::begin(time_mid), std::end(time_mid));
    time_before.insert(std::end(time_before), std::begin(time_after), std::end(time_after));
    
    for (int i=0; i<nodes.num_rows; i++) {
        nodes.time[i] = time_before[i];
    }
    
}

void sleepy_add_time(tsk_node_table_t &nodes, int by, int begin, int end = 0)
{
    std::vector<double> time_before; 
    time_before.reserve(begin);
    std::vector<double> time_mid; 
    time_before.reserve(begin+nodes.num_rows-end);
    std::vector<double> time_after;
    
    if (end == 0) end = nodes.num_rows;
    for (int i=0; i<begin; i++) time_before.emplace_back(nodes.time[i]);
    

    
    for (int i=begin; i<end; i++) {
        time_mid.emplace_back(nodes.time[i] + by);
    }
    
    for (int i=end; i<nodes.num_rows; i++) {
        time_after.emplace_back(nodes.time[i]);
    }
    
    time_before.insert(std::end(time_before), std::begin(time_mid), std::end(time_mid));
    time_before.insert(std::end(time_before), std::begin(time_after), std::end(time_after));
    
    for (int i=0; i<nodes.num_rows; i++) {
        nodes.time[i] = time_before[i];
    }
    
}

std::random_device rd;
std::mt19937 generator(rd());

void stats_random_discrete(std::vector<int> &output, const std::vector<double> &weights, int n)
{
    std::discrete_distribution<> d(weights.begin(), weights.end());
    for (int i=0; i<n; i++) { output.emplace_back(d(generator)); }
}

void stats_random_poisson(std::vector<int> &output, double lambda, int n)
{
    std::poisson_distribution<> d(lambda);
    for (int i=0; i<n; i++) {
        output.emplace_back(d(generator));
    }
}

void stats_random_poisson(int &output, double lambda)
{
    std::poisson_distribution<> d(lambda);
    output = d(generator);
}

void stats_random_real(double &output, int start, int end)
{
    std::uniform_real_distribution<> d(start, end);
    output = d(generator); 
}

void sleepy_random_gamete_switch(tsk_id_t &gamete1, tsk_id_t &gamete2){
    double mendel = 0;
    stats_random_real(mendel, 0, 1);
    if (mendel < 0.5) {
        tsk_id_t tmp = gamete1;
        gamete1 = gamete2;
        gamete2 = tmp;
    }
}


void sleepy_dormancy_weights(std::vector<double> &weights, double b, tsk_id_t m)
{
    /* assert don't work in jupyter-cpp kernel */
    assert(("germination rate needs to be in [1,0[ interval.", b <= 1 && b > 0));
    
    for (int i=0; i<m; i++) {
        weights.emplace_back(b * pow((1-b), i-1));
    }
    double sum = std::accumulate(std::begin(weights), std::end(weights), 0.0);
    
    
    for (int i=0; i<weights.size(); i++) {
        weights[i] = weights[i] / sum;
        
        /* all weights will be zero if b == 1*/
        /* otherwise  weights[0] == -nan */
        
        if (b == 1.0) weights[i] = 0;
    }
}


void sleepy_dormancy_generation(std::vector<tsk_id_t> &dormancy_generations, std::vector<double> &weights, tsk_id_t N) 
{
    stats_random_discrete(dormancy_generations, weights, N);
}

void sleepy_infsites(std::vector<double> &output, double mu, std::map<double, bool> &lookup, int L, int m_i=1)
{
    int nmut = 0;
    int accum_mutations = 0;
    for (int j=0; j<m_i; j++) {
        stats_random_poisson(accum_mutations, mu);
        nmut +=accum_mutations;
    }
        
    stats_random_poisson(nmut, nmut);
    
    
    if (nmut == 0) {
        return;
    }
    int i=0;
    while (i < nmut) {
        double pos{0};
        stats_random_real(pos, 0, L);
        while (lookup.find(pos) != lookup.end()) {
            stats_random_real(pos, 0, L);
        }
        output.emplace_back(pos);
        lookup[pos] = true;
        i++;
    }
}

void sleepy_fsites(std::vector<double> &output, double mu, double position, std::map<double, bool> &lookup)
{
    int nmut;
    stats_random_poisson(nmut, mu);
    if (nmut == 0) {
        return;
    }
    int i=0;
    while (i < nmut) {
        double pos{0};
        stats_random_real(pos, 0, 1);
        pos += position;
        while (lookup.find(pos) != lookup.end()) {
            stats_random_real(pos, 0, 1);
            pos += position;
        }
        output.emplace_back(pos);
        lookup[pos] = true;
        i++;
    }
}

void sleepy_fsites_determined(std::vector<double> &output, double mu, double position, std::map<double, bool> &lookup)
{
    int nmut = 1;
    int i=0;
    while (i < nmut) {
        double pos{0};
        stats_random_real(pos, 0, 1);
        pos += position;
        while (lookup.find(pos) != lookup.end()) {
            stats_random_real(pos, 0, 1);
            pos += position;
        }
        output.emplace_back(pos);
        lookup[pos] = true;
        i++;
    }
}

template<typename T>
std::vector<T> MultiplyVec(std::vector<T> v, int n)
{
    std::vector<T> rv;
    rv.reserve(n * v.size());
    for (int i=0; i<n; i++){
        for (auto e : v){ rv.emplace_back(e); }
    }
    return rv;
} 

template<typename T>
std::vector<T> get_keys(const std::map<T, std::vector<tsk_id_t>>& map)
{
    std::vector<T> keys;
    keys.reserve(map.size());
    for (const auto& it : map)
        keys.push_back(it.first);
    return keys;
}

struct recombination_event {
    double left;
    double right;
    tsk_id_t parent;
    tsk_id_t child;
};

void sleepy_recombination_events(std::vector<recombination_event> &recombination_events, double r,
                       std::pair<tsk_id_t, tsk_id_t> parent_idxs, tsk_id_t next_offspring_id, double L)
{
    int nbreaks;
    stats_random_poisson(nbreaks, r);

    if (nbreaks == 0) {
                
        recombination_event re = {0, L, parent_idxs.first, next_offspring_id};
        recombination_events.emplace_back(re);
        
        
    } else {
        std::vector<double> b;
        b.reserve(nbreaks);
        int i{0};
        while (i < nbreaks){
            double p;
            stats_random_real(p, 0, L);
            b.emplace_back(p);
            i++;
        }
        std::sort(b.begin(), b.end());
        if (b.back() != L) { b.emplace_back(L); }
        if (*b.begin() != 0.0) { 
            b.insert(b.begin(), 0.0);
        } else {
            int tmp_parent_idx = parent_idxs.first;
            parent_idxs.first = parent_idxs.second;
            parent_idxs.second = tmp_parent_idx;
        }
        std::vector<int> p_idxs{parent_idxs.first, parent_idxs.second};
        std::vector<int> pgams = MultiplyVec(p_idxs, b.size() / 2);
        for (int i=0; i<b.size()-1; i++) {
            recombination_event re = {b[i], b[i+1], pgams[i], next_offspring_id};
            recombination_events.emplace_back(re);
            //std::cout <<  b[i] << " " << b[i+1] << " " <<  pgams[i] << std::endl;
        }  
    }
}

/* returns the parent of the mutation */
tsk_id_t sleepy_gamete_position_recombination(double position, std::vector<recombination_event> &recombination_events)
{
    for (int i=0; i<recombination_events.size(); i++) {
        if (position >= recombination_events[i].left && position <= recombination_events[i].right) {
            return recombination_events[i].parent;
        }
    }
    return -9;
}

struct MutationMetaData {
    tsk_id_t origin;
    double position;
};

void sleepy_simplify(tsk_table_collection_t &tables, std::vector<std::pair<tsk_id_t, MutationMetaData>> temp_mutations,
              tsk_id_t sample_generation)
{
    int ret;
    std::vector<tsk_id_t> samples;
    for (int i=0; i<tables.nodes.num_rows; i++) {
        if ((int) tables.nodes.time[i] < sample_generation) {
            samples.emplace_back(i);
        }
    }
    for (int i=0; i<temp_mutations.size(); i++) {
        tsk_site_table_add_row(&tables.sites, temp_mutations[i].second.position, "0", 1, NULL, 0);
        
        tsk_mutation_table_add_row(&tables.mutations, tables.sites.num_rows-1,
                                   temp_mutations[i].first, TSK_NULL, TSK_UNKNOWN_TIME , "1", 1, NULL, 0);
    }
    
    
    tsk_id_t* samples_array = &samples[0];
    ret = tsk_table_collection_sort(&tables, NULL, 0); 
    //std::cout << "ret sort: " << std::endl;
    check_tsk_error(ret);
    ret = tsk_table_collection_simplify(&tables, samples_array, samples.size(), TSK_FILTER_SITES, NULL); 
    //std::cout << "ret simp: " << std::endl;
    check_tsk_error(ret);

}

/*
template<typename T>
std::vector<T> get_keys(const std::map<T, std::vector<tsk_id_t>>& map)
{
    std::vector<T> keys;
    keys.reserve(map.size());
    for (const auto& it : map)
        keys.push_back(it.first);
    return keys;
}
*/

typedef std::vector<std::vector<tsk_id_t>> genotype_matrix;

template<typename T>
std::vector<T> get_keys(const std::map<T, genotype_matrix>& map)
{
    std::vector<T> keys;
    keys.reserve(map.size());
    for (const auto& it : map)
        keys.push_back(it.first);
    return keys;
}

template<typename T>
std::vector<T> get_keys(const std::map<T, bool>& map)
{
    std::vector<T> keys;
    keys.reserve(map.size());
    for (const auto& it : map)
        keys.push_back(it.first);
    return keys;
}

void print_double_vector(std::vector<double> v) {
    for (int i=0; i<v.size(); i++) {
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}

void print_int_vector(std::vector<int> v) {
    for (int i=0; i<v.size(); i++) {
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}

void print_int_vector2(std::vector<int> v) {
    for (int i=0; i<v.size(); i++) {
        std::cout << v[i] << "";
    }
    std::cout << std::endl;
}

void print_double_vector_vector(std::vector<std::vector<double>> v) {
    for (int i=0; i<v.size(); i++) {
        print_double_vector(v[i]);
    }
    std::cout << std::endl;
}

void sleepy_merge_selection_genotypes(std::map<int, genotype_matrix> &selection_genotypes_output,
    std::map<double, genotype_matrix> selection_genotypes)
{
    std::vector<double> keys = get_keys(selection_genotypes);
    double current_key = keys[0];
    genotype_matrix current_genotype = selection_genotypes[current_key];
    std::vector<genotype_matrix> merged_genotypes;
    
    for (int i=1; i<keys.size(); i++) {
        if (((int) current_key) == ((int) keys[i])) {
            
            for (int k=0; k<selection_genotypes[keys[i]].size(); k++){
                for (int j=0; j<selection_genotypes[keys[i]][k].size(); j++) {
                    if (selection_genotypes[keys[i]][k][j] > 0 || current_genotype[k][j] > 0) {
                        current_genotype[k][j] = 1;
                    }
                }
            }
            
        } else {
            merged_genotypes.emplace_back(current_genotype);
            current_key = keys[i];
            current_genotype = selection_genotypes[keys[i]];
        }
    }
    
    std::vector<int> current_active_generations_int(keys.begin(), keys.end());
    current_active_generations_int.erase(std::unique(current_active_generations_int.begin(), current_active_generations_int.end()), current_active_generations_int.end());

    merged_genotypes.emplace_back(current_genotype);
    for (int i=0; i<current_active_generations_int.size(); i++) {
        selection_genotypes_output[current_active_generations_int[i]] = merged_genotypes[i];
    }
}

void sleepy_weights(std::vector<double> &weights, std::vector<tsk_id_t> selection_genotype, 
                double dominance_coefficient, double selection_coefficient)
{
    for (int i=0; i<selection_genotype.size(); i+=2) {
        tsk_id_t genotype_A = selection_genotype[i];
        tsk_id_t genotype_a = selection_genotype[i+1];
        if (genotype_A == 0 && genotype_a == 0) {
            weights.emplace_back(1);
        } else if ((genotype_A == 1 && genotype_a == 0) || (genotype_A == 0 && genotype_a == 1)) {
            weights.emplace_back(1.0 + (2.0 * dominance_coefficient * selection_coefficient));
        } else {
            weights.emplace_back(1.0 + (2.0 * selection_coefficient));
        }
    }

    // check for fixation
    /*
    if (selection_genotype.size() == std::accumulate(selection_genotype.begin(), selection_genotype.end(), 0)) {
        for (int i=0; i<weights.size(); i++) {
            weights[i] = 1;
        }
    }
    */

}

template<typename T>
void print_genotypes(std::map<T, genotype_matrix> g)
{
    std::vector<T> keys = get_keys(g);
    for (auto k : keys) {
        std::cout << k << std::endl;
        for (int i=0; g[k].size(); i++) {
            print_int_vector(g[k][i]);
        }
    }
}

void dormancy(tsk_table_collection_t &tables,
              recorder_t &recorder,
              tsk_id_t num_generations,
              tsk_id_t N, 
              tsk_id_t m, 
              double b,
              tsk_id_t gc,
              double mu,
              double r,
              double L,
              std::vector<double> mu_selection_rates,
              std::vector<double> selection_coefficients,
              std::vector<double> dominance_coefficients,
              std::vector<double> selection_positions,
              tsk_id_t selection_activation_generation = 0,
              bool stop_after_mrca = false,
              bool mutations_in_seeds = true,
              bool debug_print = false,
              double percent_fixed = 1.0,
              int generations_post_fixation_threshold = 1,
              bool add_mutations_after_fixation = true
             ) {
    
    int ret = 0;
    int reached_fixation = 0;
    int generations_post_fixation = 0;

    bool keep_adding_mutation = true;
    
    std::map<double, double> position_coefficient_mapping;
    std::map<double, double> position_dominance_mapping;
    for (int i=0; i<selection_positions.size(); i++) {
        position_coefficient_mapping[selection_positions[i]] = selection_coefficients[i];
        position_dominance_mapping[selection_positions[i]] = dominance_coefficients[i];
    }

    //mu = 4 * N * mu * L;
    //r = 4 * N * r * L;
    mu = mu * L;
    r = r * L;

    for (int i=0; i<mu_selection_rates.size(); i++) {
        mu_selection_rates[i] = 4 * N * mu_selection_rates[i];
    }
    
    /* selection genotypes stores the genotypes of mutations under selection for m generations */
    std::map<double, genotype_matrix> selection_genotypes;
    std::map<double, genotype_matrix> next_selection_genotypes;

    std::map<double, bool> lookup;
    std::map<double, bool> lookup_selection;
    std::map<double, int> fixed_mutations;
    std::map<double, int> prev_generation;

    std::vector<std::pair<tsk_id_t, MutationMetaData>> temp_mutations;
    std::map<double, bool> mutation_adding_generation;
    
    /* vectors storing position of mutations */
    std::vector<double> mpos_selection;
    std::vector<double> mpos;

    std::vector<tsk_id_t> parents;
    
    int gen = 0;
    int selection_generation = 0;
    
    // adding first 2*N*m nodes to the node table
    sleepy_init_tables(&tables, N, m, gen);

    // updating generation counter
    gen = m;

    sleepy_reverse_time(tables.nodes, 0);    
    sleepy_simplify(tables, std::vector<std::pair<tsk_id_t, MutationMetaData>>{}, m);
    
    for(int num_generation=0; num_generation<num_generations;) {

        tsk_id_t next_offspring_index = tables.nodes.num_rows;
        tsk_id_t num_non_sample = sleepy_num_non_sample_node(tables.nodes);
        tsk_id_t first_parental_index = next_offspring_index - (2*N) - num_non_sample; /* num_non_sample will be 0 during first iteration */
        
        for (int c=0; c<gc; c++) {

            if (reached_fixation == 1) {
                generations_post_fixation++;
                if (add_mutations_after_fixation) { keep_adding_mutation = true; } 
                else { keep_adding_mutation = false; }
            }
            
            /* exit condition upon fixation */
            if (reached_fixation == 1 && generations_post_fixation > generations_post_fixation_threshold) {
                sleepy_add_time(tables.nodes, c, 0, tables.nodes.num_rows-2*N*(c));
                sleepy_reverse_time(tables.nodes, tables.nodes.num_rows -2*N*(c), tables.nodes.num_rows);  
                sleepy_simplify(tables,temp_mutations, 1);
                return;
            }

            std::vector<double> dorm_weights; sleepy_dormancy_weights(dorm_weights, b, m);
            std::vector<tsk_id_t> dormancy_generations; sleepy_dormancy_generation(dormancy_generations, dorm_weights, 2*N);
            std::vector<double> keys = get_keys(selection_genotypes);

            /* selection and at least one active selective mutation */
            if (selection_activation_generation && gen > selection_activation_generation && keys.size() != 0) {
                selection_generation++;
                
                std::map<int, genotype_matrix> selection_genotypes_output;
                sleepy_merge_selection_genotypes(selection_genotypes_output, selection_genotypes);

                /* usually should be just one */
                std::vector<int> current_active_generations_int = get_keys(selection_genotypes_output);
                assert(("current active generations is not equal to one", current_active_generations_int.size() == 1));
                
                /* update recorder */
                std::vector<double> keys = get_keys(selection_genotypes);
                for (auto k : keys) {
                    genotype_matrix gm = selection_genotypes[k];
                    std::vector<tsk_id_t> genotype_current_gen = gm[m-1];
                    tsk_id_t genotype_sum = std::accumulate(std::begin(genotype_current_gen), std::end(genotype_current_gen), 0);
                    if (genotype_sum != 0) {
                        recorder.insert(k, selection_generation, genotype_sum, position_coefficient_mapping[k], position_dominance_mapping[k], keys.size(), gen);
                    }
                }

                /* check if any mutation has reached fixation */
                for (auto k : current_active_generations_int) {
                    
                    genotype_matrix gm = selection_genotypes_output[k];
                    std::vector<tsk_id_t> genotype_current_gen = gm[m-1];
                    
                    if (debug_print) {std::cout << selection_generation << " " << k << " : ";}
                    if (debug_print) {print_int_vector2(genotype_current_gen);}
                    
                    tsk_id_t genotype_sum = std::accumulate(std::begin(genotype_current_gen), std::end(genotype_current_gen), 0);
                    if (genotype_sum >= (2*N*percent_fixed)) {
                        
                        // currently only works for single sweep
                        if (prev_generation.count(k) && num_generation - prev_generation[k] == 1) { fixed_mutations[k]++; } 
                        else { fixed_mutations[k] == 1; }

                        if (fixed_mutations[k] >= 50) { 
                            //std::cout << "fix" << std::endl;
                            reached_fixation = 1; 
                        }
                        prev_generation[k] = num_generation;
                    }
                }

                /* nessessary for recovery */
                if (reached_fixation == 1) {
                    generations_post_fixation++;
                    if (add_mutations_after_fixation) { keep_adding_mutation = true; } 
                    else { keep_adding_mutation = false; }
                }
                
                /* exit condition upon fixation */
                if (reached_fixation == 1 && generations_post_fixation > generations_post_fixation_threshold) {
                    sleepy_add_time(tables.nodes, c, 0, tables.nodes.num_rows-2*N*(c));
                    sleepy_reverse_time(tables.nodes, tables.nodes.num_rows -2*N*(c), tables.nodes.num_rows);  
                    sleepy_simplify(tables,temp_mutations, 1);
                    return;
                }
                
                
                /* retrieving parents based on genotypes associated selection and dominance coefficents from the given generation*/
                parents.clear();
                for(auto d : dormancy_generations) {
                    std::vector<std::vector<double>> total_weights{};
                    for (auto k : current_active_generations_int) {
                            std::vector<double> weights;
                            sleepy_weights(weights, selection_genotypes_output[k][m-1-d],
                                      position_dominance_mapping[k],
                                      position_coefficient_mapping[k]);                      
                            total_weights.emplace_back(weights);   
                    }
                    std::vector<double> multiplied_weights;
                    for (int j=0; j<total_weights[0].size(); j++) {
                        double weight = 1;
                        for (int i=0; i<total_weights.size(); i++) {
                            weight *= total_weights[i][j];
                        }
                        multiplied_weights.emplace_back(weight);
                    }
                    
                    double sum = std::accumulate(std::begin(multiplied_weights), std::end(multiplied_weights), 0.0);
                    for (int i=0; i<multiplied_weights.size(); i++) {
                        multiplied_weights[i] = (double) multiplied_weights[i] / sum;
                        
                    }
                    std::vector<tsk_id_t> parent;
                    stats_random_discrete(parent, multiplied_weights, 1);
                    parents.emplace_back(parent[0]);
                }

            } else {

                /* no selection */
                parents.clear();
                stats_random_discrete(parents, std::vector<double>(N, 1), 2*N);
            }
            
            tsk_id_t i_genotype = 0;
            for (int p=0; p<parents.size(); p+=2) {
                
                /* way too complicated dormancy index updating dependent on state in the simplication process */
                tsk_id_t parent1 = parents[p];
                tsk_id_t parent2 = parents[p+1];
                tsk_id_t dormancy_updated_first_parental_index_1 = first_parental_index - (2*N*dormancy_generations[parent1]); 
                tsk_id_t dormancy_updated_first_parental_index_2 = first_parental_index - (2*N*dormancy_generations[parent2]); 
                if (c > 0 && (2*N*dormancy_generations[parent1]) >= c*2*N) { dormancy_updated_first_parental_index_1 -= num_non_sample; }
                if (c > 0 && (2*N*dormancy_generations[parent2]) >= c*2*N) { dormancy_updated_first_parental_index_2 -= num_non_sample; }
                tsk_id_t p1g1 = dormancy_updated_first_parental_index_1 + 2*parent1;
                tsk_id_t p1g2 = p1g1+1;
                tsk_id_t p2g1 = dormancy_updated_first_parental_index_2 + 2*parent2;
                tsk_id_t p2g2 = p2g1+1;

                /* randomly change the gamete order */
                sleepy_random_gamete_switch(p1g1, p1g2);
                sleepy_random_gamete_switch(p2g1, p2g2);

                /* adding next generation to node table */
                tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, gen, TSK_NULL, TSK_NULL, NULL, 0);
                tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, gen, TSK_NULL, TSK_NULL, NULL, 0);
                    
                std::vector<recombination_event> recombination_events;
                sleepy_recombination_events(recombination_events, r, std::pair<tsk_id_t, tsk_id_t>{p1g1,p1g2},next_offspring_index, L);
                
                if (selection_activation_generation && gen > selection_activation_generation) {     
                    
                    /* get_keys simply retries all mutations currently under selection and saves them in lookup_selection_keys */
                    std::vector<double> lookup_selection_keys = get_keys(lookup_selection);

                    for (auto k : lookup_selection_keys) {
                        
                        /* mutation_adding_generation checks if it is the first generation were mutation is added; so it doesn't get added to multiple individuals */
                        /* tldr: only runs if mutation was already introduced below */
                        if (!mutation_adding_generation[k]) {

                            int mutation_generation = m-1-dormancy_generations[parent1];
                            int mutation_parent = sleepy_gamete_position_recombination(k, recombination_events);

                            /* parent idx has to be between 0 and 2*N */
                            mutation_parent = mutation_parent - dormancy_updated_first_parental_index_1;

                            tsk_id_t genotype = selection_genotypes[k][mutation_generation][mutation_parent];

                            /* i_genotype counter variable*/
                            next_selection_genotypes[k][0][i_genotype] = genotype; 
                        }
                    }
                    
                    /* introducing new mutations under selection */
                    /* in original code multiple mutations could be introduced each mutation having a different mutation rate */
                    /* code was changed to just add a single mutation at a fixed position independet of rate */
                    for (int i=0; i<mu_selection_rates.size(); i++) {
                        mpos_selection.clear();

                        sleepy_fsites(mpos_selection, mu_selection_rates[i], selection_positions[i], lookup_selection);
                        
                        /* retrieve mutations currently under selection */
                        std::vector<double> lookup_selection_keys = get_keys(lookup_selection);


                        /* if no other mutation is present in currently add a mutation*/
                        if (lookup_selection.size() == 0) {
                            if (keep_adding_mutation) {
                                sleepy_fsites_determined(mpos_selection, mu_selection_rates[i], selection_positions[i], lookup_selection);
                            }
                        }
                        
                        for (auto mi : mpos_selection) {
                            recorder.insert(mi, selection_generation, 1, selection_coefficients[i], dominance_coefficients[i], lookup_selection.size(), gen);
                            mutation_adding_generation[mi] = true;
                            temp_mutations.emplace_back(std::pair<tsk_id_t, MutationMetaData>{next_offspring_index, {gen,mi}});

                            /* keeping track of genotypes independently of tskit */
                            /* empty genotype matrix */
                            selection_genotypes[mi] = genotype_matrix(m, std::vector(2*N, 0));

                            /* saving next selection genotypes with added genotype */
                            next_selection_genotypes[mi] = genotype_matrix(1, std::vector(2*N, 0));

                            next_selection_genotypes[mi][0][i_genotype] = 1;
                        }
                    }                    
                }
                
                /* adding edges */
                for (int re=0; re<recombination_events.size(); re++) {
                    ret = tsk_edge_table_add_row(&tables.edges, recombination_events[re].left, recombination_events[re].right, recombination_events[re].parent, recombination_events[re].child, NULL, 0);
                }
                
                /* mutations in seeds currently never used */
                mpos.clear();
                if (mutations_in_seeds) { sleepy_infsites(mpos, mu, lookup, L, dormancy_generations[parent1]+1);} 
                else { sleepy_infsites(mpos, mu, lookup, L, 1); }

                for (int mi=0; mi<mpos.size(); mi++) {
                    temp_mutations.emplace_back(std::pair<tsk_id_t, MutationMetaData>{next_offspring_index, {gen,mpos[mi]}});
                }
                i_genotype++;
                next_offspring_index++;

                /* same as above but currenly no mutation is added */
                
                recombination_events.clear();
                sleepy_recombination_events(recombination_events, r, std::pair<tsk_id_t, tsk_id_t>{p2g1,p2g2},next_offspring_index, L);

                /* no selective mutation is added on the second parent */
                if (selection_activation_generation && gen > selection_activation_generation) {
                    std::vector<double> lookup_selection_keys = get_keys(lookup_selection);
                    for (auto k : lookup_selection_keys) {
                        if (!mutation_adding_generation[k]) {
                            tsk_id_t genotype = selection_genotypes[k][m-1-dormancy_generations[parent2]][sleepy_gamete_position_recombination(k, recombination_events) - dormancy_updated_first_parental_index_2];
                            //std::cout << "genotype: " << genotype << std::endl;
                            next_selection_genotypes[k][0][i_genotype] = genotype;
                        }
                    }
                    for (int i=0; i<mu_selection_rates.size(); i++) {
                        mpos_selection.clear();

                        for (auto mi : mpos_selection) {
                            recorder.insert(mi, selection_generation, 1, selection_coefficients[i], dominance_coefficients[i], lookup_selection.size(), gen);
                            mutation_adding_generation[mi] = true;
                            temp_mutations.emplace_back(std::pair<tsk_id_t, MutationMetaData>{next_offspring_index, {gen,mi}});
                            selection_genotypes[mi] = genotype_matrix(m, std::vector(2*N, 0));
                            next_selection_genotypes[mi] = genotype_matrix(1, std::vector(2*N, 0));
                            next_selection_genotypes[mi][0][i_genotype] = 1;
                        }
                    }
                }
  
                for (int re=0; re<recombination_events.size(); re++) {
                    ret = tsk_edge_table_add_row(&tables.edges, recombination_events[re].left, recombination_events[re].right, recombination_events[re].parent, recombination_events[re].child, NULL, 0);
                }
                
                /* neutral mutations */
                mpos.clear();
                if (mutations_in_seeds) { sleepy_infsites(mpos, mu, lookup, L, dormancy_generations[parent2]+1); } 
                else { sleepy_infsites(mpos, mu, lookup, L, 1); }                
                for (int mi=0; mi<mpos.size(); mi++) { temp_mutations.emplace_back(std::pair<tsk_id_t, MutationMetaData>{next_offspring_index, {gen,mpos[mi]}}); }

                i_genotype++;
                next_offspring_index++;
            }
                            
            first_parental_index = tables.nodes.num_rows - 2*N;
            /* get keys currenlty under selection */

            std::vector<double> lookup_selection_keys = get_keys(lookup_selection);

            /* shifting  selection_genotypes one generation into the past */
            if (m > 1) {
                for (int mi=0; mi<lookup_selection_keys.size(); mi++) {
                    for (int i=1; i<m; i++) {

                        /* last selection_genotypes will stay the same */
                        /* first row will get overwritten */
                        selection_genotypes[lookup_selection_keys[mi]][i-1] =  selection_genotypes[lookup_selection_keys[mi]][i]; 
                    }
                }
            }
            
            std::vector<double> lost_mutations;
            
            /*
            for (auto k : lookup_selection_keys) {
                if (std::accumulate(std::begin(next_selection_genotypes[k][0]), std::end(next_selection_genotypes[k][0]),0) == 0) {
                    lost_mutations.emplace_back(k);

                } else {
                    selection_genotypes[k][m-1] = next_selection_genotypes[k][0];
                    next_selection_genotypes[k] = genotype_matrix(1, std::vector(2*N, 0));
                }
            }
            */

            for (auto k : lookup_selection_keys) {

                /* updating last row of genotypes */
                selection_genotypes[k][m-1] = next_selection_genotypes[k][0];

                /* zeroing out next_selection genotypes */
                next_selection_genotypes[k] = genotype_matrix(1, std::vector(2*N, 0));

                /* iterating over all dormancy generation */
                int haplotype_sum = 0;
                if (b < 1) { for (int j=0; j<m; j++) { haplotype_sum += std::accumulate(std::begin(selection_genotypes[k][j]), std::end(selection_genotypes[k][j]),0); }
                } else {
                    /* only check last generation */
                    haplotype_sum = std::accumulate(std::begin(selection_genotypes[k][m-1]), std::end(selection_genotypes[k][m-1]),0); 
                }
                    
                if (haplotype_sum == 0) { lost_mutations.emplace_back(k); }
            }
            
            for (int i=0; i<lost_mutations.size(); i++) {
                lookup_selection.erase(lost_mutations[i]);
                selection_genotypes.erase(lost_mutations[i]);
                next_selection_genotypes.erase(lost_mutations[i]);                
            }

            for (auto k : lookup_selection_keys) {
                if (fixed_mutations[(int) k] > 50) {
                    lookup_selection.erase(k);
                    selection_genotypes.erase(k);
                    next_selection_genotypes.erase(k);
                    fixed_mutations.erase((int) k);
                }
            }
            
            std::vector<double> mutation_adding_generation_keys = get_keys(mutation_adding_generation);
            for (int i=0; i<mutation_adding_generation_keys.size(); i++) {
                mutation_adding_generation[mutation_adding_generation_keys[i]] = false;
            }
            
            gen++;
            num_generation++; 
        }
        
        sleepy_add_time(tables.nodes, gc, 0, tables.nodes.num_rows-2*N*gc);
        sleepy_reverse_time(tables.nodes, tables.nodes.num_rows-2*N*gc, tables.nodes.num_rows);  
                
        if (stop_after_mrca) {
            tsk_table_collection_t tables_copy;
            ret = tsk_table_collection_copy(&tables, &tables_copy, 0);
            sleepy_simplify(tables_copy,temp_mutations, 1);
            if (sleepy_mrca_found(tables_copy)) {
                sleepy_simplify(tables,temp_mutations, 1);
                return;
            }
        }
        
        sleepy_simplify(tables,temp_mutations, m);

        temp_mutations.clear();
        //num_generations--;
        //std::cout << "num_generations: " << num_generations << std::endl;            
    }
        
    sleepy_simplify(tables, temp_mutations, 1) ;
    
}