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

#include<map>

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        fprintf(stderr, "line %d: %s", __LINE__, tsk_strerror(val));                    \
        exit(EXIT_FAILURE);                                                             \
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

std::vector<int> s3_num_roots(tsk_table_collection_t &tables)
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

bool s3_mrca_found(tsk_table_collection_t &tables)
{
    std::vector<int> num_roots = s3_num_roots(tables);
    for (auto num_root : num_roots) {
        if (num_root != 1) { 
            return false; };
    }
    return true;
}

void s3_reverse_time(tsk_node_table_t &nodes, int begin, int end = 0)
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

void s3_add_time(tsk_node_table_t &nodes, int by, int begin, int end = 0)
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

void s3_dormancy_weights(std::vector<double> &weights, double b, tsk_id_t m)
{
    for (int i=0; i<m; i++) {
        weights.emplace_back(b * pow((1-b), i-1));
    }
    double sum = std::accumulate(std::begin(weights), std::end(weights), 0.0);
    for (int i=0; i<weights.size(); i++) {
        weights[i] = weights[i] / sum;
    }
}

void s3_dormancy_generation(std::vector<tsk_id_t> &dormancy_generations, std::vector<double> &weights, tsk_id_t N) 
{
    stats_random_discrete(dormancy_generations, weights, N);
}

void s3_infsites(std::vector<double> &output, double mu, std::map<double, bool> &lookup, int L, int m_i=1)
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

void s3_fsites(std::vector<double> &output, double mu, double position, std::map<double, bool> &lookup)
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

void s3_fsites_determined(std::vector<double> &output, double mu, double position, std::map<double, bool> &lookup)
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

void s3_recombination_events(std::vector<recombination_event> &recombination_events, double r,
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

tsk_id_t s3_gamete_position_recombination(double position, std::vector<recombination_event> &recombination_events)
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

void s3_simplify(tsk_table_collection_t &tables, std::vector<std::pair<tsk_id_t, MutationMetaData>> temp_mutations,
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

void s3_merge_selection_genotypes(std::map<int, genotype_matrix> &selection_genotypes_output,
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
    
    std::vector<int> int_keys(keys.begin(), keys.end());
    int_keys.erase(std::unique(int_keys.begin(), int_keys.end()), int_keys.end());

    merged_genotypes.emplace_back(current_genotype);
    for (int i=0; i<int_keys.size(); i++) {
        selection_genotypes_output[int_keys[i]] = merged_genotypes[i];
    }
}

void s3_weights(std::vector<double> &weights, std::vector<tsk_id_t> selection_genotype, 
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
              int generations_post_fixation_threshold = 1
             ) {
    
    
    int ret;
    int reached_fixation = 0;
    int generations_post_fixation = 0;
    
    std::map<double, double> position_coefficient_mapping;
    std::map<double, double> position_dominance_mapping;
    for (int i=0; i<selection_positions.size(); i++) {
        position_coefficient_mapping[selection_positions[i]] = selection_coefficients[i];
        position_dominance_mapping[selection_positions[i]] = dominance_coefficients[i];
    }
    
    //tsk_node_table_t nodes = tables.nodes;
    //tsk_edge_table_t edges = tables.edges;
    //tsk_site_table_t sites = tables.sites;
    //tsk_mutation_table_t mutations = tables.mutations;
    
    
    //mu = 4 * N * mu * L;
    //r = 4 * N * r * L;
    mu = mu * L;
    r = r * L;

    for (int i=0; i<mu_selection_rates.size(); i++) {
        mu_selection_rates[i] = 4 * N * mu_selection_rates[i];
    }
    
    
    std::map<double, genotype_matrix> selection_genotypes;
    std::map<double, genotype_matrix> next_selection_genotypes;

    std::map<double, bool> lookup;
    std::map<double, bool> lookup_selection;
    
    
    std::map<double, int> fixed_mutations;
    
    std::vector<std::pair<tsk_id_t, MutationMetaData>> temp_mutations;
    std::map<double, bool> mutation_adding_generation;
    
    std::vector<double> mpos_selection;
    std::vector<double> mpos;
    
    std::vector<tsk_id_t> parents;
    
    int gen = 0;
    int selection_generation = 0;
    
    for (int i=0; i<m; i++) {
        for (int j=0; j<2*N; j++) {
            tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, (double) gen, TSK_NULL, TSK_NULL, NULL, 0);
        }
        gen++;
    }
    
    s3_reverse_time(tables.nodes, 0);    
    s3_simplify(tables, std::vector<std::pair<tsk_id_t, MutationMetaData>>{}, m);
    
    for(int num_generation=0; num_generation<num_generations;) {

        tsk_id_t next_offspring_index = tables.nodes.num_rows;
        
        tsk_id_t num_non_sample{0};
        for (int i=0; i<tables.nodes.num_rows; i++) {
            if (tables.nodes.flags[i] != 1) {
                num_non_sample++;
            }
        }
        
        tsk_id_t first_parental_index = next_offspring_index - (2*N) - num_non_sample;
        
        for (int c=0; c<gc; c++) {
            
            std::vector<double> dorm_weights;
            s3_dormancy_weights(dorm_weights, b, m);

            std::vector<tsk_id_t> dormancy_generations;
            s3_dormancy_generation(dormancy_generations, dorm_weights, 2*N);
            
            std::vector<double> keys = get_keys(selection_genotypes);
            if (selection_activation_generation && gen > selection_activation_generation && keys.size() != 0) {
                selection_generation++;
                
                // probably should be a function
                std::map<int, genotype_matrix> selection_genotypes_output;
                s3_merge_selection_genotypes(selection_genotypes_output, selection_genotypes);
                std::vector<int> int_keys = get_keys(selection_genotypes_output);
                
                
                // update recorder
                std::vector<double> keys = get_keys(selection_genotypes);
                for (auto k : keys) {
                    genotype_matrix gm = selection_genotypes[k];
                    std::vector<tsk_id_t> genotype_current_gen = gm[m-1];
                    tsk_id_t genotype_sum = std::accumulate(std::begin(genotype_current_gen), std::end(genotype_current_gen), 0);
                    if (genotype_sum != 0) {
                        recorder.insert(k, selection_generation, genotype_sum, position_coefficient_mapping[k], position_dominance_mapping[k], keys.size(), gen);
                    }
                }
                
                
                
                for (auto k : int_keys) {
                    // check if any mutation has reached fixation
                    genotype_matrix gm = selection_genotypes_output[k];
                    std::vector<tsk_id_t> genotype_current_gen = gm[m-1];
                    
                    if (debug_print) {std::cout << selection_generation << " " << k << " : ";}
                    if (debug_print) {print_int_vector2(genotype_current_gen);}
                    
                    tsk_id_t genotype_sum = std::accumulate(std::begin(genotype_current_gen), std::end(genotype_current_gen), 0);
                    if (genotype_sum >= (2*N*percent_fixed)) {
                        
                        //s3_add_time(tables.nodes, c, 0, tables.nodes.num_rows-2*N*(c));
                        //s3_reverse_time(tables.nodes, tables.nodes.num_rows -2*N*(c), tables.nodes.num_rows);  
                        
                        //s3_simplify(tables,temp_mutations, 1);
                        //std::cout << "Fixation reached" << std::endl;
                        reached_fixation = 1;
                        fixed_mutations[k]++;
                        //return;
                    }
                }
                if (reached_fixation == 1) {
                    generations_post_fixation++;
                }
                
                
                // exit condition upon fixation
                if (reached_fixation == 1 && generations_post_fixation > generations_post_fixation_threshold) {
                    s3_add_time(tables.nodes, c, 0, tables.nodes.num_rows-2*N*(c));
                    s3_reverse_time(tables.nodes, tables.nodes.num_rows -2*N*(c), tables.nodes.num_rows);  
                    s3_simplify(tables,temp_mutations, 1);
                    return;
                }
                
                
                
                
                
                
                parents.clear();
                
                for(auto d : dormancy_generations) {
                    
                    std::vector<std::vector<double>> total_weights{};
                    for (auto k : int_keys) {
                            std::vector<double> weights;
                            s3_weights(weights, selection_genotypes_output[k][m-1-d],
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
                parents.clear();
                stats_random_discrete(parents, std::vector<double>(N, 1), 2*N);
            }
            
            tsk_id_t i_genotype = 0;
            for (int p=0; p<parents.size(); p+=2) {
                
                
                
                tsk_id_t parent1 = parents[p];
                tsk_id_t parent2 = parents[p+1];
                
                tsk_id_t dormancy_updated_first_parental_index_1 = first_parental_index - (2*N*dormancy_generations[parent1]); 
                tsk_id_t dormancy_updated_first_parental_index_2 = first_parental_index - (2*N*dormancy_generations[parent2]); 
                
                if (c > 0 && (2*N*dormancy_generations[parent1]) >= c*2*N) {
                       dormancy_updated_first_parental_index_1 -= num_non_sample;
                }
                
                if (c > 0 && (2*N*dormancy_generations[parent2]) >= c*2*N) {
                       dormancy_updated_first_parental_index_2 -= num_non_sample;
                }

                tsk_id_t p1g1 = dormancy_updated_first_parental_index_1 + 2*parent1;
                tsk_id_t p1g2 = p1g1+1;
                tsk_id_t p2g1 = dormancy_updated_first_parental_index_2 + 2*parent2;
                tsk_id_t p2g2 = p2g1+1;
                    
                double mendel1 = 0;
                double mendel2 = 0;
                stats_random_real(mendel1, 0, 1);
                stats_random_real(mendel2, 0, 1);
                
                if (mendel1 < 0.5) {
                    tsk_id_t tmp = p1g1;
                    p1g1 = p1g2;
                    p1g2 = tmp;
                }
 
                if (mendel2 < 0.5) {
                    tsk_id_t tmp = p2g1;
                    p2g1 = p2g2;
                    p2g2 = tmp;
                }
                
                
                tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, gen, TSK_NULL, TSK_NULL, NULL, 0);
                tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, gen, TSK_NULL, TSK_NULL, NULL, 0);
                    
                std::vector<recombination_event> recombination_events;
                s3_recombination_events(recombination_events, r, std::pair<tsk_id_t, tsk_id_t>{p1g1,p1g2},next_offspring_index, L);
                
                if (selection_activation_generation && gen > selection_activation_generation) {     
                    
                    
                    std::vector<double> lookup_selection_keys = get_keys(lookup_selection);
                    for (auto k : lookup_selection_keys) {
                        if (!mutation_adding_generation[k]) {
                            tsk_id_t genotype = selection_genotypes[k][m-1-dormancy_generations[parent1]][s3_gamete_position_recombination(k, recombination_events) - dormancy_updated_first_parental_index_1];
                            //std::cout << "genotype: " << genotype << std::endl;
                            next_selection_genotypes[k][0][i_genotype] = genotype;
                        }
                    }
                    for (int i=0; i<mu_selection_rates.size(); i++) {
                        mpos_selection.clear();

                        //s3_fsites(mpos_selection, mu_selection_rates[i], selection_positions[i], lookup_selection);
                        
                        
                        std::vector<double> lookup_selection_keys = get_keys(lookup_selection);
                        if (lookup_selection.size() == 0) {
                        
                            //for (auto k : lookup_selection_keys) {
                                s3_fsites_determined(mpos_selection, mu_selection_rates[i], selection_positions[i], lookup_selection);
                            //}
                        }
                        


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
                
                mpos.clear();
                if (mutations_in_seeds) {
                     s3_infsites(mpos, mu, lookup, L, dormancy_generations[parent1]+1);
                } else {
                     s3_infsites(mpos, mu, lookup, L, 1);
                }
                for (int mi=0; mi<mpos.size(); mi++) {
                    temp_mutations.emplace_back(std::pair<tsk_id_t, MutationMetaData>{next_offspring_index, {gen,mpos[mi]}});
                }
                i_genotype++;
                next_offspring_index++;
                
                
                recombination_events.clear();
                s3_recombination_events(recombination_events, r, std::pair<tsk_id_t, tsk_id_t>{p2g1,p2g2},next_offspring_index, L);
                if (selection_activation_generation && gen > selection_activation_generation) {
                    std::vector<double> lookup_selection_keys = get_keys(lookup_selection);
                    for (auto k : lookup_selection_keys) {
                        if (!mutation_adding_generation[k]) {
                            tsk_id_t genotype = selection_genotypes[k][m-1-dormancy_generations[parent2]][s3_gamete_position_recombination(k, recombination_events) - dormancy_updated_first_parental_index_2];
                            //std::cout << "genotype: " << genotype << std::endl;
                            next_selection_genotypes[k][0][i_genotype] = genotype;
                        }
                    }
                    for (int i=0; i<mu_selection_rates.size(); i++) {
                        mpos_selection.clear();

                        //s3_fsites(mpos_selection, mu_selection_rates[i], selection_positions[i], lookup_selection);

                        /*
                        std::vector<double> lookup_selection_keys = get_keys(lookup_selection);
                        for (auto k : lookup_selection_keys) {
                            if (lookup_selection.size() == 0) {
                                s3_fsites_determined(mpos_selection, mu_selection_rates[i], selection_positions[i], lookup_selection);
                            }
                        }
                        */


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
                
                mpos.clear();
                if (mutations_in_seeds) {
                     s3_infsites(mpos, mu, lookup, L, dormancy_generations[parent2]+1);
                } else {
                     s3_infsites(mpos, mu, lookup, L, 1);
                }                
                
                
                for (int mi=0; mi<mpos.size(); mi++) {
                    temp_mutations.emplace_back(std::pair<tsk_id_t, MutationMetaData>{next_offspring_index, {gen,mpos[mi]}});
                }
                i_genotype++;
                next_offspring_index++;
                
            }
            
            
            //tsk_node_table_print_state(&tables.nodes, stdout);
            
    
            
            first_parental_index = tables.nodes.num_rows - 2*N;
            
            std::vector<double> lookup_selection_keys = get_keys(lookup_selection);
            if (m > 1) {
                for (int mi=0; mi<lookup_selection_keys.size(); mi++) {
                    for (int i=1; i<m; i++) {
                        selection_genotypes[lookup_selection_keys[mi]][i-1] =  selection_genotypes[lookup_selection_keys[mi]][i]; 
                    }
                }
            }
            
            
            
            std::vector<double> lost_mutations;
            //std::vector<double> fixed_mutations;
            
            for (auto k : lookup_selection_keys) {
                //if (std::accumulate(std::begin(next_selection_genotypes[k][0]), std::end(next_selection_genotypes[k][0]),0) == 2*N) {
                    
                    //std::cout << "loop: ";
                    //print_int_vector(next_selection_genotypes[k][0]);
                    
                    //fixed_mutations.emplace_back(k);
                //}
                //else 
                if (std::accumulate(std::begin(next_selection_genotypes[k][0]), std::end(next_selection_genotypes[k][0]),0) == 0) {
                    lost_mutations.emplace_back(k);
                } else {
                    selection_genotypes[k][m-1] = next_selection_genotypes[k][0];
                    next_selection_genotypes[k] = genotype_matrix(1, std::vector(2*N, 0));
                }
            }
            
            for (int i=0; i<lost_mutations.size(); i++) {
                lookup_selection.erase(lost_mutations[i]);
                selection_genotypes.erase(lost_mutations[i]);
                next_selection_genotypes.erase(lost_mutations[i]);
            }

            
            for (auto k : lookup_selection_keys) {
                //std::cout << k << std::endl;
                //std::cout << fixed_mutations[(int) k] << std::endl;
                if (fixed_mutations[(int) k] > 50) {
                    lookup_selection.erase(k);
                    selection_genotypes.erase(k);
                    next_selection_genotypes.erase(k);
                    fixed_mutations.erase((int) k);
                }
            }

            
            /*
            for (int i=0; i<fixed_mutations.size(); i++) {
                
                std::cout << "fixed: " << fixed_mutations[i] << std::endl;
                print_int_vector(selection_genotypes[fixed_mutations[i]][m-1]);
                
                s3_add_time(tables.nodes, c+1, 0, tables.nodes.num_rows-2*N*(c+1));
                s3_reverse_time(tables.nodes, tables.nodes.num_rows -2*N*(c+1), tables.nodes.num_rows);  
                s3_simplify(tables,temp_mutations, 1);
                return;
            }
            */
            
            
            std::vector<double> mutation_adding_generation_keys = get_keys(mutation_adding_generation);
            for (int i=0; i<mutation_adding_generation_keys.size(); i++) {
                mutation_adding_generation[mutation_adding_generation_keys[i]] = false;
            }
            
            gen++;
            num_generation++; 
        }
        
        s3_add_time(tables.nodes, gc, 0, tables.nodes.num_rows-2*N*gc);
        s3_reverse_time(tables.nodes, tables.nodes.num_rows-2*N*gc, tables.nodes.num_rows);  
                
        
        
        /*
        ret = tsk_table_collection_build_index(&tables, 0);
        check_tsk_error(ret);
        tsk_treeseq_t ts;
        ret = tsk_treeseq_init(&ts, &tables, 0);
        check_tsk_error(ret);
        int num_mutations = tsk_treeseq_get_num_mutations(&ts);
        int num_samples = tsk_treeseq_get_num_samples(&ts);
        tsk_vargen_t vargen;
        tsk_variant_t *variant;
        tsk_vargen_init(&vargen, &ts, ts.samples, num_samples, NULL , 0);
        int i=0;
        while(tsk_vargen_next(&vargen, &variant) != 0) {

            double position = tables.sites.position[i];
            std::cout << i << " " << position << " : ";
            int population_abundance = 0;
            int n = 0;
            for (int s=0; s<2*N; s++) {
                int state = (int) vargen.variant.genotypes.i8[s];
                if (state == 1) {
                    n++;
                }
                std::cout << state;
            }
            std::cout << std::endl;
            if (n == 2*N) return;
            i++;
        }
        
        std::vector<double> lookup_selection_keys = get_keys(lookup_selection);
        if (lookup_selection_keys.size() != 0) {
            //print_genotypes(selection_genotypes);
    
            for (auto k : lookup_selection_keys) {
                std::cout << "k " << k << " : ";
                print_int_vector2(selection_genotypes[k][m-1]);
                std::cout << std::endl;
            }
        }
        
        */
        
        
        if (stop_after_mrca) {
            tsk_table_collection_t tables_copy;
            ret = tsk_table_collection_copy(&tables, &tables_copy, 0);
            s3_simplify(tables_copy,temp_mutations, 1);
            if (s3_mrca_found(tables_copy)) {
                s3_simplify(tables,temp_mutations, 1);
                return;
            }
        }
        
        s3_simplify(tables,temp_mutations, m);


        temp_mutations.clear();
        //num_generations--;
        //std::cout << "num_generations: " << num_generations << std::endl;
            
            
        
    }
        

    

    s3_simplify(tables, temp_mutations, 1) ;
    
}