#include <iostream>

#include <string>
#include <cstdlib>


#include "dormancy.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

void create_dir(std::string &outdir)
{
    std::string cmd_create_outdir = "mkdir -p " + outdir;
    //std::cout << cmd_create_outdir << std::endl;
    int ret = system(cmd_create_outdir.c_str());

}



int main(int argc, char *argv[])
{
    po::options_description description("\nUsage");
    description.add_options()
        ("help", "produce help message")
        ("num_generations", po::value<int>()->default_value(1000), "num_generations = num_generations * gc")
        ("N", po::value<int>()->default_value(50), "population size -> 2N haplotypes")
        ("m", po::value<int>()->default_value(1), "m")
        ("b", po::value<double>()->default_value(1), "b")
        ("gc", po::value<int>()->default_value(20), "gc")
        ("generations_post_fixation_threshold", po::value<int>()->default_value(1), "generations_post_fixation_threshold")
        ("add_mutations_after_fixation", po::value<bool>()->default_value(true), "add_mutations_after_fixation")
        ("mu", po::value<double>()->default_value(0), "mu")
        ("r", po::value<double>()->default_value(0.5e-6), "r")
        ("L", po::value<double>()->default_value(5000), "L")
        ("mu_selection_rates", po::value<std::vector<double> >()->multitoken()->default_value(std::vector<double>{5e-8}, ""), "mu_selection_rates")
        ("selection_coefficients", po::value<std::vector<double> >()->multitoken()->default_value(std::vector<double>{2}, ""), "selection_coefficients")
        ("dominance_coefficients", po::value<std::vector<double> >()->multitoken()->default_value(std::vector<double>{0.5}, ""), "dominance_coefficients")
        ("selection_positions", po::value<std::vector<double> >()->multitoken()->default_value(std::vector<double>{2500}, ""), "selection_positions")
        ("selection_activation_generation", po::value<int>()->default_value(500), "selection_activation_generation")
        ("mutation_in_seeds", po::value<bool>()->default_value(false), "mutation_in_seeds")
        ("stop_after_mrca", po::value<bool>()->default_value(false), "stop_after_mrca")
        ("percent_fixed", po::value<double>()->default_value(0), "percent_fixed")
        ("debug_print", po::value<bool>()->default_value(false), "debug_print")
        ("output_name", po::value<std::string>()->default_value("run"), "output_name")
        ("output_directory", po::value<std::string>()->default_value("./"), "output_directory")

    ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(description).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);
    po::notify(vm);

    if (vm.count("help") || argc<2) {
        std::cout << "sleepy (selection) simulator - (s3)" << std::endl;
        std::cout << description << "\n";
        return 1;
    }






    tsk_id_t num_generations = vm["num_generations"].as<int>();
    tsk_id_t N = vm["N"].as<int>();


    tsk_id_t m = vm["m"].as<int>();
    int generations_post_fixation_threshold = vm["generations_post_fixation_threshold"].as<int>();
    bool add_mutations_after_fixation = vm["add_mutations_after_fixation"].as<bool>();


    double b =  vm["b"].as<double>();
    tsk_id_t gc = vm["gc"].as<int>();
    double mu = vm["mu"].as<double>();
    double r = vm["r"].as<double>();
    double L = vm["L"].as<double>();
    double percent_fixed = vm["percent_fixed"].as<double>();

    std::vector<double> mu_selection_rates = vm["mu_selection_rates"].as<std::vector<double>>();
    std::vector<double> selection_coefficients = vm["selection_coefficients"].as<std::vector<double>>();
    std::vector<double> dominance_coefficients = vm["dominance_coefficients"].as<std::vector<double>>();
    std::vector<double> selection_positions = vm["selection_positions"].as<std::vector<double>>();
    tsk_id_t selection_activation_generation = vm["selection_activation_generation"].as<int>();
    bool mutation_in_seeds = vm["mutation_in_seeds"].as<bool>();
    bool stop_after_mrca = vm["stop_after_mrca"].as<bool>();
    bool debug_print = vm["debug_print"].as<bool>();
    std::string output_name = vm["output_name"].as<std::string>();
    std::string output_directory = vm["output_directory"].as<std::string>();

    create_dir(output_directory);

    recorder_t recorder;

    int ret;
    tsk_table_collection_t tables;

    ret = tsk_table_collection_init(&tables, 0);
    tables.sequence_length = L;

    check_tsk_error(ret);
    //tsk_table_collection_free(&tables);


    dormancy(tables,recorder,num_generations, N, m, b, gc, mu, r, L,
            mu_selection_rates,
            selection_coefficients,
            dominance_coefficients,
            selection_positions,
            selection_activation_generation,
            stop_after_mrca,
            mutation_in_seeds,
            debug_print,
            percent_fixed,
            generations_post_fixation_threshold,
            add_mutations_after_fixation
            );
        

    tsk_table_collection_build_index(&tables, 0);
    tsk_table_collection_dump(&tables, (output_directory + output_name + ".trees").c_str(), 0);
    recorder.save(output_directory + output_name + ".csv");


    return 0;
}