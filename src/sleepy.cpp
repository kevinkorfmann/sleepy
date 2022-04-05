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
        ("help", "produce help message.")
        ("num_generations", po::value<int>()->default_value(100000000), "maximum generation run-time (num_generations = num_generations * gc)")
        ("N", po::value<int>()->default_value(500), "population size -> 2N haplotypes")
        ("m", po::value<int>()->default_value(100), "upper-bound for seed resusication")
        ("b", po::value<double>()->default_value(1), "germination rate")
        ("gc", po::value<int>()->default_value(50), "garbage collection")
        ("generations_post_fixation_threshold", po::value<int>()->default_value(1), "run for n number of generations after fixation event has occured")
        ("add_mutations_after_fixation", po::value<bool>()->default_value(true), "further adding selective mutations after fixation event; has to be false when studying recovery of sweeps")
        //("mu", po::value<double>()->default_value(0), "mutation rate")
        ("r", po::value<double>()->default_value(5e-5), "recombination rate")
        ("L", po::value<double>()->default_value(10000), "mapping length")
        //("mu_selection_rates", po::value<std::vector<double> >()->multitoken()->default_value(std::vector<double>{5e-8}, ""), "mu_selection_rates")
        ("selection_coefficient", po::value<std::vector<double> >()->multitoken()->default_value(std::vector<double>{1}, ""), "selection coefficient")
        ("dominance_coefficient", po::value<std::vector<double> >()->multitoken()->default_value(std::vector<double>{0.5}, ""), "dominance coefficient")
        ("selection_position", po::value<std::vector<double> >()->multitoken()->default_value(std::vector<double>{5000}, ""), "selection position")
        ("selection_activation_generation", po::value<int>()->default_value(500), "generation at which point to start introductin selective mutations")
        //("mutation_in_seeds", po::value<bool>()->default_value(false), "mutation in seeds possible or not")
        ("stop_after_mrca", po::value<bool>()->default_value(false), "stop simulation after mrca is found")
        //("percent_fixed", po::value<double>()->default_value(0), "percent_fixed")
        ("debug_print", po::value<bool>()->default_value(false), "used for debugging during developement")
        ("output_name", po::value<std::string>()->default_value("run"), "output file nume part")
        ("output_directory", po::value<std::string>()->default_value("./"), "output directory")

    ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(description).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);
    po::notify(vm);

    if (vm.count("help") || argc<2) {
        std::cout << "sleepy (selection) simulator" << std::endl;
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
    double mu = 0; //vm["mu"].as<double>();
    double r = vm["r"].as<double>();
    double L = vm["L"].as<double>();

    /* percent fixed set to 1.0 only full sweeps are subject of this study */
    double percent_fixed = 1.0; //vm["percent_fixed"].as<double>();

    /* rates are not used anymore, therefore initialized as empty */
    std::vector<double> mu_selection_rates{0};// vm["mu_selection_rates"].as<std::vector<double>>();

    std::vector<double> selection_coefficients = vm["selection_coefficient"].as<std::vector<double>>();
    std::vector<double> dominance_coefficients = vm["dominance_coefficient"].as<std::vector<double>>();
    std::vector<double> selection_positions = vm["selection_position"].as<std::vector<double>>();
    tsk_id_t selection_activation_generation = vm["selection_activation_generation"].as<int>();
    bool mutation_in_seeds = true; // vm["mutation_in_seeds"].as<bool>();
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