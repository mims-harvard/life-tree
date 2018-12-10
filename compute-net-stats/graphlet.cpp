#include "common.h"
#include <unordered_map>
#include <vector>

using namespace std;

// Creates the necessary directories
void make_orca_dir() {
    make_dir_rooted(ORCA_FILE_DIR);
    make_dir_rooted(ORCA_STAT_DIR);
    make_dir_rooted(ORCA_LOG_DIR);
}

// Same as create_graph in analyze.cpp
// TODO merge
PGraph create_graph(int organism, string input, vector<string> &proteins) {
    string data_path = (input == "") ? get_path(EXTRACTED_DATA_DIR, organism) : input;
    ifstream data(data_path.c_str());
    string line;
    char first_protein[MAX_LENGTH];
    char second_protein[MAX_LENGTH];

    PGraph G = PGraph::TObj::New();
    // We need a map in order to retrieve the correct node for a protein
    unordered_map<string, TInt> protein_to_node;

    while (getline(data, line)) {
        sscanf(line.c_str(), "%[^ ] %s", first_protein, second_protein);

        if (protein_to_node.count(first_protein) == 0) {
            protein_to_node[first_protein] = G->AddNode();
            proteins.push_back(first_protein);
        }
        if (protein_to_node.count(second_protein) == 0) {
            protein_to_node[second_protein] = G->AddNode();
            proteins.push_back(second_protein);
        }
        G->AddEdge(protein_to_node[first_protein], protein_to_node[second_protein]);
    }

    IAssert(G->IsOk());
    return G;
}

// Generates a file compatible with the Orca graphlet-counting program
void generate_orca_file(int organism, int instance, string input, bool simple_output, vector<string> &proteins) {
    ofstream orca_file(get_output_path(ORCA_FILE_DIR, organism, simple_output).c_str());
    PGraph G = create_graph(organism, input, proteins);
    if (instance > 0) {
        TRnd Rnd(instance);
        G = TSnap::GenRewire(G, REWIRE_Q, Rnd);
    }

    orca_file << G->GetNodes() << " " << G->GetEdges() << endl;

    for (PGraph::TObj::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
        orca_file << EI.GetSrcNId() << " " << EI.GetDstNId() << endl;
    }
}

// Makes a system call to run the Orca program
void run_orca(int organism, int graphlet_size, bool simple_output) {
    // Use local executable for local runs
    string orca_exec = simple_output ? "./orca" : ORCA_EXEC_LOCATION;

    string command = orca_exec + " "
        + to_string((long long)graphlet_size) + " "
        + get_output_path(ORCA_FILE_DIR, organism, simple_output) + " "
        + get_output_path(ORCA_STAT_DIR, organism, simple_output) + "-temp > "
        + get_output_path(ORCA_LOG_DIR, organism, simple_output);
    system(command.c_str());

    // Remove edge list file
    string file_path = get_output_path(ORCA_FILE_DIR, organism, simple_output);
    system(("rm " + file_path).c_str());

    // Remove log file for local runs
    if (simple_output) {
        string log_path = get_output_path(ORCA_LOG_DIR, organism, simple_output);
        system(("rm " + log_path).c_str());
    }
}

// Prepends the protein name corresponding to each line of the Orca output
void prepend_protein_name(int organism, bool simple_output, const vector<string> &proteins) {
    string path = get_output_path(ORCA_STAT_DIR, organism, simple_output);
    ifstream old_file((path + "-temp").c_str());
    ofstream new_file(path.c_str());
    string line;

    for (int i = 0; getline(old_file, line); i++) {
        new_file << proteins.at(i) << " " << line << endl;
    }
    system(("rm " + path + "-temp").c_str());
}

// Analyzes graphlets of organism
void analyze_organism_graphlets(int organism, int instance, string input, bool simple_output, int graphlet_size) {
    vector<string> proteins;
    generate_orca_file(organism, instance, input, simple_output, proteins);
    run_orca(organism, graphlet_size, simple_output);
    prepend_protein_name(organism, simple_output, proteins);
    cout << "Finished organism " << organism << endl;
}

// Runs desired function in parallel
void run_in_parallel() {
    thread threads[NUM_ORG];
    ifstream org_list(ORG_LOCATION.c_str());
    string line;
    for (int i = 0; getline(org_list, line); i++) {
        if (i >= CONCURRENCY) threads[i - CONCURRENCY].join();

        int organism = stoi(line);
        threads[i] = thread(analyze_organism_graphlets, organism, 0, "", false, 4);
    }
    for (int i = NUM_ORG - CONCURRENCY; i < NUM_ORG; i++) {
        threads[i].join();
    }
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        make_orca_dir();
        run_in_parallel();
    }
    else if (argc > 4) {
        cout << "Use one of the following formats:" << endl
            << "[executable] to run instance 0 of all organisms in parallel" << endl
            << "[executable] [organism] to run instance 0 of [organism]" << endl
            << "[executable] [organism] [instance]" << endl
            << "[executable] [organism] [instance] [input]" << endl;
    }
    else {
        int organism = atoi(argv[1]);
        int instance = (argc >= 3) ? atoi(argv[2]) : 0;
        string input = (argc == 4) ? argv[3]: "";
        analyze_organism_graphlets(organism, instance, input, true, 4);
    }
    return 0;
}
