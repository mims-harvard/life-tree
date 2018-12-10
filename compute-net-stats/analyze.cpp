#include "common.h"
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <cmath>

using namespace std;

// Creates the necessary directories
void make_analyze_dir() {
    make_dir_rooted(STAT_DIR);
    make_dir_rooted(KCORE_DIR);
    make_dir_rooted(DEGREE_DIR);
    make_dir_rooted(CURVE_DIR);
    make_dir_rooted(CLUSTER_DIR);
    make_dir_rooted(LORENZ_DIR);
    make_dir_rooted(DEGREECTR_DIR);
    make_dir_rooted(BETWEEN_DIR);
    make_dir_rooted(CLOSE_DIR);
    //make_dir_rooted(ADJACENCY_DIR);
    make_dir_rooted(FRAGMENT_DIR);
}

// Creates a graph out of the organism's data, where vertices are proteins and edges are interactions
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

// Calculates average degree
double average_degree(PGraph G) {
    int total_degree = 0;
    int nodes = G->GetNodes();
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        total_degree += NI.GetDeg();
    }
    return (double)total_degree / nodes;
}

// Calculates maximum degree
int maximum_degree(PGraph G) {
    int max_degree = 0;
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        int curr_degree = NI.GetDeg();
        if (curr_degree > max_degree) max_degree = curr_degree;
    }
    return max_degree;
}

// Calculates density, (number of edges present) / (number of edges possible)
double density(PGraph G) {
    int edges_present = G->GetEdges();
    int nodes = G->GetNodes();
    int edges_possible = nodes * (nodes - 1) / 2;
    return (double)edges_present / edges_possible;
}

// Calculates number of connected components
int num_components(PGraph G) {
    TCnComV CnComV;
    TSnap::GetSccs(G, CnComV);
    return CnComV.Len();
}

// Calculates global clustering coefficient, (number of closed triads) / (number of triads)
double global_clustering(PGraph G) {
    int64 closed = 0, open = 0;
    TSnap::GetTriads(G, closed, open);
    return (double)closed / (closed + open);
}

// Returns n choose k
int64 choose(int n, int k) {
    int64 result = 1;
    for (int i = 0; i < k; i++) {
        // Separate into two lines to prevent truncation during integer division
        result *= n - i;
        result /= (i + 1);
    }
    return result;
}

// Calculates density of k-stars, (number of k-stars present) / (number of k-stars possible)
double star_density(PGraph G, int k) {
    int64 count = 0;
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        count += choose(NI.GetDeg(), k);
    }
    int nodes = G->GetNodes();
    // Divide in two steps to prevent overflow
    double density = (double)count / nodes;
    return density / choose(nodes, k);
}

// Calculates the Gini coefficient, a measure of degree inequality
double gini_coefficient(PGraph G) {
    TIntPrV DegToCntV;
    TSnap::GetDegCnt(G, DegToCntV);

    int nodes = G->GetNodes();
    int64 moment_degree_sum = 0; // sum of i * d_i over a sorted list of degrees d_1, d_2, ..., d_n
    int64 degree_sum = 2 * G->GetEdges();

    int nodes_so_far = 0;

    for (TIntPrV::TIter TI = DegToCntV.BegI(); TI < DegToCntV.EndI(); TI++) {
        TInt degree = TI->GetVal1();
        TInt frequency = TI->GetVal2();

        // We want to sum the series (nodes_so_far + 1) + (nodes_so_far + 2) + ... + (nodes_so_far + frequency)
        int64 series_sum = frequency * nodes_so_far + frequency * (frequency + 1) / 2;

        moment_degree_sum += series_sum * degree;
        nodes_so_far += frequency;
    }

    double result = 2 * moment_degree_sum / ((double)nodes * degree_sum);
    result -= ((double)nodes + 1) / nodes;
    return result;
}

// Calculates the normalized edge distribution entropy, a measure of degree equality
double edge_entropy(PGraph G) {
    TIntPrV DegToCntV;
    TSnap::GetDegCnt(G, DegToCntV);

    int nodes = G->GetNodes();
    int edges = G->GetEdges();
    double result = 0;

    for (TIntPrV::TIter TI = DegToCntV.BegI(); TI < DegToCntV.EndI(); TI++) {
        TInt degree = TI->GetVal1();
        TInt frequency = TI->GetVal2();
        double degree_fraction = (double)degree / (2 * edges);
        result -= frequency * degree_fraction * log(degree_fraction);
    }

    result /= log(nodes);
    return result;
}

// Calculates the Pearson correlation of the degrees at the ends of edges, a measure of assortative mixing
double assortative_mixing(PGraph G) {
    double inverse_of_edges = (double)1 / G->GetEdges();
    long long sum_of_products = 0, sum_of_sums = 0, sum_of_squares = 0;
    for (PGraph::TObj::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
        int src = G->GetNI(EI.GetSrcNId()).GetDeg();
        int dst = G->GetNI(EI.GetDstNId()).GetDeg();
        sum_of_products += src * dst;
        sum_of_sums += src + dst;
        sum_of_squares += src * src + dst * dst;
    }

    double sum_of_products_term = inverse_of_edges * sum_of_products;
    double sum_of_sums_term = inverse_of_edges * sum_of_sums / 2;
    sum_of_sums_term *= sum_of_sums_term;
    double sum_of_squares_term = inverse_of_edges * sum_of_squares / 2;
    return (sum_of_products_term - sum_of_sums_term) / (sum_of_squares_term - sum_of_sums_term);
}

// Analyzes basic summary statistics of the graph
void analyze_statistics(PGraph G, int organism, bool simple_output) {
    string output_path = get_output_path(STAT_DIR, organism, simple_output);
    //ofstream statistics(output_path.c_str(), ios_base::app);
    ofstream statistics(output_path.c_str());

    statistics << fixed << setprecision(PRECISION)
        << G->GetNodes() << endl
        << G->GetEdges() << endl
        << average_degree(G) << endl
        << maximum_degree(G) << endl
        << density(G) << endl
        << num_components(G) << endl
        << TSnap::GetMxSccSz(G) << endl
        << TSnap::GetBfsFullDiam(G, 100) << endl
        << TSnap::GetBfsEffDiam(G, 100) << endl
        << global_clustering(G) << endl
        << TSnap::GetClustCf(G) << endl
        << star_density(G, 2) << endl
        << star_density(G, 3) << endl
        << gini_coefficient(G) << endl
        << edge_entropy(G) << endl
        << assortative_mixing(G) << endl;
}

// Analyzes the k-core distribution of the graph
void analyze_kcores(PGraph G, int organism, bool simple_output) {
    string output_path = get_output_path(KCORE_DIR, organism, simple_output);
    //ofstream kcores(output_path.c_str(), ios_base::app);
    ofstream kcores(output_path.c_str());
    TKCore<PGraph> KCore(G);

    while (KCore.GetNextCore()) {
        kcores << KCore.GetCurK() << " "
            << KCore.GetCoreNodes() << " "
            << KCore.GetCoreEdges() << endl;
    }
}

// Analyzes the degree distribution of the graph
void analyze_degrees(PGraph G, int organism, bool simple_output) {
    string output_path = get_output_path(DEGREE_DIR, organism, simple_output);
    //ofstream degrees(output_path.c_str(), ios_base::app);
    ofstream degrees(output_path.c_str());
    TIntPrV DegToCntV;
    TSnap::GetDegCnt(G, DegToCntV);

    for (TIntPrV::TIter TI = DegToCntV.BegI(); TI < DegToCntV.EndI(); TI++) {
        degrees << TI->GetVal1() << " "
            << TI->GetVal2() << endl;
    }
}

// Analyzes the node curvature distribution of the graph
void analyze_curvature(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins) {
    string output_path = get_output_path(CURVE_DIR, organism, simple_output);
    //ofstream curvature(output_path.c_str(), ios_base::app);
    ofstream curvature(output_path.c_str());

    curvature << fixed << setprecision(PRECISION);

    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        int NId = NI.GetId();
        double node_curvature = 1 - (double)NI.GetDeg() / 2 + (double)TSnap::GetNodeTriads(G, NId) / 3;
        curvature << proteins.at(NId) << " "
            << node_curvature << endl;

    }
}

// Analyzes the clustering coefficient distribution of the graph
void analyze_clustering(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins) {
    string output_path = get_output_path(CLUSTER_DIR, organism, simple_output);
    //ofstream clustering(output_path.c_str(), ios_base::app);
    ofstream clustering(output_path.c_str());
    TIntFltH NIdCCfH;
    TSnap::GetNodeClustCf(G, NIdCCfH);

    clustering << fixed << setprecision(PRECISION);

    for (TIntFltH::TIter TI = NIdCCfH.BegI(); TI < NIdCCfH.EndI(); TI++) {
        clustering << proteins.at(TI.GetKey()) << " "
            << TI.GetDat() << endl;
    }
}

// Analyzes the degree inequality of the graph with Lorenz curves, using units of parts per thousand
void analyze_lorenz(PGraph G, int organism, bool simple_output) {
    static const int NUM_POINTS = 50;
    string input_path = get_output_path(DEGREE_DIR, organism, simple_output);
    string output_path = get_output_path(LORENZ_DIR, organism, simple_output);
    ifstream degrees(input_path.c_str(), ios_base::app);
    //ofstream lorenz(output_path.c_str(),f ios_base::app);
    ofstream lorenz(output_path.c_str());

    int64 total_degree_sum = 2 * G->GetEdges();
    int64 total_nodes = G->GetNodes();

    lorenz << fixed << setprecision(PRECISION);

    for (int i = 1; i <= NUM_POINTS; i++) {
        int64 degree_sum = 0;
        int64 nodes = total_nodes * i / NUM_POINTS;
        string line;
        int degree, frequency;

        degrees.clear();
        degrees.seekg(0);

        while (getline(degrees, line)) {
            sscanf(line.c_str(), "%d %d", &degree, &frequency);
            if (nodes > frequency) {
                nodes -= frequency;
                degree_sum += frequency * degree;
            }
            else {
                degree_sum += nodes * degree;
                break;
            }
        }

        lorenz << (double)i / NUM_POINTS << " "
            << (double)degree_sum / total_degree_sum << endl;
    }
}

// Analyzes the degree centrality distribution of the graph
void analyze_degreecentrality(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins, TIntFltH &NIdDegH) {
    string output_path = get_output_path(DEGREECTR_DIR, organism, simple_output);
    //ofstream degreecentrality(output_path.c_str(), ios_base::app);
    ofstream degreecentrality(output_path.c_str());

    degreecentrality << fixed << setprecision(PRECISION);

    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        int NId = NI.GetId();
        double Deg = TSnap::GetDegreeCentr(G, NId);
        NIdDegH.AddDat(NId, Deg);
        degreecentrality << proteins.at(NId) << " "
            << Deg << endl;
    }
}

// Analyzes the normalized node betweenness centrality distribution of the graph
void analyze_betweenness(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins, TIntFltH &NIdBtwH) {
    string output_path = get_output_path(BETWEEN_DIR, organism, simple_output);
    //ofstream betweenness(output_path.c_str(), ios_base::app);
    ofstream betweenness(output_path.c_str());
    TSnap::GetBetweennessCentr(G, NIdBtwH);
    int nodes = G->GetNodes();
    long long max_nodes = ((long long)nodes - 1) * (nodes - 2) / 2;

    betweenness << fixed << setprecision(PRECISION);

    for (TIntFltH::TIter TI = NIdBtwH.BegI(); TI < NIdBtwH.EndI(); TI++) {
        betweenness << proteins.at(TI.GetKey()) << " "
            << TI.GetDat() / max_nodes << endl;
    }
}

// Analyzes the closeness centrality distribution of the graph
void analyze_closeness(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins, TIntFltH &NIdCloH) {
    string output_path = get_output_path(CLOSE_DIR, organism, simple_output);
    //ofstream closeness(output_path.c_str(), ios_base::app);
    ofstream closeness(output_path.c_str());

    closeness << fixed << setprecision(PRECISION);

    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        int NId = NI.GetId();
        double Clo = TSnap::GetClosenessCentr<PGraph>(G, NId, false);
        NIdCloH.AddDat(NId, Clo);
        closeness << proteins.at(NId) << " "
            << Clo << endl;
    }
}

// Analyzes the eigenvalues of the adjacency matrix of the graph
// NOTE: TSnap's GetEigVec and GetEigVals currently do not work due to numerical instability.
void analyze_adjacency(PGraph G, int organism, bool simple_output) {
    string output_path = get_output_path(ADJACENCY_DIR, organism, simple_output);
    //ofstream adjacency(output_path.c_str(), ios_base::app);
    ofstream adjacency(output_path.c_str());
    TFltV EigValV;
    TVec<TFltV> EigVecV;
    int nodes = G->GetNodes();
    TSnap::GetEigVec(G, nodes, EigValV, EigVecV);

    adjacency << fixed << setprecision(PRECISION);

    for (TFltV::TIter TI = EigValV.BegI(); TI < EigValV.EndI(); TI++) {
        adjacency << *TI << endl;
    }
}

// Helper that analyzes fragmentation for a single node removal strategy
// Calculates distribution of strongly connected component sizes
    void node_removal(PGraph G, int organism, string strategy, bool simple_output, TIntFltH &Hash) {
    static const int NUM_POINTS = 100;
    string output_path = get_output_path(FRAGMENT_DIR, organism, simple_output);
    output_path += strategy;
    //ofstream fragmentation(output_path.c_str(), ios_base::app);
    ofstream fragmentation(output_path.c_str());

    fragmentation << fixed << setprecision(PRECISION);

    int nodes = G->GetNodes();

    Hash.SortByDat(false);
    TIntFltH::TIter TI = Hash.BegI();
    int num_deleted = 0;
    PGraph H = new PGraph::TObj(*G);
    int current_nodes;

    for (int i = 1; i <= NUM_POINTS; i++) {
        for ( ; num_deleted < nodes * i / NUM_POINTS; TI++, num_deleted++) {
            H->DelNode(TI.GetKey());
        }
        TIntPrV SccSzCnt;
        TSnap::GetSccSzCnt(H, SccSzCnt);
        current_nodes = H->GetNodes();
        fragmentation << (double)i / NUM_POINTS << " " << current_nodes << ", ";
        for (TIntPrV::TIter TI1 = SccSzCnt.BegI(); TI1 < SccSzCnt.EndI(); TI1++) {
            TInt size = TI1->GetVal1();
            TInt frequency = TI1->GetVal2();
            fragmentation << size << " " << frequency << ", ";
        }
        fragmentation << endl;
    }
}

// Analyzes the fragmentation of the graph to node and edge removal
void analyze_fragmentation(PGraph G, int organism, bool simple_output,
        TIntFltH &NIdDegH, TIntFltH &NIdBtwH, TIntFltH &NIdCloH) {
    TIntFltH NIdRndH;
    TRnd random(0); // 0 means seed from clock
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        NIdRndH.AddDat(NI.GetId(), random.GetUniDev());
    }

    node_removal(G, organism, "random", simple_output, NIdRndH); // remove random
    node_removal(G, organism, "degree", simple_output, NIdDegH); // remove highest degree
}

// Creates graph of organism's data and analyzes it
void analyze_organism(int organism, string input, bool simple_output) {
    cout << "Started organism " << organism << endl;
    vector<string> proteins;
    TIntFltH NIdDegH, NIdBtwH, NIdCloH;
    PGraph G = create_graph(organism, input, proteins);

    // analyze_statistics(G, organism, simple_output);
    G = TSnap::GetMxScc(G);
    analyze_statistics(G, organism, simple_output); // append lcc statistics to end of file
    analyze_kcores(G, organism, simple_output);
    analyze_degrees(G, organism, simple_output);
    analyze_curvature(G, organism, simple_output, proteins);
    analyze_clustering(G, organism, simple_output, proteins);
    analyze_lorenz(G, organism, simple_output);
    analyze_degreecentrality(G, organism, simple_output, proteins, NIdDegH);
//    analyze_betweenness(G, organism, simple_output, proteins, NIdBtwH);
//    analyze_closeness(G, organism, simple_output, proteins, NIdCloH);
//    analyze_adjacency(G, organism, simple_output);
    analyze_fragmentation(G, organism, simple_output, NIdDegH, NIdBtwH, NIdCloH);
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
        threads[i] = thread(analyze_organism, organism, "", false);
    }
    for (int i = NUM_ORG - CONCURRENCY; i < NUM_ORG; i++) {
        threads[i].join();
    }
}

// Generates list of commands
void generate_job_list(int num_instances) {
    ifstream org_list(ORG_LOCATION.c_str());
    ofstream job_list(JOB_LIST_LOCATION.c_str());
    string line;

    while(getline(org_list, line)) {
        for (int j = 0; j <= num_instances; j++) {
            job_list << "./analyze " << line << " " << j << " "
                << EXTRACTED_DATA_DIR << "/" << line << endl;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        make_analyze_dir();
        run_in_parallel();
    }
    else if (argc > 3) {
        cout << "Use one of the following formats:" << endl
            << "[executable] to run all organisms in parallel" << endl
            << "[executable] [organism] to run [organism]" << endl
            << "[executable] [organism]" << endl
            << "[executable] [organism] [input]" << endl;
    }
    else {
        int organism = atoi(argv[1]);
        string input = (argc == 3) ? argv[2]: "";
        analyze_organism(organism, input, false);
    }
    return 0;
}
