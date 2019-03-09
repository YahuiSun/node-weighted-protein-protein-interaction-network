#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <boost/config.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <chrono>
#include <boost/heap/pairing_heap.hpp> // pairing_heap uses less memory
#include <tuple>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <typeinfo>
#include <stdlib.h>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::numeric::ublas; // there is a vector file inside which conflicts with std; so use std::vector to specify the vector definition
using namespace std;
using namespace boost::heap;

#pragma region 
// define an adjacency list with edge weights
typedef boost::property<boost::edge_weight_t, float> EdgeWeightProperty; // define edge weight property
typedef boost::property<boost::vertex_name_t, float> VertexWeightProperty; // define node weight property; note that: vertex_index_t is not mutable
typedef boost::adjacency_list<boost::setS, boost::vecS,
	boost::undirectedS, VertexWeightProperty, EdgeWeightProperty> graph; // define all the graph properties
typedef boost::graph_traits<graph>::adjacency_iterator AdjacencyIterator;
#pragma endregion define graph property 2018年3月23日20:15:27


#pragma region
struct node {
	int index;
	double priority_value;
}; // define the node in the queue
bool operator<(node const& x, node const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the mean heap; PriorityQueue is expected to be a max-heap of integer values
} // redefine the operator since there are multiple values in each node
typedef typename pairing_heap<node>::handle_type handle_t; // define the handle type for pairing_heap<node>
#pragma endregion define heaps


#pragma region 
graph read_PPI_data(string file_name) {

	string file_ID = file_name + ".txt";
	int v1;
	int v2;
	double cost;
	double experimental;
	double database;
	double new_cost;
	graph input_graph; // define the adjacency list of the input graph;
	string line_content;
	ifstream myfile(file_ID); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int line_num = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			line_num++;
			if (line_num != 1) {
				// parse the sting：line_content
				std::vector<string> Parsed_content;
				std::string delimiter = " "; // the delimiter
				size_t pos = 0;
				std::string token;
				while ((pos = line_content.find(delimiter)) != std::string::npos) {
					// find(const string& str, size_t pos = 0) function returns the position of the first occurrence of str in the string, or npos if the string is not found.
					token = line_content.substr(0, pos);
					// The substr(size_t pos = 0, size_t n = npos) function returns a substring of the object, starting at position pos and of length npos
					Parsed_content.push_back(token); // store the subtr to the list
					line_content.erase(0, pos + delimiter.length()); // remove the front substr and the first delimiter
				}
				Parsed_content.push_back(line_content); // store the subtr to the list

				// split Parsed_content[0] by "_"
				delimiter = "P"; // the delimiter
				pos = 0;
				while ((pos = Parsed_content[0].find(delimiter)) != std::string::npos) {
					token = line_content.substr(0, pos);
					Parsed_content[0].erase(0, pos + delimiter.length()); // remove the front substr and the first delimiter
				}
				// split Parsed_content[1] by "_"
				pos = 0;
				while ((pos = Parsed_content[1].find(delimiter)) != std::string::npos) {
					token = line_content.substr(0, pos);
					Parsed_content[1].erase(0, pos + delimiter.length()); // remove the front substr and the first delimiter
				}
				v1 = stoi(Parsed_content[0]); // protein 1
				v2 = stoi(Parsed_content[1]); // protein 2
				experimental = stod(Parsed_content[6]);
				database = stod(Parsed_content[6]);
				cost = sqrt(experimental*experimental);
				if (cost > 0) {
					boost::add_edge(v1, v2, cost, input_graph); // add edge; edge cost: combined score
				}

				//for (int i = 0; i < Parsed_content.size(); i++) {
				//	cout << Parsed_content[i] << " ";
				//}
				//getchar();

			}
		}

		myfile.close(); //close the file

		return input_graph;
	}
	else
	{
		std::cout << "Unable to open file " << file_ID << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}

}
#pragma endregion read_PPI_data 2018年3月23日20:14:52


#pragma region
std::vector<std::vector<int>> read_csv_to_path(string staynerd_dimacs_file) {
	std::vector<std::vector<string>> input;

	string line_content;
	ifstream myfile(staynerd_dimacs_file); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		while (getline(myfile, line_content)) // read file line by line
		{
			// parse the sting：line_content
			list<string> Parsed_content;
			std::string delimiter = ","; // the delimiter
			size_t pos = 0;
			std::string token;
			while ((pos = line_content.find(delimiter)) != std::string::npos) {
				// find(const string& str, size_t pos = 0) function returns the position of the first occurrence of str in the string, or npos if the string is not found.
				token = line_content.substr(0, pos);
				// The substr(size_t pos = 0, size_t n = npos) function returns a substring of the object, starting at position pos and of length npos
				Parsed_content.push_back(token); // store the subtr to the list
												 //cout << token << std::endl; // print the front substr
				line_content.erase(0, pos + delimiter.length()); // remove the front substr and the first delimiter
			}
			//std::cout << line_content << std::endl; // this is the last substr in this line
			Parsed_content.push_back(line_content); // store the subtr to the list

			std::vector<string> line(Parsed_content.size());
			while (Parsed_content.size() > 0) {
				line.insert(line.end(), Parsed_content.front().c_str());
				Parsed_content.pop_front();
			}
			input.insert(input.end(), line);
		}
		myfile.close(); //close the file

		std::vector<std::vector<int>> output;
		for (int i = 0; i < input.size(); i++) {
			std::vector<int> new_line;
			for (int j = 0; j < input[i].size(); j++) {
				// split input[i][j] by "P"
				std::string delimiter = "P"; // the delimiter
				size_t pos = 0;
				std::string token;
				while ((pos = input[i][j].find(delimiter)) != std::string::npos) {
					token = line_content.substr(0, pos);
					input[i][j].erase(0, pos + delimiter.length()); // remove the front substr and the first delimiter
				}
				if (input[i][j].size() > 0) { // empty may be inside
					new_line.insert(new_line.end(), stoi(input[i][j]));
				}
			}
			output.insert(output.end(), new_line);
		}

		std::vector<int> edge;
		std::vector<std::vector<int>> path_edges;
		for (int i = 0; i < output.size() - 1; i++) {
			for (int j = 0; j < output[i].size(); j++) {
				for (int k = 0; k < output[i + 1].size(); k++) {
					int v1 = output[i][j];
					int v2 = output[i + 1][k];
					edge.clear();
					edge.insert(edge.end(), v1);
					edge.insert(edge.end(), v2);
					path_edges.insert(path_edges.end(), edge);
				}
			}
		}

		return path_edges;
	}
	else
	{
		std::cout << "Unable to open file " << staynerd_dimacs_file << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}
}

#pragma endregion  read_csv_to_path


#pragma region 

std::vector<std::vector<int>> read_path() {

	std::vector<std::vector<int>> path_edges;

	std::vector<std::vector<int>> sub_path = read_csv_to_path("path_edges_1.csv");
	for (int i = 0; i < sub_path.size(); i++) {
		path_edges.insert(path_edges.end(), sub_path[i]);
	}
	sub_path = read_csv_to_path("path_edges_2.csv");
	for (int i = 0; i < sub_path.size(); i++) {
		path_edges.insert(path_edges.end(), sub_path[i]);
	}
	sub_path = read_csv_to_path("path_edges_3.csv");
	for (int i = 0; i < sub_path.size(); i++) {
		path_edges.insert(path_edges.end(), sub_path[i]);
	}
	sub_path = read_csv_to_path("path_edges_4.csv");
	for (int i = 0; i < sub_path.size(); i++) {
		path_edges.insert(path_edges.end(), sub_path[i]);
	}
	sub_path = read_csv_to_path("path_edges_5.csv");
	for (int i = 0; i < sub_path.size(); i++) {
		path_edges.insert(path_edges.end(), sub_path[i]);
	}
	sub_path = read_csv_to_path("path_edges_6.csv");
	for (int i = 0; i < sub_path.size(); i++) {
		path_edges.insert(path_edges.end(), sub_path[i]);
	}
	sub_path = read_csv_to_path("path_edges_7.csv");
	for (int i = 0; i < sub_path.size(); i++) {
		path_edges.insert(path_edges.end(), sub_path[i]);
	}


	ofstream outputFile;
	outputFile.open("path_edges.csv"); // stp file
	for (int i = 0; i < path_edges.size(); i++) {
		outputFile << path_edges[i][0] << "," << path_edges[i][1] << endl;
	}


	return path_edges;
}

#pragma endregion read_path


#pragma region
/* Matrix inversion routine.
Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix(const matrix<T>& input, matrix<T>& inverse)
{
	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix<T> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<T>(A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}
void AB_equal_C_toObtain_B(matrix<double> A, matrix<double>& B, matrix<double> C) {
	matrix<double> A_inverse;
	A_inverse.resize(A.size2(), A.size1());
	InvertMatrix(A, A_inverse);
	B.resize(A.size2(), C.size2());
	axpy_prod(A_inverse, C, B);
}
graph PO_node_weight_edge_cost(graph g, std::vector<int> terminal) {

	int N = num_vertices(g); // the total number of vertices in the initial graph
	int terminal_num = 0; // total number of terminals
	for (int i = 0; i < N; i++) {
		if (terminal[i] == 1) {
			terminal_num++;
		}
	}
	double max_node_weight = 0;
	for (int i = 0; i < N; i++) {
		if (max_node_weight < get(boost::vertex_name_t(), g, i)) {
			max_node_weight = get(boost::vertex_name_t(), g, i);
		}
	}

	// parameters in PSIA
	double I = 1;
	double Inner_iteration_time = 1;  // iteration times to identify each subnetwork

	// flux matrix
	graph total_flux_edge(N);
	for (int i = 0; i < N; i++) {
		graph::out_edge_iterator eit, eend;
		tie(eit, eend) = boost::out_edges(i, g); // adjacent_vertices of i
		for_each(eit, eend,
			[&g, &total_flux_edge, &i](graph::edge_descriptor it)
		{ int j = boost::target(it, g);
		if (i < j) {
			boost::add_edge(i, j, 0, total_flux_edge);
		}
		});
	}
	std::vector<double> total_flux_vertex;
	total_flux_vertex.resize(N);

	// the outer iteration
	for (int sink_terminal_order = 1; sink_terminal_order < terminal_num + 1; sink_terminal_order++) {

		std::vector<int> exist(N); // 1 means the vertex is in the subnetwork; 0 means not

		// the inner iteration
		for (int inner = 0; inner < Inner_iteration_time; inner++) {

			int included_vertex_num = 0;
			for (int i = 0; i < N; i++) {
				if (in_degree(i, g) > 0) {
					exist[i] = 1;
					included_vertex_num++;
				}
				else {
					exist[i] = 0;
				}
			}

			// choose a sink node
			int sink_vertex_order;
			int j = 1;
			for (int i = 0; i < N; i++) {
				if (terminal[i] == 1) {
					if (j == sink_terminal_order) {
						sink_vertex_order = i;
						break;
					}
					else {
						j++;
					}
				}
			}


			// calculate the pressures using the network Poisson equation
			std::vector<int> vertex2row(N);
			matrix<double> A;
			A.resize(included_vertex_num - 1, included_vertex_num - 1); // size: (included_vertex_num - 1, included_vertex_num - 1)
			int row_num = 0;
			for (int i = 0; i < N; i++) {
				if (exist[i] == 1) { // only calculate the existing vertice
					if (i != sink_vertex_order) { //  neglect the vertex if it's the sink point
						int column_num = 0;
						vertex2row[i] = row_num; // record the row number's relationship with vertex number, which will be used when using pressure
						for (int j = 0; j < N; j++) {
							if (exist[j] == 1) { // only calculate the existing vertice
								if (j != sink_vertex_order) { //  neglect the vertex if it's the sink point
									if (i == j) { // the row vertex is the colume vertex
										graph::out_edge_iterator eit, eend;
										tie(eit, eend) = boost::out_edges(i, g); // adjacent_vertices of i
										for_each(eit, eend,
											[&g, &A, &row_num, &column_num, &i, &max_node_weight](graph::edge_descriptor it)
										{ int adjacent = boost::target(it, g);
										A(row_num, column_num) = A(row_num, column_num) - 1 /
											(get(boost::edge_weight_t(), g, boost::edge(i, adjacent, g).first) -
												get(boost::vertex_name_t(), g, i) / in_degree(i, g) -
												get(boost::vertex_name_t(), g, adjacent) / in_degree(adjacent, g) + 2 * max_node_weight);
										});
									}
									else { // the row vertex is not the colume vertex
										typedef graph::edge_descriptor Edge;
										pair<Edge, bool> ed = boost::edge(i, j, g);
										if (ed.second) { // edge i,j exists
											A(row_num, column_num) = 1 / get(boost::edge_weight_t(), g, boost::edge(i, j, g).first) -
												get(boost::vertex_name_t(), g, i) / in_degree(i, g) -
												get(boost::vertex_name_t(), g, j) / in_degree(j, g) + 2 * max_node_weight;
										}
									}
									column_num++;
								}
							}
						}
						row_num++;
					}
				}
			}
			matrix<double> X;
			X.resize(included_vertex_num - 1, 1);
			for (int i = 0; i < N; i++) {
				if (terminal[i] == 1 & i != sink_vertex_order) {
					X(vertex2row[i], 0) = -I;
				}
			}
			matrix<double> P;

			//cout << vertex2row[0] << "," << vertex2row[1] << "," << vertex2row[2] << endl;
			//cout << A(0, 0) << "," << A(0, 1) << "," << A(1, 0) << "," << A(1, 1) << endl;
			//cout << X(0, 0) << "," << X(1, 0) << endl;
			//cout << P(0, 0) << "," << P(1, 0) << endl;
			//getchar();
			AB_equal_C_toObtain_B(A, P, X);

			// sequence pressures
			std::vector<double> pressure;
			pressure.resize(N);
			for (int i = 0; i < N; i++) {
				if (exist[i] == 1) { // only existing vertices have pressures
					if (i != sink_vertex_order) { // sink's pressure is 0
						pressure[i] = P(vertex2row[i], 0);
					}
					else {
						pressure[i] = 0;
					}
				}
				else {
					pressure[i] = 0;
				}
			}

			// accumulate fluxes
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < i; j++) {
					typedef graph::edge_descriptor Edge;
					pair<Edge, bool> ed = boost::edge(i, j, g);
					if (ed.second) {
						double flux = abs(1 * (pressure[i] - pressure[j]) / (get(boost::edge_weight_t(), g, boost::edge(i, j, g).first) -
							get(boost::vertex_name_t(), g, i) / in_degree(i, g) -
							get(boost::vertex_name_t(), g, j) / in_degree(j, g) + 2 * max_node_weight));
						typedef graph::edge_descriptor Edge;
						pair<Edge, bool> ed = boost::edge(i, j, total_flux_edge);
						boost::put(boost::edge_weight_t(), total_flux_edge, ed.first, get(boost::edge_weight_t(), total_flux_edge, boost::edge(i, j, total_flux_edge).first) + flux);
						total_flux_vertex[i] = total_flux_vertex[i] + flux;
						total_flux_vertex[j] = total_flux_vertex[j] + flux;
						cout << total_flux_vertex[i] << total_flux_vertex[j] << endl;
					}
				}
			}
		}
	}

	return total_flux_edge;
}
#pragma endregion PO_related  2017年10月10日12:42:53


#pragma region 
std::vector<int> SP_Yahui(graph g, int source, int destionation, double& SP_distance) {

	// Create things for Dijkstra
	std::vector<int> predecessors(boost::num_vertices(g)); // To store parents
	std::vector<double> distances(boost::num_vertices(g)); // To store distances
	typedef boost::property_map < graph, boost::vertex_index_t >::type IndexMap;
	IndexMap indexMap = boost::get(boost::vertex_index, g);
	typedef boost::iterator_property_map < int*, IndexMap, int, int& > PredecessorMap;
	PredecessorMap predecessorMap(&predecessors[0], indexMap);
	typedef boost::iterator_property_map < double*, IndexMap, double, double& > DistanceMap;
	DistanceMap distanceMap(&distances[0], indexMap);

	// Compute shortest paths from source to all vertices, and store the output in predecessors and distances
	boost::dijkstra_shortest_paths(g, source, boost::distance_map(distanceMap).predecessor_map(predecessorMap));

	// Extract a shortest path
	std::vector<int> path;
	path.insert(path.begin(), destionation);
	int v = destionation; // We want to start at the destination and work our way back to the source
	for (int u = predecessorMap[v]; // Start by setting 'u' to the destintaion node's predecessor
		u != v; // Keep tracking the path until we get to the source
		v = u, u = predecessorMap[v]) // Set the current vertex to the current predecessor, and the predecessor to one level up
	{
		path.insert(path.begin(), u);
	}

	SP_distance = distanceMap[destionation];
	return path;
}
#pragma endregion SP_Yahui


#pragma region
void source_terminal(std::vector<int>& path_source, std::vector<int>& path_terminal) {
	// source and terminal;
	path_source.insert(path_source.end(), 393312);
	path_source.insert(path_source.end(), 257290);
	path_source.insert(path_source.end(), 275493);
	path_source.insert(path_source.end(), 269571);
	path_source.insert(path_source.end(), 410294);
	path_source.insert(path_source.end(), 357178);
	path_source.insert(path_source.end(), 261799);
	path_source.insert(path_source.end(), 268035);
	// for the first big path
	path_terminal.insert(path_terminal.end(), 330237);
	path_terminal.insert(path_terminal.end(), 309103);
	path_terminal.insert(path_terminal.end(), 368880);
	path_terminal.insert(path_terminal.end(), 244741);
	path_terminal.insert(path_terminal.end(), 228872);
	path_terminal.insert(path_terminal.end(), 354558);
	path_terminal.insert(path_terminal.end(), 226574);
	path_terminal.insert(path_terminal.end(), 384273);
	path_terminal.insert(path_terminal.end(), 269305);
	path_terminal.insert(path_terminal.end(), 340347);
	path_terminal.insert(path_terminal.end(), 282111);
	path_terminal.insert(path_terminal.end(), 444972);
	path_terminal.insert(path_terminal.end(), 265165);
	// for the second big path
	path_terminal.insert(path_terminal.end(), 363822);
}
#pragma endregion source_terminal


#pragma region
graph customizing_graph_single_complex_path(graph input_graph, std::vector<int>& initial_2_new, std::vector<int>& new_2_initial, std::vector<std::vector<int>> path,
	std::vector<int> path_source, std::vector<int> path_terminal) {

	cout << "|V| before = " << num_vertices(input_graph) << endl;
	cout << "|E| before = " << num_edges(input_graph) << endl;

	int N = num_vertices(input_graph);
	std::vector<int> component(N); // vertex i is in component[i]; No.component from 0
	int cpn_num = connected_components(input_graph, &component[0]); // the number of component; decrease component

	// component sizes
	std::vector<int> cpn_size(cpn_num);
	for (int i = 0; i < N; i++) {
		cpn_size[component[i]] = cpn_size[component[i]] + 1;
	}

	// the No. of the biggest cpn
	int max_cpn = 0;
	int max_size = cpn_size[0];
	for (int i = 1; i < cpn_num; i++) {
		if (cpn_size[i] > max_size) {
			max_size = cpn_size[i];
			max_cpn = i;
		}
	}

	//// old terminal
	std::vector<int> old_terminal;
	old_terminal.resize(N);
	//for (int i = 0; i < path.size(); i++) {
	//	for (int j = 0; j < path[i].size(); j++) {
	//		old_terminal[path[i][j]] = 1;
	//	}
	//}
	for (int i = 0; i < path_source.size(); i++) {
		old_terminal[path_source[i]] = 1;
	}
	for (int i = 0; i < path_terminal.size(); i++) {
		old_terminal[path_terminal[i]] = 1;
	}



	// check existence of edges in the single complex path
	for (int i = 0; i < path.size(); i++) {
		int v1 = path[i][0];
		int v2 = path[i][1];
		typedef graph::edge_descriptor Edge;
		pair<Edge, bool> ed = boost::edge(v1, v2, input_graph);
		if (ed.second) {
			double combined_score = get(boost::edge_weight_t(), input_graph, boost::edge(v1, v2, input_graph).first);
			boost::put(boost::edge_weight_t(), input_graph, ed.first, combined_score * 2); // give edges in the pathways big scores
		}
		//else {
		//	cout << "Edge " << v1 << "," << v2 << " does not exist!" << endl;
		//	//getchar();
		//}
	}

	// output a component
	int target_vertex; // the interested vertex; note that, this terminal number is the initial protein number
	for (int i = 0; i < N; i++) {
		if (old_terminal[i] == 1) {
			target_vertex = i;
			break;
		}
	}
	int target_cpn = component[target_vertex];

	//// check terminals are inside
	//for (int i = 0; i < N; i++) {
	//	if (old_terminal[i] == 1 && component[i] != target_cpn) {
	//		cout << "Error! Terminal not inside!" << endl;
	//		getchar();
	//	}
	//}


	auto begin_time = std::chrono::high_resolution_clock::now(); // start time
	// removde degree 1 non-terminals
	int change = 1;
	while (change == 1) {
		change = 0;
		for (int i = 0; i < N; i++) {
			if (component[i] == target_cpn && in_degree(i, input_graph) == 1 && old_terminal[i] == 0) {
				component[i] = N + 1; // remove vertex
				cpn_size[target_cpn]--;
				graph::out_edge_iterator eit, eend;
				tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of i
				for_each(eit, eend,
					[&input_graph, &i](graph::edge_descriptor it)
				{ int j = boost::target(it, input_graph);
				boost::remove_edge(i, j, input_graph); // remove edge
				});
				change = 1;
			}
		}
	}
	for (int i = 0; i < N; i++) {
		if (old_terminal[i] == 1 && in_degree(i, input_graph) == 1) {
			cout << "For Terminal degree 1 test" << endl;
		}
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count(); // Nanosecond
	cout << "RT time:" << runningtime / 1e6 << " ms" << endl;


	//// remove faraway vertices
	//std::vector<int> initial_terminal_number;
	//initial_terminal_number = { 1, 2, 3, 4, 5 }; // select terminals manually;  this is the initial protein number
	//int upperbound = 10;
	//for (int i = 0; i < N; i++) {
	//	graph::out_edge_iterator eit, eend;
	//	tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of i
	//	for_each(eit, eend,
	//		[&input_graph, &i](graph::edge_descriptor it)
	//	{ int j = boost::target(it, input_graph);
	//	if (i < j) {
	//		typedef graph::edge_descriptor Edge;
	//		pair<Edge, bool> ed = boost::edge(i, j, input_graph);
	//		boost::put(boost::edge_weight_t(), input_graph, ed.first, 1); // change all the edge cost to 1
	//	}
	//	});
	//}
	//for (int i = 0; i < N; i++) {
	//	if (component[i] == target_cpn) {
	//		double max_distance = 0;
	//		for (int j = 0; j < initial_terminal_number.size(); j++) {
	//			// Create things for Dijkstra
	//			std::vector<int> predecessors(boost::num_vertices(input_graph)); // To store parents
	//			std::vector<double> distances(boost::num_vertices(input_graph)); // To store distances
	//			typedef boost::property_map < graph, boost::vertex_index_t >::type IndexMap;
	//			IndexMap indexMap = boost::get(boost::vertex_index, input_graph);
	//			typedef boost::iterator_property_map < int*, IndexMap, int, int& > PredecessorMap;
	//			PredecessorMap predecessorMap(&predecessors[0], indexMap);
	//			typedef boost::iterator_property_map < double*, IndexMap, double, double& > DistanceMap;
	//			DistanceMap distanceMap(&distances[0], indexMap);
	//			// Compute shortest paths from source to all vertices, and store the output in predecessors and distances
	//			boost::dijkstra_shortest_paths(input_graph, i, boost::distance_map(distanceMap).predecessor_map(predecessorMap));
	//			if (distanceMap[initial_terminal_number[j]] > max_distance) {
	//				max_distance = distanceMap[j]; // the max SP to terminals
	//			}
	//		}
	//		if (max_distance > upperbound) {
	//			component[i] = N + 1; // remove the farawat vertex
	//			cpn_size[target_cpn]--;
	//		}
	//	}
	//}
	//cout << "get there! " << cpn_size[target_cpn] << endl;

	auto generation_begin_time = std::chrono::high_resolution_clock::now(); // start time

	graph output_graph(cpn_size[target_cpn]);
	initial_2_new.resize(N);
	new_2_initial.resize(cpn_size[target_cpn]);
	int x = 0;
	for (int i = 0; i < N; i++) {
		if (component[i] == target_cpn) {
			initial_2_new[i] = x; // i in input_graph is initial_2_new[i] in output_graph
			new_2_initial[x] = i;
			x++;
		}
	}

	for (int i = 0; i < N; i++) {
		if (component[i] == target_cpn) {
			graph::out_edge_iterator eit, eend;
			tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of i
			for_each(eit, eend,
				[&input_graph, &i, &output_graph, &initial_2_new](graph::edge_descriptor it)
			{
				int j = boost::target(it, input_graph);
				if (j > i) {
					// add edge and the initial edge cost
					double combined_score = get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first);
					double cost = 2e6 / pow(combined_score, 2); // initially edge costs are the reverse of combined_score
					boost::add_edge(initial_2_new[i], initial_2_new[j], cost, output_graph);
				}
			});
		}
	}

	// add node weight
	N = num_vertices(output_graph);
	std::vector<int> terminal;
	terminal.resize(N);
	for (int i = 0; i < N; i++) {
		double degree = in_degree(i, output_graph); // degree
		double weight = -5e0 / degree; // big negatives help to purefy the subnetworks
		if (old_terminal[new_2_initial[i]]) { // note that, this terminal number is the new vertex number, not the initial protein number
			terminal[i] = 1;
			boost::put(boost::vertex_name_t(), output_graph, i, 1e8); // give terminals big node weights
		}
		else {
			terminal[i] = 0;
			boost::put(boost::vertex_name_t(), output_graph, i, weight);
		}

	}

	auto generation_end_time = std::chrono::high_resolution_clock::now();
	double generation_runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(generation_end_time - generation_begin_time).count(); // Nanosecond
	cout << "generation time:" << generation_runningtime / 1e6 << " ms" << endl;

	// add PO edge costs and node weights (flux)
	// PO_node_weight_edge_cost(output_graph, terminal);

	cout << "|V| after = " << num_vertices(output_graph) << endl;
	cout << "|E| after = " << num_edges(output_graph) << endl;

	return output_graph;
}
#pragma endregion customizing_graph_single_complex_path


#pragma region
void save_data(string instance_name, graph result_graph, std::vector<int> new_2_initial) {

	string save_name = instance_name; // save_name
	ofstream outputFile;
	outputFile.precision(4);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	outputFile.open("NWSTP_graph_" + save_name + ".stp"); // stp file

	outputFile << "33D32945 STP File, STP Format Version 1.0" << endl;
	outputFile << endl;

	// comments
	outputFile << "SECTION Comments" << endl;
	outputFile << "Name \"" << save_name << "\"" << endl;
	outputFile << "Creator \"Yahui Sun\"" << endl;
	outputFile << "Problem \"Node-weighted Steiner Problem in Graphs\"" << endl;
	outputFile << "END" << endl;
	outputFile << endl;

	// graph
	outputFile << "SECTION Graph" << endl;
	outputFile << "Nodes " << num_vertices(result_graph) << endl;
	outputFile << "Edges " << num_edges(result_graph) << endl;
	graph::out_edge_iterator eit, eend;
	for (int i = 0; i < num_vertices(result_graph); i++) {
		outputFile << "V " << new_2_initial[i] << " " << get(boost::vertex_name_t(), result_graph, i) << endl;
	}
	for (int i = 0; i < num_vertices(result_graph); i++) {
		tie(eit, eend) = boost::out_edges(i, result_graph); // adjacent_vertices of 2
		for_each(eit, eend,
			[&result_graph, &i, &outputFile, &new_2_initial](graph::edge_descriptor it)
		{
			int j = boost::target(it, result_graph);
			if (i < j) {
				outputFile << "E " << new_2_initial[i] << " " << new_2_initial[j] << " " << get(boost::edge_weight_t(), result_graph, boost::edge(i, j, result_graph).first) << endl;
			}
		});
	}
	outputFile << "END" << endl;
	outputFile << endl;

	//// TP
	//outputFile << "SECTION Non-Compulsory Terminals" << endl;
	//int p = 0;
	//for (int i = 0; i < num_vertices(result_graph); i++) {
	//	if (get(boost::vertex_name_t(), result_graph, i) > 0) {
	//		p++;
	//	}
	//}
	//outputFile << "Terminals " << p << endl;
	//for (int i = 0; i < num_vertices(result_graph); i++) {
	//	if (get(boost::vertex_name_t(), result_graph, i) > 0) {
	//		outputFile << "TP " << i + 1 << " " << get(boost::vertex_name_t(), result_graph, i) << endl;
	//	}
	//}
	//outputFile << "END" << endl;
	//outputFile << endl;

	//// CTP
	//outputFile << "SECTION Compulsory Terminals" << endl;
	//p = 0;
	//for (int i = 0; i < num_vertices(result_graph); i++) {
	//	if (terminal[i] == 1) {
	//		p++;
	//	}
	//}
	//outputFile << "CompulsoryTerminals " << p << endl;
	//for (int i = 0; i < num_vertices(result_graph); i++) {
	//	if (terminal[i] == 1) {
	//		outputFile << "CTP " << i + 1 << endl;
	//	}
	//}
	//outputFile << "END" << endl;
	//outputFile << endl;


	outputFile << "EOF" << endl;

}
#pragma endregion save_data


#pragma region

graph FGW(graph input_graph, double& growth_time, double distribution_ratio) {

	double Global_time = 0; // global time
	int Active_C_num = 0; // the number of active clusters
	int N = num_vertices(input_graph); // number of vertices
	int ep_num = num_edges(input_graph) * 2; // total number of edge parts
	int ep_order = 0;
	node node0;

	// Clusters: the number of clusters is always N
	std::vector<bool> C_activity(N); // activity value of each C; false means inactive; initial value is false
	std::vector<double> C_event_time(N); // the event time for each C
	std::vector<double> C_deactivate_time(N); // the deactivate time for each C
	std::vector<std::vector<int>> C_V_list(N); // record the vertices in each C
	std::vector<pairing_heap<node>> C_ep_PQ(N); // the PQ for edge parts in each C; node index: ep order in ep_list
	std::vector<int> V_locator(N); // the index of the maximal cluster containing the vertex
							  // edge parts: PQ and their handles
	std::vector<int> ep_v1_list(ep_num); // class may be slow, so I seperate the ep_list
	std::vector<int> ep_v2_list(ep_num);
	std::vector<double> ep_EventTime_list(ep_num);
	std::vector<int> ep_ep2_order_list(ep_num);
	std::vector<handle_t> handle_ep(ep_num); // store the handle for each edge part
										// the event PQ and their handles
	pairing_heap<node> C_event_PQ; // PQ storing the event time of the active clusters; node index: cluster order
	std::vector<handle_t> handle_Cevent(N);
	pairing_heap<node> E_event_PQ; // PQ storing the active clusters; node index: cluster order
	std::vector<handle_t> handle_Eevent(N);

	graph::out_edge_iterator eit, eend;

	// initialize the clusters
	for (int i = 0; i < N; i++)
	{
		C_V_list[i].insert(C_V_list[i].end(), i); // insert a vertex into the rear of C_V_list[i]
		V_locator[i] = i; // the maximal cluster containing vertex i
						  // add edge parts into C
		tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of i
		for_each(eit, eend,
			[&input_graph, &ep_v1_list, &ep_v2_list, &ep_ep2_order_list, &handle_ep, &C_ep_PQ, &node0, // the & above is the capture-list: the variables you can use inside
			&ep_order, &i, &ep_EventTime_list, &distribution_ratio](graph::edge_descriptor it) // for each adjacenct vertex boost::target(it, input_graph)
		{
			int j = boost::target(it, input_graph); // the adjacent vetex to i
			if (j > i) { // don't overcheck an edge
						 // the first ep
				ep_v1_list[ep_order] = i;
				ep_v2_list[ep_order] = j;
				ep_EventTime_list[ep_order] = get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first) / distribution_ratio; // halve the edge cost
				ep_ep2_order_list[ep_order] = ep_order + 1; // it points to the ep below
				node0.index = ep_order; // node index: ep order
				node0.priority_value = ep_EventTime_list[ep_order]; // priority: ep_EventTime
				handle_ep[ep_order] = C_ep_PQ[i].push(node0); // put this ep into cluster i
				ep_order++;
				// the second ep
				ep_v1_list[ep_order] = j;
				ep_v2_list[ep_order] = i;
				ep_EventTime_list[ep_order] = ep_EventTime_list[ep_order - 1] * (distribution_ratio - 1); // halve the edge cost
				ep_ep2_order_list[ep_order] = ep_order - 1; // it points to the ep above
				node0.index = ep_order; // node index: ep order
				node0.priority_value = ep_EventTime_list[ep_order]; // priority: ep_EventTime
				handle_ep[ep_order] = C_ep_PQ[j].push(node0); // put this ep into cluster j
				ep_order++;
			}
		});
		// for active cluster
		if (get(boost::vertex_name_t(), input_graph, i) > 0) {
			Active_C_num++; // the number of active clusters
			C_activity[i] = true; // this cluster is active
			C_event_time[i] = get(boost::vertex_name_t(), input_graph, i); // the event time is the node weight
																		   // push node into C_event_PQ
			node0.index = i; // node index: cluster order
			node0.priority_value = C_event_time[i]; // priority: node weight
			handle_Cevent[i] = C_event_PQ.push(node0); // into PQ
													   // all the ep for cluster i have been inserted into C_ep_PQ[i]; Note that, it's only true when i starts from 0 and j>i above
													   // push node into E_event_PQ
			node0.priority_value = C_ep_PQ[i].top().priority_value; // priority: the closest ep time
			handle_Eevent[i] = E_event_PQ.push(node0); // into PQ
		}
		else if (get(boost::vertex_name_t(), input_graph, i) < 0) {
			C_event_time[i] = get(boost::vertex_name_t(), input_graph, i); // the C_event_time is below 0 for negatively weighted nodes
		}
	}


	// FGW growth starts!
	graph output_graph(N); // the output graph
	for (int i = 0; i < N; i++) {
		double new_weight = get(boost::vertex_name_t(), input_graph, i);
		boost::put(boost::vertex_name_t(), output_graph, i, new_weight); // input node weights
	}
	int C0;
	int C1;
	int C2;
	int ep1;
	int ep2;
	double Tc;
	double Te;
	double r;
	double lowerbound = 1e-7;  // d is not used in this coding

	auto begin_time = std::chrono::high_resolution_clock::now(); // start time

	while (Active_C_num > 1) // stop the loop when there is only one active cluster left
	{
		// find the closest event
		Tc = C_event_PQ.top().priority_value; // this cluster event time
		Te = E_event_PQ.top().priority_value; // this edge event time

		if (Tc >= Te) { // edge event
			C1 = E_event_PQ.top().index; // the cluster C1 for this edge event
			ep1 = C_ep_PQ[C1].top().index; // the top ep in C1
			C2 = V_locator[ep_v2_list[ep1]]; // the cluster C2 for this edge event

			if (C1 == C2) { // the inside ep is triggered
				C_ep_PQ[C1].pop(); // pop out the inside ep
								   // decrease E_event_PQ for the change of event_C1
				node0.index = C1;
				node0.priority_value = C_ep_PQ[C1].top().priority_value; // theoretically C_ep_PQ[event_C1] should not be empty
				E_event_PQ.decrease(handle_Eevent[C1], node0);
			}
			else { // the outside ep is triggered
				Global_time = Te;
				ep2 = ep_ep2_order_list[ep1];

				if (C_activity[C2] == true) { // C2 is active
					r = ep_EventTime_list[ep2] - Global_time; // the slack of the responsible edge

					if (r > lowerbound) { // r is big; d is not used in this coding
										  // change two ep event time
						ep_EventTime_list[ep1] = Global_time + r / 2;
						ep_EventTime_list[ep2] = Global_time + r / 2;
						// update C_ep_PQ in C1
						node0.index = ep1;
						node0.priority_value = ep_EventTime_list[ep1];
						C_ep_PQ[C1].decrease(handle_ep[ep1], node0);
						// update C_ep_PQ in C2
						node0.index = ep2;
						node0.priority_value = ep_EventTime_list[ep2];
						C_ep_PQ[C2].increase(handle_ep[ep2], node0);
						// update E_event_PQ for the change of C1
						node0.index = C1;
						node0.priority_value = C_ep_PQ[C1].top().priority_value;
						E_event_PQ.decrease(handle_Eevent[C1], node0);
						// update E_event_PQ for the change of C2
						if (C_ep_PQ[C2].top().index == ep2) {
							node0.index = C2;
							node0.priority_value = C_ep_PQ[C2].top().priority_value;
							E_event_PQ.increase(handle_Eevent[C2], node0);
						}
					}
					else { // r is small; merge event
						   // add edge with the original cost
						boost::add_edge(ep_v1_list[ep1], ep_v1_list[ep2],
							get(boost::edge_weight_t(), input_graph,
								boost::edge(ep_v1_list[ep1], ep_v1_list[ep2], input_graph).first), output_graph);
						// merge V_list of C2 into C1
						C_V_list[C1].insert(end(C_V_list[C1]), begin(C_V_list[C2]), end(C_V_list[C2]));
						//decrease V_locator
						for (int i = 0; i < C_V_list[C2].size(); i++) {
							V_locator[C_V_list[C2][i]] = C1;
						}
						// update event time of C1
						C_event_time[C1] = C_event_time[C1] + C_event_time[C2] - Global_time;
						// minus one active cluster
						Active_C_num--;
						C_activity[C2] = false;
						C_deactivate_time[C2] = Global_time;
						// merge two C_ep_PQ
						C_ep_PQ[C1].pop(); // pop out the responsible ep
						C_ep_PQ[C1].merge(C_ep_PQ[C2]);
						// update C1 in C_event_time and E_event_time
						node0.index = C1;
						node0.priority_value = C_event_time[C1];
						C_event_PQ.decrease(handle_Cevent[C1], node0);
						node0.priority_value = C_ep_PQ[C1].top().priority_value;
						E_event_PQ.decrease(handle_Eevent[C1], node0);
						// remove C2 from C_event_time and E_event_time
						C_event_PQ.erase(handle_Cevent[C2]);
						E_event_PQ.erase(handle_Eevent[C2]);
					}
				}
				else { // C2 is inactive
					r = ep_EventTime_list[ep2] - C_deactivate_time[C2]; // the slack of the responsible edge

					if (r > lowerbound) { // r is big; d is not used in this coding
										  // change two ep event time
						ep_EventTime_list[ep1] = Global_time + r;
						ep_EventTime_list[ep2] = C_event_time[C2];
						// update C_ep_PQ in C1
						node0.index = ep1;
						node0.priority_value = ep_EventTime_list[ep1];
						C_ep_PQ[C1].decrease(handle_ep[ep1], node0);
						// update C_ep_PQ in C2
						node0.index = ep2;
						node0.priority_value = ep_EventTime_list[ep2];
						C_ep_PQ[C2].increase(handle_ep[ep2], node0);
						// update E_event_PQ for the change of C1
						node0.index = C1;
						node0.priority_value = C_ep_PQ[C1].top().priority_value;
						E_event_PQ.decrease(handle_Eevent[C1], node0);
					}
					else { // r is small; merge event
						   // add edge
						boost::add_edge(ep_v1_list[ep1], ep_v1_list[ep2],
							get(boost::edge_weight_t(), input_graph,
								boost::edge(ep_v1_list[ep1], ep_v1_list[ep2], input_graph).first), output_graph);
						// merge V_list of C2 into C1
						C_V_list[C1].insert(end(C_V_list[C1]), begin(C_V_list[C2]), end(C_V_list[C2]));
						//decrease V_locator
						for (int i = 0; i < C_V_list[C2].size(); i++) {
							V_locator[C_V_list[C2][i]] = C1;
						}
						// merge two C_ep_PQ
						C_ep_PQ[C1].pop(); // pop out the responsible ep		   
						typename pairing_heap<node>::iterator begin = C_ep_PQ[C2].begin();
						typename pairing_heap<node>::iterator end = C_ep_PQ[C2].end();
						for (typename pairing_heap<node>::iterator it = begin; it != end; ++it)
						{
							node0 = *it;
							if (V_locator[ep_v2_list[node0.index]] != C1) { // only push outside nodes into C_ep_PQ[event_C1]; it's a little faster than not do that
								node0.priority_value = node0.priority_value + Global_time - C_event_time[C2]; // decrease priority values
								handle_ep[node0.index] = C_ep_PQ[C1].push(node0); // push; decrease handle
							}
						}

						// update event time of C1
						C_event_time[C1] = C_event_time[C1] - Global_time + C_event_time[C2] - C_deactivate_time[C2] + Global_time;
						if (C_event_time[C1] > Global_time) { // new C1 is active
															  // update C1 in C_event_time and E_event_time
							node0.index = C1;
							node0.priority_value = C_event_time[C1];
							C_event_PQ.decrease(handle_Cevent[C1], node0);
							node0.priority_value = C_ep_PQ[C1].top().priority_value;
							E_event_PQ.decrease(handle_Eevent[C1], node0);
						}
						else { // new C1 is inactive; deactivate C1
							Active_C_num--; // minus one active cluster
							C_event_PQ.erase(handle_Cevent[C1]);
							E_event_PQ.erase(handle_Eevent[C1]);
							C_activity[C1] = false; // deactivate it
							C_deactivate_time[C1] = Global_time;
						}
					}
				}
			}
		}
		else { // cluster event
			Global_time = Tc; // decrease time
			C0 = C_event_PQ.top().index; // the cluster for this cluster event
			Active_C_num--; // minus one active cluster
			C_event_PQ.pop(); // remove the cluster from C_event_PQ
			E_event_PQ.erase(handle_Eevent[C0]); // remove the cluster from E_event_PQ
			C_activity[C0] = false; // deactivate it
			C_deactivate_time[C0] = Global_time;
		}
	}

	// remove disconnected parts
	std::vector<int> component(N); // vertex i is in component[i]; No.component from 0
	int cpn_num = connected_components(output_graph, &component[0]); // the number of component; decrease component
	int R_cpn;
	if (C_event_PQ.size() == 0) {
		R_cpn = component[C_V_list[C1][0]]; // C1 is the last active cluster
	}
	else {
		R_cpn = component[C_V_list[C_event_PQ.top().index][0]]; // it throw exception when TP=0 and C_event_PQ.size()=0
	}
	for (int i = 0; i < N; i++) {
		if (component[i] != R_cpn && in_degree(i, output_graph) > 0) { // disconnected vertex
			clear_vertex(i, output_graph); // clear_vertex removes adjacent vertices, but not node weight
		}
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count(); // Nanosecond
	growth_time = runningtime / 1e6;

	return output_graph;

}
#pragma endregion FGW_NWSTP


#pragma region
graph GPrA(graph input_graph, double& GPrA_time) {

	int N = num_vertices(input_graph); // number of vertices
	std::vector<double> nw1(N); // the nw value for finding R, it will be decreased after the check
	std::vector<double> nw2(N); // the nw value for pruning
	for (int i = 0; i < N; i++) {
		nw1[i] = get(boost::vertex_name_t(), input_graph, i);
		nw2[i] = get(boost::vertex_name_t(), input_graph, i);
	}
	std::vector<bool> pcheck1(N); // true means it has been checked 
	int num_check1 = 0; // the number of checked vertices to find R
	std::vector<bool> pcheck2(N); // true means it has been checked 
	int num_check2 = 0; // the number of checked vertices for pruning
	std::vector<int> pdegree1(N); // the processing degree to find R
	std::vector<int> pdegree2(N); // the processing degree for pruning
	std::vector<int> leaf1; // the leaves for finding R
	std::vector<int> leaf2; // the leaves for pruning

	for (int i = 0; i < N; i++) {
		pdegree1[i] = in_degree(i, input_graph); // decrease pdegree
		pdegree2[i] = pdegree1[i];
		if (pdegree1[i] == 0) {
			pcheck1[i] = true; // check disconnected vertices
			num_check1++;
			pcheck2[i] = true;
			num_check2++;
		}
		else if (pdegree1[i] == 1) {
			leaf1.insert(leaf1.end(), i);
			leaf2.insert(leaf2.end(), i);
		}
	}

	graph::out_edge_iterator eit, eend;
	AdjacencyIterator ai, a_end;
	int leaf_num1 = N - num_check1 - 1; // the number of leaves you need to process
	int leaf_num2 = N - num_check1 - 1; // the number of leaves you need to process

	auto begin_time = std::chrono::high_resolution_clock::now(); // start time

																 //this version is similar to the version below
	int k = 0;
	while (k < leaf_num1) {
		int i = leaf1[k];
		k++;
		tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of i
		for_each(eit, eend,
			[&input_graph, &pcheck1, &i, &nw1, &pdegree1, &num_check1, &leaf1](graph::edge_descriptor it)
		{
			int j = boost::target(it, input_graph);
			if (pcheck1[j] == false) {
				double cost = get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first);
				if (cost < nw1[i]) {
					nw1[j] = nw1[j] + nw1[i] - cost; // decrease nw1[j]
				}
				pcheck1[i] = true; // i has been checked, there is no need to delete it in this phase
				pdegree1[j]--;// decrease pdegree[j]
				if (pdegree1[j] == 1) {
					leaf1.insert(leaf1.end(), j); // it's fine to insert in the end, but not in the biginning
				}
				// break; // how to break a for_each???
			}
		});
	}

	//// find the Root, which is the mark of the optimal prunning result
	//while (num_check1 < N - 1) { 
	//	// there will be a single vertex left unchecked (but its nw value will be decreased)
	//	//// the version below is slower
	//	//while (leaf1.size() > 0) {
	//	//	int i = leaf1[0];
	//	//	leaf1.erase(leaf1.begin());
	//	//	tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of i
	//	//	for_each(eit, eend,
	//	//		[&input_graph, &pcheck1, &i, &nw1, &pdegree1, &num_check1, &leaf1](graph::edge_descriptor it)
	//	//	{
	//	//		int j = boost::target(it, input_graph);
	//	//		if (pcheck1[j] == false) {
	//	//			double cost = get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first);
	//	//			if (cost < nw1[i]) {
	//	//				nw1[j] = nw1[j] + nw1[i] - cost; // decrease nw1[j]
	//	//			}
	//	//			pcheck1[i] = true; // i has been checked, there is no need to delete it in this phase
	//	//			num_check1++;
	//	//			pdegree1[j]--;// decrease pdegree[j]
	//	//			if (pdegree1[j] == 1) {
	//	//				leaf1.insert(leaf1.end(), j); // it's fine to insert in the end, but not in the biginning; note for (int k = 0; k < leaf1.size(); k++)
	//	//			}
	//	//			// break; // how to break a for_each???
	//	//		}
	//	//	});
	//	//}
	//	//// the version below is fast
	//	for (int k = 0; k < leaf1.size(); k++)
	//	{
	//		int i = leaf1[k];
	//		if (pdegree1[i] == 1) {
	//			tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of i
	//			for_each(eit, eend,
	//				[&input_graph, &pcheck1, &i, &nw1, &pdegree1, &num_check1, &leaf1, &k](graph::edge_descriptor it)
	//			{
	//				int j = boost::target(it, input_graph);
	//				if (pcheck1[j] == false) {
	//					double cost = get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first);
	//					if (cost < nw1[i]) {
	//						nw1[j] = nw1[j] + nw1[i] - cost; // decrease nw1[j]
	//					}
	//					pcheck1[i] = true; // i has been checked, there is no need to delete it in this phase
	//					num_check1++;
	//					pdegree1[i] = 0; // it's not the leaf any more
	//					pdegree1[j]--; // decrease pdegree[j]
	//					if (pdegree1[j] == 1) {
	//						leaf1.insert(leaf1.end(), j); // it's fine to insert in the end, but not in the biginning; note for (int k = 0; k < leaf1.size(); k++)
	//					}
	//					// break; // how to break a for_each???
	//				}
	//			});
	//		}
	//		//// the version below is slower than that above
	//		//boost::tie(ai, a_end) = boost::adjacent_vertices(i, input_graph);
	//		//for (; ai != a_end; ai++) {
	//		//	int j = *ai;
	//		//	if (pcheck1[j] == false) {
	//		//		double cost = get(boost::edge_weight_t(), input_graph, boost::edge(i, j, input_graph).first);
	//		//		if (cost < nw1[i]) {
	//		//			nw1[j] = nw1[j] + nw1[i] - cost; // decrease nw1[j]
	//		//		}
	//		//		pcheck1[i] = true; // i has been checked, there is no need to delete it in this phase
	//		//		num_check1++;
	//		//		pdegree1[i] = 0; // it's not the leaf any more
	//		//		pdegree1[j]--;// decrease pdegree[j]
	//		//if (pdegree1[j] == 1) {
	//		//	leaf1.insert(leaf1.end(), j);
	//		//}
	//		//		break; // how to break a for_each???
	//		//	}
	//		//}
	//	}
	//}

	// R is the vertex with the biggest nw
	int R = 0;
	double max = nw1[0];
	for (int i = 1; i < N; i++) {
		if (nw1[i] > max) {
			max = nw1[i];
			R = i;
		}
	}

	// Strong pruning tree
	graph output_graph = input_graph; // the output graph

									  //this version is similar to the version below
	k = 0;
	while (k < leaf_num2 + 1) { // since R is ignored, it must be leaf_num2+1
		int i = leaf2[k];
		k++;
		if (i != R) {
			tie(eit, eend) = boost::out_edges(i, output_graph); // adjacent_vertices of i
			for_each(eit, eend,
				[&output_graph, &pcheck2, &i, &nw2, &pdegree2, &num_check2, &leaf2, &k](graph::edge_descriptor it)
			{
				int j = boost::target(it, output_graph);
				if (pcheck2[j] == false) {
					double cost = get(boost::edge_weight_t(), output_graph, boost::edge(i, j, output_graph).first);
					if (cost < nw2[i]) {
						nw2[j] = nw2[j] + nw2[i] - cost; // decrease nw2[j]
					}
					else {
						boost::remove_edge(i, j, output_graph); // remove edge(i,j)	
					}
					pcheck2[i] = true; // i has been checked
					pdegree2[j]--;// decrease pdegree[j]
					if (pdegree2[j] == 1) {
						leaf2.insert(leaf2.end(), j);
					}
					// break; // how to break a for_each???
				}
			});
		}
	}

	//while (num_check2 < N - 1) {
	//	for (int k = 0; k < leaf2.size(); k++)
	//	{
	//		int i = leaf2[k];
	//		if (pdegree2[i] == 1 && i != R) {
	//			tie(eit, eend) = boost::out_edges(i, output_graph); // adjacent_vertices of i
	//			for_each(eit, eend,
	//				[&output_graph, &pcheck2, &i, &nw2, &pdegree2, &num_check2, &leaf2, &k](graph::edge_descriptor it)
	//			{
	//				int j = boost::target(it, output_graph);
	//				if (pcheck2[j] == false) {
	//					double cost = get(boost::edge_weight_t(), output_graph, boost::edge(i, j, output_graph).first);
	//					if (cost < nw2[i]) {
	//						nw2[j] = nw2[j] + nw2[i] - cost; // decrease nw2[j]
	//					}
	//					else {
	//						boost::remove_edge(i, j, output_graph); // remove edge(i,j)	
	//					}
	//					pcheck2[i] = true; // i has been checked
	//					num_check2++;
	//					pdegree2[i] = 0; // it's not the leaf any more
	//					pdegree2[j]--;// decrease pdegree[j]
	//					if (pdegree2[j] == 1) {
	//						leaf2.insert(leaf2.end(), j);
	//					}
	//					// break; // how to break a for_each???
	//				}
	//			});
	//			// the version below is slower than that above
	//			//boost::tie(ai, a_end) = boost::adjacent_vertices(i, input_graph);
	//			//for (; ai != a_end; ai++) {
	//			//	int j = *ai;
	//			//	if (pcheck2[j] == false) {
	//			//		double cost = get(boost::edge_weight_t(), output_graph, boost::edge(i, j, output_graph).first);
	//			//		if (cost < nw2[i]) {
	//			//			nw2[j] = nw2[j] + nw2[i] - cost; // decrease nw2[j]
	//			//		}
	//			//		else {
	//			//			boost::remove_edge(i, j, output_graph); // remove edge(i,j)
	//			//			//// check
	//			//			//std::cout << "output_graph net_cost: " << net_cost(output_graph) << endl; // this line causes errors; becuase edge_descriptor? why?
	//			//		}
	//			//		pcheck2[i] = true; // i has been checked
	//			//		num_check2++;
	//			//		pdegree2[i] = 0; // it's not the leaf any more
	//			//		pdegree2[j]--;// decrease pdegree[j]
	//			//if (pdegree2[j] == 1) {
	//			//	leaf2.insert(leaf2.end(), j);
	//			//}
	//			//		break; 
	//			//	}
	//			//}
	//		}
	//	}
	//}

	// deleted disconnected parts
	std::vector<int> component(N); // vertex i is in component[i]; No.component from 0
	int cpn_num = connected_components(output_graph, &component[0]); // the number of component; decrease component
	for (int i = 0; i < N; i++) {
		if (component[i] != component[R]) { // disconnected vertex
			clear_vertex(i, output_graph); // clear_vertex removes adjacent vertices, but not node weight
		}
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count(); // Nanosecond
	GPrA_time = runningtime / 1e6;

	return output_graph;

}
#pragma endregion  GPrA  2016年11月26日18:06


#pragma region
void PO_evaluation(graph GPrA_graph, std::vector<int> initial_2_new, std::vector<std::vector<int>> path, std::vector<int> new_2_initial,
	std::vector<int> path_source, std::vector<int> path_terminal) {

	std::vector<int> terminal;
	terminal.resize(new_2_initial.size());
	for (int i = 0; i < path_source.size(); i++) {
		terminal[initial_2_new[path_source[i]]] = 1;
	}
	for (int i = 0; i < path_terminal.size(); i++) {
		terminal[initial_2_new[path_terminal[i]]] = 1;
	}

	int PO_N = 0;
	std::vector<int> PO_2_new;
	std::vector<int> new_2_PO;
	new_2_PO.resize(num_vertices(GPrA_graph));
	for (int i = 0; i < num_vertices(GPrA_graph); i++) {
		if (in_degree(i, GPrA_graph) > 0) {
			PO_2_new.insert(PO_2_new.end(), i);
			new_2_PO[i] = PO_N;
			PO_N++;
		}
	}
	std::vector<int> PO_terminal;
	PO_terminal.resize(PO_N);
	for (int i = 0; i < PO_N; i++) {
		if (terminal[PO_2_new[i]] == 1) {
			PO_terminal[i] = 1;
		}
	}

	graph total_flux_edge(PO_N);
	graph PO_graph(PO_N);
	for (int i = 0; i < num_vertices(GPrA_graph); i++) {
		if (in_degree(i, GPrA_graph) > 0) {
			int v1 = new_2_PO[i];
			boost::put(boost::vertex_name_t(), PO_graph, v1, get(boost::vertex_name_t(), GPrA_graph, i));
			graph::out_edge_iterator eit, eend;
			tie(eit, eend) = boost::out_edges(i, GPrA_graph); // adjacent_vertices of 2
			for_each(eit, eend,
				[&GPrA_graph, &i, &PO_graph, &new_2_PO, &v1](graph::edge_descriptor it)
			{ int j = boost::target(it, GPrA_graph);
			if (j > i) {
				int v2 = new_2_PO[j];
				boost::add_edge(v1, v2, get(boost::edge_weight_t(), GPrA_graph, boost::edge(i, j, GPrA_graph).first), PO_graph);
			}
			});
		}
	}

	total_flux_edge = PO_node_weight_edge_cost(PO_graph, PO_terminal);
	std::vector<int> vertex_PO_centrality;
	vertex_PO_centrality.resize(initial_2_new.size());
	for (int i = 0; i < PO_N; i++) {
		double flux_vertex;
		graph::out_edge_iterator eit, eend;
		tie(eit, eend) = boost::out_edges(i, total_flux_edge); // adjacent_vertices of 2
		for_each(eit, eend,
			[&total_flux_edge, &i, &flux_vertex](graph::edge_descriptor it)
		{ int j = boost::target(it, total_flux_edge);
		flux_vertex = flux_vertex + get(boost::edge_weight_t(), total_flux_edge, boost::edge(i, j, total_flux_edge).first);
		});
		vertex_PO_centrality[new_2_initial[PO_2_new[i]]] = flux_vertex;
	}

	ofstream outputFile;
	outputFile.open("PO_centrality.csv"); // stp file
	outputFile << "vertex_PO_centrality" << endl;
	outputFile << "Protein,Centrality" << endl;
	int a = 1;
	while (a > 0) {
		int max_num = distance(vertex_PO_centrality.begin(), max_element(vertex_PO_centrality.begin(), vertex_PO_centrality.end()));
		int max_centrality = vertex_PO_centrality[max_num];
		if (max_centrality > 0) {
			outputFile << max_num << "," << max_centrality << endl;
			vertex_PO_centrality[max_num] = 0;
		}
		else {
			a = 0;
		}
	}
	outputFile << "edge_PO_centrality" << endl;
	outputFile << "Vertex 1,Vertex 2,Centrality" << endl;
	cout << num_edges(total_flux_edge) << endl;
	while (num_edges(total_flux_edge) > 0) {
		int max_centrality = -1;
		int max_v1;
		int max_v2;
		for (int i = 0; i < num_vertices(total_flux_edge); i++) {
			graph::out_edge_iterator eit, eend;
			tie(eit, eend) = boost::out_edges(i, total_flux_edge); // adjacent_vertices of i
			for_each(eit, eend,
				[&total_flux_edge, &i, &max_centrality, &max_v1, &max_v2](graph::edge_descriptor it)
			{ int v2 = boost::target(it, total_flux_edge);
			int centrality = get(boost::edge_weight_t(), total_flux_edge, boost::edge(i, v2, total_flux_edge).first);
			if (centrality >= max_centrality) {
				max_centrality = centrality;
				max_v1 = i;
				max_v2 = v2;
			}
			});
		}
		boost::remove_edge(max_v1, max_v2, total_flux_edge);
		outputFile << new_2_initial[PO_2_new[max_v1]] << "," << new_2_initial[PO_2_new[max_v2]] << "," << max_centrality << endl;
	}

}
#pragma endregion PO_evaluation


#pragma region
void evaluate_subnetworks_single_complex_path(graph GPrA_graph, std::vector<int> initial_2_new, std::vector<std::vector<int>> path, std::vector<int> new_2_initial,
	std::vector<int> path_source, std::vector<int> path_terminal) {

	int N = num_vertices(GPrA_graph);
	graph identify_graph(N);
	double SP_distance;

	ofstream outputFile;
	outputFile.open("Identified_edges_in_single_complex_path.csv"); // stp file

	for (int i = 0; i < path.size(); i++) {
		int v1 = path[i][0];
		int v2 = path[i][1];
		typedef graph::edge_descriptor Edge;
		pair<Edge, bool> ed = boost::edge(initial_2_new[v1], initial_2_new[v2], GPrA_graph);
		if (ed.second) {
			boost::add_edge(initial_2_new[v1], initial_2_new[v2], 1, identify_graph);
		}
	}
	outputFile.close();
	cout << "Subnetwork edges' overlap percentage: " << double(num_edges(identify_graph)) / num_edges(GPrA_graph) * 100 << "%" << endl;

	auto evaluate_begin_time = std::chrono::high_resolution_clock::now(); // start time

	std::vector<int> vertex_betweenness_centrality;
	vertex_betweenness_centrality.resize(initial_2_new.size());
	graph edge_betweenness_centrality(initial_2_new.size()); // edge cost is the centrality
	for (int i = 0; i < path_source.size(); i++) {
		for (int j = 0; j < path_terminal.size(); j++) {
			std::vector<int> pathSP = SP_Yahui(GPrA_graph, initial_2_new[path_source[i]], initial_2_new[path_terminal[j]], SP_distance);
			for (int k = 0; k < pathSP.size(); k++) {
				int v1 = new_2_initial[pathSP[k]];
				vertex_betweenness_centrality[v1] = vertex_betweenness_centrality[v1] + 1;
			}
			for (int k = 0; k < pathSP.size() - 1; k++) {
				int v1 = new_2_initial[pathSP[k]];
				int v2 = new_2_initial[pathSP[k + 1]];
				typedef graph::edge_descriptor Edge;
				pair<Edge, bool> ed = boost::edge(v1, v2, edge_betweenness_centrality);
				if (ed.second) {
					int new_weight = get(boost::edge_weight_t(), edge_betweenness_centrality, boost::edge(v1, v2, edge_betweenness_centrality).first) + 1;
					boost::put(boost::edge_weight_t(), edge_betweenness_centrality, ed.first, new_weight);
				}
				else { // cout << "This edge does not exist!" << endl;
					boost::add_edge(v1, v2, 1, edge_betweenness_centrality);
				}
			}
		}
	}

	auto evaluate_end_time = std::chrono::high_resolution_clock::now();
	double evaluate_runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(evaluate_end_time - evaluate_begin_time).count(); // Nanosecond
	cout << "evaluate time:" << evaluate_runningtime / 1e6 << " ms" << endl;


	outputFile.open("Subnetworks.csv"); // stp file
	graph::out_edge_iterator eit, eend;
	outputFile << "Protein 1,Protein 2,Edge B" << endl;
	for (int i = 0; i < N; i++) {
		tie(eit, eend) = boost::out_edges(i, GPrA_graph); // adjacent_vertices of 2
		for_each(eit, eend,
			[&GPrA_graph, &i, &outputFile, &new_2_initial, &vertex_betweenness_centrality, &edge_betweenness_centrality](graph::edge_descriptor it)
		{
			int j = boost::target(it, GPrA_graph);
			if (i < j) {
				outputFile << new_2_initial[i] << "," << new_2_initial[j] << ","
					<< vertex_betweenness_centrality[new_2_initial[i]] << "," << vertex_betweenness_centrality[new_2_initial[j]] << ","
					<< get(boost::edge_weight_t(), edge_betweenness_centrality, boost::edge(new_2_initial[i], new_2_initial[j], edge_betweenness_centrality).first) << endl;
			}
		});
	}
	outputFile << "Protein,Protein B" << endl;
	for (int i = 0; i < initial_2_new.size(); i++) {
		if (vertex_betweenness_centrality[i] > 0) {
			outputFile << i << "," << vertex_betweenness_centrality[i] << endl;
		}
	}
	outputFile.close();


	outputFile.open("Betweenness_centrality.csv"); // stp file
	outputFile << "vertex_betweenness_centrality" << endl;
	outputFile << "Protein,Centrality" << endl;
	int a = 1;
	while (a > 0) {
		int max_num = distance(vertex_betweenness_centrality.begin(), max_element(vertex_betweenness_centrality.begin(), vertex_betweenness_centrality.end()));
		int max_centrality = vertex_betweenness_centrality[max_num];
		if (max_centrality > 0) {
			outputFile << max_num << "," << max_centrality << endl;
			vertex_betweenness_centrality[max_num] = 0;
		}
		else {
			a = 0;
		}
	}
	outputFile << "edge_betweenness_centrality" << endl;
	outputFile << "Vertex 1,Vertex 2,Centrality" << endl;
	while (num_edges(edge_betweenness_centrality) > 0) {
		int max_centrality = 0;
		int max_v1;
		int max_v2;
		for (int i = 0; i < initial_2_new.size(); i++) {
			graph::out_edge_iterator eit, eend;
			tie(eit, eend) = boost::out_edges(i, edge_betweenness_centrality); // adjacent_vertices of i
			for_each(eit, eend,
				[&edge_betweenness_centrality, &i, &max_centrality, &max_v1, &max_v2](graph::edge_descriptor it)
			{ int v2 = boost::target(it, edge_betweenness_centrality);
			int centrality = get(boost::edge_weight_t(), edge_betweenness_centrality, boost::edge(i, v2, edge_betweenness_centrality).first);
			if (centrality > max_centrality) {
				max_centrality = centrality;
				max_v1 = i;
				max_v2 = v2;
			}
			});
		}
		boost::remove_edge(max_v1, max_v2, edge_betweenness_centrality);
		outputFile << max_v1 << "," << max_v2 << "," << max_centrality << endl;
	}


}
#pragma endregion  evaluate_subnetworks_single_complex_path


#pragma region
void save_subnetworks(graph input_graph, std::vector<int> new_2_initial) {

	ofstream outputFile;
	outputFile.open("Subnetworks.csv"); // stp file
	int N = num_vertices(input_graph);
	graph::out_edge_iterator eit, eend;
	for (int i = 0; i < N; i++) {
		tie(eit, eend) = boost::out_edges(i, input_graph); // adjacent_vertices of 2
		for_each(eit, eend,
			[&input_graph, &i, &outputFile, &new_2_initial](graph::edge_descriptor it)
		{
			int j = boost::target(it, input_graph);
			if (i < j) {
				outputFile << new_2_initial[i] << "," << new_2_initial[j] << endl;
			}
		});
	}
}
#pragma endregion save_subnetworks


#pragma region 
void save_connected_graph_for_Grephi(graph raw_graph, std::vector<int> new_2_initial) {
	ofstream outputFile;
	outputFile.open("Edges.csv"); // stp file

	int N = num_vertices(raw_graph);
	graph::out_edge_iterator eit, eend;
	for (int i = 0; i < N; i++) {
		tie(eit, eend) = boost::out_edges(i, raw_graph); // adjacent_vertices of 2
		for_each(eit, eend,
			[&raw_graph, &i, &outputFile, &new_2_initial](graph::edge_descriptor it)
		{
			int j = boost::target(it, raw_graph);
			if (i < j) {
				outputFile << "9606.ENSP" + to_string(new_2_initial[i]) << "," << "9606.ENSP" + to_string(new_2_initial[j]) << endl;
			}
		});
	}
	outputFile.close();

	outputFile.open("IDs.csv"); // stp file
	outputFile << "Id,Label" << endl;
	for (int i = 0; i < N; i++) {
		outputFile << i << "," << "9606.ENSP" + to_string(new_2_initial[i]) << endl;
	}
	outputFile.close();
}
#pragma endregion save_connected_graph_for_Grephi


#pragma region
void save_NWSTP_graph(string instance_name, graph result_graph) {

	string save_name = instance_name; // save_name
	ofstream outputFile;
	outputFile.precision(4);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	outputFile.open(save_name + ".stp"); // stp file

	outputFile << "33D32945 STP File, STP Format Version 1.0" << endl;
	outputFile << endl;

	// comments
	outputFile << "SECTION Comments" << endl;
	outputFile << "Name \"" << save_name << "\"" << endl;
	outputFile << "Creator \"Yahui Sun\"" << endl;
	outputFile << "Problem \"Node - Weighted Steiner Problem in Graphs\"" << endl;
	outputFile << "END" << endl;
	outputFile << endl;

	// graph
	outputFile << "SECTION Graph" << endl;
	outputFile << "Nodes " << num_vertices(result_graph) << endl;
	outputFile << "Edges " << num_edges(result_graph) << endl;
	graph::out_edge_iterator eit, eend;
	for (int i = 0; i < num_vertices(result_graph); i++) {
		tie(eit, eend) = boost::out_edges(i, result_graph); // adjacent_vertices of 2
		for_each(eit, eend,
			[&result_graph, &i, &outputFile](graph::edge_descriptor it)
		{
			int j = boost::target(it, result_graph);
			if (i < j) {
				outputFile << "E " << i + 1 << " " << j + 1 << " " << get(boost::edge_weight_t(), result_graph, boost::edge(i, j, result_graph).first) << endl;
			}
		});
	}
	outputFile << "END" << endl;
	outputFile << endl;

	// TP
	outputFile << "SECTION Node Weights" << endl;
	for (int i = 0; i < num_vertices(result_graph); i++) {
		outputFile << "TP " << i + 1 << " " << get(boost::vertex_name_t(), result_graph, i) << endl;
	}
	outputFile << "END" << endl;
	outputFile << endl;

	outputFile << "EOF" << endl;

}
#pragma endregion save_NWSTP_graph  2017年6月25日












int main()
{
	string file_name = "9606.protein.links.detailed.v10";

	graph raw_graph = read_PPI_data(file_name); // vertex number is directly the protein index (not starts from 0); edge cost is the combined score
	// cout << in_degree(0, raw_graph) << endl; // prove that vertex 0 in raw_graph is not a protein


	std::vector<std::vector<int>> path = read_path();  // a single complex path
	std::vector<int> path_source;
	std::vector<int> path_terminal;
	source_terminal(path_source, path_terminal);

	graph connected_graph;
	std::vector<int> initial_2_new;
	std::vector<int> new_2_initial;
	connected_graph = customizing_graph_single_complex_path(raw_graph, initial_2_new, new_2_initial, path, path_source, path_terminal);
	cout << num_vertices(connected_graph) << " " << num_edges(connected_graph) << endl;
	//save_data(file_name, connected_graph, new_2_initial); // save initial graph
	//save_connected_graph_for_Grephi(connected_graph, new_2_initial);
	//save_NWSTP_graph("PPI_Network", connected_graph);

	//double growth_time = 0;
	//double distribution_ratio = 2;
	//graph Growth_graph = FGW(connected_graph, growth_time, distribution_ratio);
	//double GPrA_time = 0;
	//graph GPrA_graph = GPrA(Growth_graph, GPrA_time);
	////// check results
	////cout << "Number of Edges: " << num_edges(GPrA_graph) << endl;
	////int num_v = 0;
	////for (int i = 0; i < num_vertices(GPrA_graph); i++) {
	////	if (in_degree(i, GPrA_graph) > 0) {
	////		num_v++;
	////	}
	////}
	////cout << "Number of Vertices: " << num_v << endl;
	//cout << "GW time: " << growth_time + GPrA_time << "ms" << endl;

	//evaluate_subnetworks_single_complex_path(GPrA_graph, initial_2_new, path, new_2_initial, path_source, path_terminal);
	////PO_evaluation(GPrA_graph, initial_2_new, path, new_2_initial, path_source, path_terminal);







	cout << "END" << endl;

	getchar();
}