#include <iostream>
#include <stdlib.h>
#include <climits> // INT_MAX
using namespace std;

#include <fstream>
#include <string>
#include <string.h>
#include <sstream>  
#include <map>
#include <vector>
#include <time.h>
#include <algorithm>
#include <cmath>

#define print_flag 0
#define read_flag 0
#define topology_flag 0
#define label_flag 0
#define map_flag 0
#define output_flag 0
#define decomposition_flag 0



string model_name;
vector<string> G_node, G_input, G_output, G_intermediate, G_always_0;
vector<string> topology_order;

map<string, int> topology_label;
map<string, int> k_feasible_label;
map<string, int>::iterator it;

map<string, vector<string> > Nt;
map<string, vector<string> > Nt_PI;


map<string, vector<string> > X, X_bar, inputs_of_X_bar;
map<string, int> cut_size;
map<string, int> record;

int K;

class Edge{

public:
	friend class Gate;
	Edge(string a, string b) : in(a), out(b){}

	string in;
	string out;

};
vector<Edge> G_edge;

class Gate{
	
public:
	friend class Edge;
	Gate(string n,int size) : fanout(n), fanin_size(size){
		this->is_two_fanin = 0;
		this->is_mapped = 0;
	}

	void set_type(int type){
		this->gate_type = type;
	}

	void check_label(void){

		bool flag = 1;
		int max_label = 0;
		for(int i = 0 ; i < fanin_size ; i++){
			
			string target = this->fanin[i];
			//cout << "target : " << target << " ";
			it = topology_label.find(target);

			if( it == topology_label.end()){
				// not found 
				// fail
				//cout<< "not found" << endl;
				flag = 0;
				i = fanin_size - 1;
			}
			else{
				// found 
				// compare with max label
				//cout<< "found "<< label[target] << endl;
				if(topology_label[target] > max_label){
					max_label = topology_label[target];
				}
			}
			
		}

		// if all fanin are labeled, label the fanout
		if( flag == 1){
			// write label of the fanout
			//cout<< "label fanout with " << max_label + 1 << endl;
			topology_label[this->fanout] = max_label + 1;
		}
	}

	void create_edges(void){

		for(int i = 0 ; i < this->fanin_size ; i++){
			Edge e(this->fanin[i] , this->fanout);
			this->edge.push_back(e);
			G_edge.push_back(e);
		}
	}

	bool generate_output(bool in1, bool in2){
		
		bool out;
		if(this -> gate_type == 1){
			out = !in1;
		}
		else if(this ->gate_type == 2){
			out = in1 | in2;
		}
		else if(this ->gate_type == 3){
			out = in1 & in2;
		}
		else if(this ->gate_type == 4){
			out = false;
		}
		else if(this ->gate_type == 5){
			out = true;
		}
		else if(this ->gate_type == 6){
			out = in1;
		}
		return out;
	}
	
	string fanout;

	// 1 for inverter, 
	// 2 for or-gate, 
	// 3 for and-gate,
	// 4 for always 0,
	// 5 for always 1,
	// 6 for buffer
	int gate_type;
	int fanin_size;
	bool is_two_fanin;
	bool is_mapped;
	
	vector<string> fanin;
	vector<Edge> edge;
};
vector<Gate> gate;
map<string, int> gate_index;

class LUT{

public :
	friend class Gate;
	LUT(int x, string s) : k(x), output(s) {
		this->level_count = 0;
	}

	void node_topology_sort(void){

		this->node_topology_label.clear();
	
		// label the input
		for(int i = 0 ; i < this->input.size(); i++){
			string s = this->input[i];
			this->node_topology_label[s] = 0;
			//cout<< "input " << s << "'s label : "<<this->node_topology_label[s] << endl;
		}

		// label all nodes
		//cout<< this->input.size() << endl;
		//cout<< this->node.size() << endl;
		

		while( this->node_topology_label.size() < ( this->node.size() + this->input.size() )  ){
			
			// test all nodes in LUT
			for(int i = 0 ; i < this->node.size() ; i++){

				//cout<< this->node_topology_label.size() << endl;
				string current_node = this->node[i];
				//cout << "current node "<< current_node;
				// check current node is labeled or not
				it = this->node_topology_label.find(current_node);

				// has not been labeled
				// check whether it can be labeled
				if(it == this->node_topology_label.end()){
					
					//cout<< " has not been labeled" << endl;

					bool flag = 1;
					int max_label = 0;
					// test all the fanin of current gate
					for(int j = 0 ; j < gate[gate_index[current_node]].fanin_size ; j++){

						string target = gate[gate_index[current_node]].fanin[j];
						//cout<< "fanin of current node "<< target;
						it = this->node_topology_label.find(target);

						if( it == this->node_topology_label.end()){
							// not found 
							// fail
							//cout<< " not found" << endl;
							flag = 0;

							// end the for loop
							j = gate[gate_index[current_node]].fanin_size - 1;
						}
						else{
							// found 
							// compare with max label
							//cout<< " found " << endl;
							if(this->node_topology_label[target] > max_label){
								max_label = this->node_topology_label[target];
							}
						}

					}

					// if all of fanin is labeled, label the current node
					if(flag){
						//cout<< "label " << current_node << " with "<<max_label + 1 << endl;
						this->node_topology_label[current_node] = max_label + 1;
					}
				}
				else{
					//cout<< " has been labeled" << endl;
				}
			}
		
		}
		
		
		
		// find the max label of topology sort
		int max = 0;
		for(int i = 0 ; i < this->node.size() ; i++){
			string current_node = this->node[i];
			//cout<< "current node : " << current_node << endl;
			if(this->node_topology_label[current_node] > max){
				max = this->node_topology_label[current_node];
			}
		}
		//cout<< "max : " << max << endl;

		// store the sorting order from label = 1
		for(int i = 1 ; i <= max ; i++){
			for(int j = 0 ; j < this->node.size() ; j++){
				string current_node = this->node[j];

				// if current node's label = current label
				if(this->node_topology_label[current_node] == i){
					this->node_topology_order.push_back(current_node);
				}
			}
		}
		
		//cout<< "topology sort order : ";
		//for(int i = 0 ; i < this->node_topology_order.size() ; i++){
		//	cout<< this->node_topology_order[i] << " ";
		//}
		//cout<< endl;

	}

	bool test_pattern(bool pattern[]){

		// put pattern into input node
		for(int i = 0 ; i < this->input.size() ; i++){
			node_signal[input[i]] = pattern[i];
		}

		//for(int i = this->input.size() - 1 ; i >= 0 ; i--){
		//	cout << node_signal[input[i]];
		//}
		
		
		for(int i = 0 ; i < node_topology_order.size() ; i++){
			
			string target_node = node_topology_order[i];
			//cout<< "current target node : " << node_topology_order[i] << endl;

			if(gate[gate_index[target_node]].gate_type == 2 || gate[gate_index[target_node]].gate_type == 3){
				
				bool in0 = node_signal[gate[gate_index[target_node]].fanin[0]];
				bool in1 = node_signal[gate[gate_index[target_node]].fanin[1]];
				node_signal[gate[gate_index[target_node]].fanout] = gate[gate_index[target_node]].generate_output(in0, in1);
			}
			else{
				bool in0 = node_signal[gate[gate_index[target_node]].fanin[0]];
				node_signal[gate[gate_index[target_node]].fanout] = gate[gate_index[target_node]].generate_output(in0, 0);
			}

		}
		
		//cout<< " " << node_signal[this->output] << endl;
		return node_signal[this->output];
	}
	
	
	int k;
	int level_count;

	string output;
	
	// nodes in X_bar
	vector<string> node;
	vector<string> input;

	// topology sort
	vector<string> node_topology_order;
	map<string, int> node_topology_label;

	// gate's connection
	map<string , bool > node_signal;
	

};
vector<LUT> K_LUT;
map<string, int> LUT_index;



void ReadFile(char filename[]);
void Decomposition(void);
void trace_back(string n);
void create_new_gate(int i);
void Decompose_Output(char filename[]);


void Show_info(void);

void calculate_gate_index(void);
int is_PI(string s);
int is_PO(string s);

// labeling_phase
void labeling_phase(void);
void Topology_sort(void);
void group_Nt(string t,int index);
void find_cut(string t,int index);
int calculate_p(string t,int index);
int calculate_cut_size(string t);

// mapping phase
void mapping_phase(void);
void generate_LUT(string v);

void calculate_level(void);
void propagate(string s,int level);

void Output(char filename[]);


int new_node_count;

// cut info
int node_cut_size;
int edge_cut_size;
int volume;
int height;

int max_level;
int PI_level;

int main(int argc, char *argv[]){
	
	// set k
	char *p;
	K = strtol (argv[2],&p,10);
	//cout<< "K = " << K << endl;

	char file[] = "2_input.blif";

	ReadFile(argv[3]);
	Decomposition();

	//Decompose_Output(file);
	
	//Show_info();
	calculate_gate_index();
	
	Topology_sort();
	
	labeling_phase();

	
	mapping_phase();
	
	if(G_always_0.size() > 0){
		for(int i = 0 ; i < G_always_0.size() ; i++){
			
			for(int j = 0 ; j < G_input.size() ; j++){
				if( G_input[j] == G_always_0[i]){
					G_input.erase(G_input.begin() + j);
					j = G_input.size() - 1;
				}
			}
		}

	}
	
	
	calculate_level();
	Output(argv[4]);
	
	cout<< "The circuit level is "<< max_level << "."<<endl;
	cout<< "The number of LUT is "<< K_LUT.size() << "."<<endl;
	return 0;
}

void ReadFile(char filename[]){
	
	//clock_t start, finish;
	//double duration;

	//cout<< "Start reading blif file" << endl;
	//start = clock();
	if(read_flag){
		cout<< "Open blif file : " << filename <<endl;
	}
	

	ifstream fin;
	fin.open(filename,ios::in);
	if( !fin ){
		cout<<"[Error] Fail to open the file"<<endl;
		return;
	}

	string temp, port, node;
	// model name
	fin >> temp >> model_name;
	if(read_flag){
		cout<< "Model name : "<< model_name << endl;
	}
	

	// input port
	fin >> temp;
	while(1){
		fin >> port;
		if( strcmp(port.c_str(),".outputs") == 0){
			break;
		}
		else{
			if( port[0] != '\\'){
				// create new node
				// store into input vector
				G_node.push_back(port);
				G_input.push_back(port);
			}
		}	
	}

	// output port
	while(1){
		fin >> port;
		if( strcmp(port.c_str(),".names") == 0){
			break;
		}
		else{
			if( port[0] != '\\'){
				// create new node
				// store into output vector
				G_node.push_back(port);
				G_output.push_back(port);
			}
		}
	}

	int node_count = 0;
	string node_buf[20];

	// truth table
	while(1){
		fin >> node;
		// stop at ".end"
		if( strcmp(node.c_str(),".end") == 0){
			break;
		}
		else{
			//cout<< node << endl;
			// truth table part
			if( node[0] == '0' || node[0] == '1' || node[0] == '-'){
				//cout<< "truth table part" << endl;
				//cout<< "node count : " << node_count << endl;

				if( node_count != 0){
					
					// show nodes of the gate
					//for(int i = 0 ; i < node_count ; i++){		
					//	cout<< node_buf[i] << " ";
					//}
					//cout<<endl;

					// check the last node is the output or not

					int flag = 0;
					for(int i = 0 ; i < G_output.size() ; i++){

						if( node_buf[node_count-1] == G_output[i] ){
							flag = 1;
							i = G_output.size() - 1;
						}
					}

					// not found
					if( flag == 0){
						
						// create the new node(last node)
						if(read_flag){
							cout<< "create intermediate node : " << node_buf[node_count-1] << endl;
						}
						
						// store into intermediate node vector
						G_node.push_back(node_buf[node_count-1]);
						G_intermediate.push_back(node_buf[node_count-1]);
						
					}
					else {
						if(read_flag){
							cout<< "node " << node_buf[node_count-1] << " is a output" << endl;
						}
						
					}
					

					// create a gate
					// fanout = the last node of buffer
					Gate my_gate(node_buf[node_count-1], node_count - 1);

					// invertor : 1
					// 0 1
					// always 1 : 5
					// - 1
					// buffer : 6
					// 1 1
					if (node_count == 2){
						
						if(node[0] == '0'){
							if(read_flag){
								cout<< "This gate is a not gate" << endl;
							}
							my_gate.set_type(1);
						}
						else if(node[0] == '-'){
							if(read_flag){
								cout<< "This gate is always 1" << endl;
							}
							
							my_gate.set_type(5);
						}
						else if(node[0] == '1'){
							if(read_flag){
								cout<< "This gate is a buffer" << endl;
							}

							my_gate.set_type(6);
						}
						
					}

					// or gate : 2
					else if (node[1] == '-'){
						if(read_flag){
							cout<< "This gate is an or gate" << endl;
						}
						
						my_gate.set_type(2);
						
					}
					// and gate : 3
					else{
						if(read_flag){
							cout<< "This gate is an and gate" <<endl;
						}
						
						my_gate.set_type(3);
						
					}

					// create edge between node_buf[i] and the last node(node_buf[node_count-1])
					// store fanin of the gate
					for(int i = 0 ; i < node_count - 1 ; i++){
						my_gate.fanin.push_back(node_buf[i]);
					}
					//cout<<endl;

					// push my_gate into vector gate
					my_gate.create_edges();
					gate.push_back(my_gate);

					node_count = 0;
				}
			}
			
			// nodes part
			else if( node[0] != '.' ){
				
				//cout<<"node part" << endl;
				//cout<< node <<endl;
				node_buf[node_count] = node;
				node_count++;
			}
			// always 0
			else if( node[0] == '.' && node_count == 1){
				
				if(read_flag){
					cout<< "This gate is always 0" << endl;
				}
				// create a new node
				
				G_input.push_back(node_buf[0]);
				G_always_0.push_back(node_buf[0]);
				
				// always 0 : 4
				Gate my_gate(node_buf[0], 0);
				my_gate.set_type(4);
				gate.push_back(my_gate);

				node_count = 0;
			}
		}
	}
	fin.close();

	//cout<< "Read blif file finished" << endl;
	//finish = clock();
	//cout<< "Time duration : " << (double)(finish - start) / CLOCKS_PER_SEC << " sec" << endl;
	//cout<< "====================================================================="  << endl;
}


void Decomposition(void){

	//clock_t start, finish;

	//cout<< "Start decomposing" << endl;
	//start = clock();

	// start from outputs
	for(int i = 0 ; i < G_output.size() ; i++){
		
		// output node : G_output[i]
		trace_back(G_output[i]);
	}
	//cout<< "Decomposition finished" << endl;
	//finish = clock();
	//cout<< "Time duration : " << (double)(finish - start) / CLOCKS_PER_SEC << " sec" << endl;
	//cout<< "====================================================================="  << endl;
}


void trace_back(string n){

	// search the gate which fanout is G_output[i]
	int i = 0;
	while(gate[i].fanout.compare(n) != 0) {
		i++;
		if( i == gate.size()){
			if(decomposition_flag){
				cout<< "reach input " << n << endl;
			}
			return;
		}
	}

	if( gate[i].is_two_fanin){
		if(decomposition_flag){
			cout<< n << " and it's fanin gates are two fanin gates" << endl;
		}
		return;
	}
	else if(gate[i].gate_type == 1){	
		if(decomposition_flag){
			cout<< n << " is a not gate" << endl;
		}
		
		trace_back(gate[i].fanin[0]);
				
	}
	else if(gate[i].gate_type == 4 ){
		if(decomposition_flag){
			cout<< n << " is always 0" << endl;
		}
		
		return;
	}
	else if(gate[i].gate_type == 5 ){
		if(decomposition_flag){
			cout<< n << " is always 1" << endl;
		}
		
		return;
	}
	else if(gate[i].gate_type == 6){
		if(decomposition_flag){
			cout<< n <<" is a buffer" << endl;
		}
		
		trace_back(gate[i].fanin[0]);
	}
	else{
		// if the number of fanin is 1 or 2, skip decomposition
		// Otherwise, decompose the gate

		if( gate[i].fanin_size > 2 ){
			
			if(decomposition_flag){
				cout<< "Decompose the gate " << n << endl;
			}

			new_node_count = 1;
			while( gate[i].fanin_size > 2 ){
				create_new_gate(i);
			}
			
			// set the gate is two a two fanin gate
			gate[i].is_two_fanin = 1;

			trace_back(gate[i].fanin[0]);
			trace_back(gate[i].fanin[1]);
		}
		else{

			if(decomposition_flag){
				cout<< n << " skip decomposition" << endl;
			}
			
			// set the gate is two a two fanin gate
			gate[i].is_two_fanin = 1;
			
			trace_back(gate[i].fanin[0]);
			trace_back(gate[i].fanin[1]);
				
		}
	}
}


void create_new_gate(int i){
	// stituation 1 : fanin = 3
	// stituation 2 : fanin > 3
	
	if( gate[i].fanin_size == 3){
		// create a node as the fanout of new gate
		stringstream ss;
		ss << new_node_count;
		string num;
		ss >> num;
		string name = gate[i].fanout+ "-" + num;
		G_intermediate.push_back(name);
		G_node.push_back(name);

		// create a new gate, select 1st and 2nd fanin of g as fanins of new gate
		Gate my_gate(name , 2);
		my_gate.set_type(gate[i].gate_type);
		my_gate.fanin.push_back(gate[i].fanin[0]);
		my_gate.fanin.push_back(gate[i].fanin[1]);
		gate.push_back(my_gate);

		// delete 1st and 2nd fanin 
		// push the new node into g's fanin
		gate[i].fanin.erase(gate[i].fanin.begin(), gate[i].fanin.begin() + 2);
		gate[i].fanin.push_back(name);

		new_node_count++;
		gate[i].fanin_size--;
	}
	else{
		// create two nodes as the fanout of two new gate
		stringstream ss1;
		ss1 << new_node_count;
		string num1;
		ss1 >> num1;

		stringstream ss2;
		ss2 << new_node_count+1;
		string num2;
		ss2 >> num2;

		string name1 = gate[i].fanout + "-"  + num1;
		string name2 = gate[i].fanout + "-"  + num2;
		G_intermediate.push_back(name1);
		G_intermediate.push_back(name2);
		G_node.push_back(name1);
		G_node.push_back(name2);

		// create a new gate(name1), , select 1st and 2nd fanin of g as fanins of new gate
		Gate my_gate1(name1, 2);
		my_gate1.set_type(gate[i].gate_type);
		my_gate1.fanin.push_back(gate[i].fanin[0]);
		my_gate1.fanin.push_back(gate[i].fanin[1]);
		gate.push_back(my_gate1);

		// create a new gate(name2), , select 3rd and 4th fanin of g as fanins of new gate
		Gate my_gate2(name2, 2);
		my_gate2.set_type(gate[i].gate_type);
		my_gate2.fanin.push_back(gate[i].fanin[2]);
		my_gate2.fanin.push_back(gate[i].fanin[3]);
		gate.push_back(my_gate2);

		// delete 1st ~ 4th fanin from g
		// push two new nodes into g's fanin
		gate[i].fanin.erase(gate[i].fanin.begin(), gate[i].fanin.begin() + 4);
		gate[i].fanin.push_back(name1);
		gate[i].fanin.push_back(name2);

		new_node_count = new_node_count + 2;
		gate[i].fanin_size = gate[i].fanin_size - 2;
	}
	
}


void Decompose_Output(char filename[]){

	//<<"Output blif file : "<< filename <<endl;

	ofstream fout;
	fout.open(filename,ios::out);

	// model name
	fout << ".model " << model_name << endl;

	// inputs
	fout << ".inputs ";
	for(int i = 0 ; i < G_input.size() ; i++){
		fout << G_input[i] << " ";
	}
	fout << endl;

	// outputs
	fout << ".outputs ";
	for(int i = 0 ; i < G_output.size() ; i++){
		fout << G_output[i] << " ";
	}
	fout << endl;



	// gates
	for(int i = 0 ; i < gate.size() ; i++){
		fout << ".names ";

		if( gate[i].fanin.size() > 2){
			cout<< "[Error] The gate has more fanin than 2" << endl;
			return;
		}

		if( gate[i].gate_type == 4){
			// always 0
			fout << gate[i].fanout << endl;
		}
		else{

			for(int j = 0 ; j < gate[i].fanin.size() ; j++){
				fout << gate[i].fanin[j] << " ";
			}
			fout << gate[i].fanout <<endl;
			
			// inverter
			if( gate[i].gate_type == 1){
				fout << 0 << " " << 1 << endl;
			}
			// or gate
			else if( gate[i].gate_type == 2){
				fout << 1 << "- " << 1 << endl;
				fout << "-" << 1 << " " << 1 << endl;
				
			} 
			// and gate
			else if( gate[i].gate_type == 3){
				fout << 11 << " " <<1 << endl;
			}
			// always 1
			else if( gate[i].gate_type == 5){
				fout << "- " << 1 << endl;
			}
			// buffer 
			else if( gate[i].gate_type == 6){
				fout << 1 << " " << 1 << endl;
			}

		}
		
	}
	fout << ".end";
	//cout<< "====================================================================="  << endl;
}


void Show_info(void){

	cout<< "Number of inputs : " << G_input.size() << endl;
	if(print_flag){
		for(int i = 0 ; i < G_input.size() ; i++){
			cout<< G_input[i] << " ";
		}
		cout<< endl;
	}
	

	cout<< "Number of outputs : " << G_output.size() << endl;
	if(print_flag){
		for(int i = 0 ; i < G_output.size() ; i++){
			cout<< G_output[i] << " ";
		}
		cout<< endl;
	}

	cout<< "Number of intermediate nodes : " << G_intermediate.size() << endl;
	if(print_flag){
		for(int i = 0 ; i < G_intermediate.size() ; i++){
			cout<< G_intermediate[i] << " ";
		}
		cout<< endl;
	}

	cout<< "Number of gates : " << gate.size() << endl;
	if(print_flag){
		for(int i = 0 ; i < gate.size() ; i++){
			cout<< "gate's fanout : " << gate[i].fanout << endl;
			cout<< "edges : "<<endl;
			for(int j = 0 ; j < gate[i].edge.size() ; j++){
				cout<< gate[i].edge[j].in << " " << gate[i].edge[j].out << endl;
			} 
		}
	}
	cout<< "====================================================================="  << endl;
}

void calculate_gate_index(void){

	// intermediate node
	for(int i = 0 ; i < G_intermediate.size(); i++){
		string s = G_intermediate[i];
		int index = 0;
		while(gate[index].fanout.compare(s) != 0) {
			index++;
			if( index == gate.size() ){
				cout<< "[Error] Can't find the gate which fanout is "<< s << endl;
				return;
			}
		}
		gate_index[s] = index;
	}



	// output node
	for(int i = 0 ; i < G_output.size(); i++){
		string s = G_output[i];
		int index = 0;
		while(gate[index].fanout.compare(s) != 0) {
			index++;
			if( index == gate.size() ){
				cout<< "[Error] Can't find the gate which fanout is "<< s << endl;
				return;
			}
		}
		gate_index[s] = index;
	}


}

int is_PI(string s){

	// is PI : flag = 1
	 // otherwise, flag = 0
	int flag = 0;
	for(int k = 0 ; k < G_input.size() ; k++){
		// found
		if(s == G_input[k]){
			flag = 1;
			k = G_input.size() - 1;
		}
	}

	return flag;
}

int is_PO(string s){

	// is PO : flag = 1
	 // otherwise, flag = 0
	int flag = 0;
	for(int k = 0 ; k < G_output.size() ; k++){
		// found
		if(s == G_output[k]){
			flag = 1;
			k = G_output.size() - 1;
		}
	}

	return flag;
}



void Topology_sort(void){

	//clock_t start, finish;
	//double duration;

	//cout<< "Start topology sorting"<<endl;
	//start = clock();

	// label the PI
	for(int i = 0 ; i < G_input.size() ; i++){
		string s = G_input[i];
		topology_label[s] = 0;
	}

	// label all nodes
	while( topology_label.size() < G_node.size()){

		for(int i = 0 ; i < gate.size() ; i++){
			//cout << "current fanout " << gate[i].fanout << endl;
			// if all of gate's input is labeled, label the fanout
			it = topology_label.find(gate[i].fanout);

			// the gate's fanout has not been labeled
			// check
			if( it == topology_label.end()){
				gate[i].check_label();
			}
			
		}
	}

	// find max label
	int max_label = 0;
	for(int i = 0 ; i < G_node.size() ; i++){
		if(topology_label[G_node[i]] > max_label){
			max_label = topology_label[G_node[i]];
		}
	}

	// don't have to put PI into topology sort order
	// start from label = 1
	for(int i = 1 ; i <= max_label ; i++){
		if(topology_flag){
			cout<< "label "<< i << " : ";
		}
		for( int j = 0 ; j < G_node.size() ; j++){
			if(topology_label[G_node[j]] == i ){
				if(topology_flag){
					cout<< G_node[j] << " ";
				}
				topology_order.push_back(G_node[j]);
			}
		}
		if(topology_flag){
			cout<< endl;
		}
		
	}

	if(topology_flag){
		cout<< "topology sort order : ";
		for(int i = 0 ; i < topology_order.size() ; i++){
			cout<< topology_order[i] << " ";
		}
		cout<< endl;
	}
	
	//cout<< "Topology sort finished"<<endl;
	//finish = clock();
	//cout<< "Time duration : " << (double)(finish - start) / CLOCKS_PER_SEC << " sec" << endl;
	//cout<< "====================================================================="  << endl;
}

void labeling_phase(void){

	//clock_t start, finish;
	//double duration;

	//cout<< "Start labeling phase"<<endl;
	//start = clock();
	
	
	// create source for all the PI
	G_node.push_back("Source");

	// set all PI's label as 0
	for(int i = 0 ; i < G_input.size() ; i++){
		k_feasible_label[G_input[i]] = 0;
	}

	for(int i = 0 ; i < topology_order.size() ; i++){

		string t = topology_order[i];

		int index = gate_index[t];
		group_Nt(t,index);

		// delete duplicate elements in Nt and Nt_PI
		sort(Nt[t].begin(), Nt[t].end());
		sort(Nt_PI[t].begin(), Nt_PI[t].end());
		Nt[t].erase(unique(Nt[t].begin(),Nt[t].end()),Nt[t].end());
		Nt_PI[t].erase(unique(Nt_PI[t].begin(),Nt_PI[t].end()),Nt_PI[t].end());
		
		find_cut(t,index);
	}

	//cout<< "Labeling phase finished"<<endl;
	//finish = clock();
	//cout<< "Time duration : " << (double)(finish - start) / CLOCKS_PER_SEC << " sec" << endl;
	//cout<< "====================================================================="  << endl;
}

void group_Nt(string t, int index){
	
	// push t into Nt[t]
	Nt[t].push_back(t);
	//cout << t << endl;
	// go to previous nodes 
	for(int j = 0 ; j < gate[index].fanin_size ; j++){
		
		// fanin_node = the fanin of the gate[i]
		string fanin_node = gate[index].fanin[j];
		//cout<< "fanin node : " << fanin_node << endl; 

		// check fanin_node is a PI or not
		int flag = is_PI(fanin_node);

		// t is a PI
		if(flag == 1){
			//cout<<"t is a PI"<<endl;
			Nt_PI[t].push_back(fanin_node);
		}
		// t is not a PI
		else{
			//cout<<"t is not a PI"<<endl;
			Nt[t].insert(Nt[t].end(), Nt[fanin_node].begin(), Nt[fanin_node].end());
			Nt_PI[t].insert(Nt_PI[t].end(), Nt_PI[fanin_node].begin(), Nt_PI[fanin_node].end());
		}
		
	}
}

void find_cut(string t, int index){

	// Initialize with 
	// X : PI in the Nt
	// X_bar : rest of nodes
	X[t].clear();
	X_bar[t].clear();

	// find the max label in Nt[t]
	int max_label = 1;
	for(int i = 0 ; i < Nt[t].size() ; i++){
		if(k_feasible_label[Nt[t][i]] > max_label){
			max_label = k_feasible_label[Nt[t][i]];
		}
	}

	// move Nt_PI[t] into X
	for(int i = 0 ; i < Nt_PI[t].size() ; i++){
		X[t].push_back(Nt_PI[t][i]);
		
	}

	// If label of Nt[t][i] < max_label, move into X (except t, t's label is 0)
	// Otherwise, move into X_bar 
	for(int i = 0 ; i < Nt[t].size() ; i++){
		if( k_feasible_label[Nt[t][i]] < max_label && k_feasible_label[Nt[t][i]] > 0){
			X[t].push_back(Nt[t][i]);
		}
		else{
			X_bar[t].push_back(Nt[t][i]);
		}	
	}

	height = max_label - 1;


	node_cut_size = calculate_cut_size(t);

	int max_size;
	int max_index;
	while( node_cut_size > K ){

		// move the node with largest fanin size in X_bar to X
		max_size = 0;
		max_index = 0;
		int found = 0;
		for(int i = 0 ; i < X_bar[t].size() ; i++){

			int size = gate[gate_index[X_bar[t][i]]].fanin_size;

			if( size > max_size){
				if( X_bar[t][i] != t){
					found = 1;
					max_size = size;
					max_index = i;
				}
				
			}
		}

		if(found){
			if(label_flag){
				//cout<< "Move " << X_bar[t][min_index] << "(" << min_size << ") to X"<< endl;
			}
			// minimum label of node in X_bar is the maximum label of nodes in X
			height = max_label;

			X[t].push_back(X_bar[t][max_index]);
			X_bar[t].erase(X_bar[t].begin() + max_index);
			
			node_cut_size = calculate_cut_size(t);	
			//cout << node_cut_size << endl;
		}
		else{
			cout<<"[Error] Find cut fail "<<endl;
			break;
		}
	}
	
	int p = calculate_p(t,index);
	// height = p - 1 exist, label of t is p
	// Otherwise, label of t is p + 1
	//cout<< "p = " << p <<endl;
	//cout<< "height = " << height << endl;
	if( height < p ){
		// t can be group with nodes in X_Bar
		//cout<< "height = p - 1" << endl;
		k_feasible_label[t] = p;
	}
	else{
		// t is the next group
		//cout<< "height = p " << endl;
		k_feasible_label[t] = p + 1;
	}

	if(label_flag){
		cout<< "nodes in X : ";
		for( int i = 0 ; i < X[t].size() ; i++){
			cout<< X[t][i] << " ";
		}
		cout<< endl;
		cout<< "nodes in X_bar : ";
		for( int i = 0 ; i < X_bar[t].size() ; i++){
			cout<< X_bar[t][i] << " ";
		}
		cout<< endl;
		cout<< "inputs to X_bar : ";
		for( int i = 0 ; i < inputs_of_X_bar[t].size() ; i++){
			cout<< inputs_of_X_bar[t][i] << " ";
		}
		cout<< endl;
		cout<< "cut size : "<< node_cut_size << endl;	
	}

	//cout<< t << "'s label : " <<k_feasible_label[t] << endl << endl; 
	
	
}

int calculate_p(string t, int index){


	int max_label = 0;
	for(int j = 0 ; j < gate[index].fanin_size ; j++){
		string fanin_node = gate[index].fanin[j];

		if(k_feasible_label[fanin_node] > max_label){
			max_label = k_feasible_label[fanin_node];
		}
	}

	return max_label;

}

int calculate_cut_size(string t){

	// compute cut size
	// calculate how many nodes in X has connection to nodes in X_bar
	
	// create a record for nodes in X
	// 0 : the node has no connection to nodes in X_bar
	// 1 : the node has one or more connection to nodes in X_bar
	record.clear();
	inputs_of_X_bar[t].clear();
	for(int i = 0 ; i < X[t].size() ; i++){
		record[X[t][i]] = 0;
	}
	
	// test every gates which fanout is in X_bar
	for(int i = 0 ; i < X_bar[t].size() ; i++){

		string X_bar_node = X_bar[t][i];

		// j is the index of the gate which fanout is X_bar[i]
		int j = 0;
		while(gate[j].fanout.compare(X_bar_node) != 0) {
			j++;
			if( j == gate.size() ){
				cout<< "[Error] No matching gate for node in X bar " << endl;
				return 0;
			}
		}
		//cout<< "current gate's fanout " << X_bar_node << endl;
		// check the every fanin_node is in X or not
		for(int k = 0 ; k < gate[j].fanin_size ; k++){
			string fanin_node = gate[j].fanin[k];
			//cout<< "test fanin "<< fanin_node << endl;
			for(int l = 0 ; l < X[t].size() ; l++){
				// found
				if(fanin_node == X[t][l]){
					//cout<< fanin_node << " is in X " << endl;
					record[fanin_node]++;
				}

			}
		}
	}

	// calculate how many nodes in X has connection with nodes on X_bar
	int count = 0;
	for(int i = 0 ; i < X[t].size() ; i++){

		if(record[X[t][i]] > 0){
			if(count <= K){
				inputs_of_X_bar[t].push_back(X[t][i]);
			}
			
			count++;
		}
	}
	//cout<< "count = " << count << endl;
	return count;
}

void mapping_phase(void){

	//clock_t start, finish;
	//double duration;

	//cout<< "Start mapping phase"<<endl;
	//start = clock();

	
	for(int i = 0 ; i < G_output.size() ; i++){
		
		string current_node = G_output[i];
		generate_LUT(current_node);
		
	}

	//cout<< "Mapping phase finished"<<endl;
	//finish = clock();
	//cout<< "Time duration : " << (double)(finish - start) / CLOCKS_PER_SEC << " sec" << endl;
	//cout<< "====================================================================="  << endl;

}

void generate_LUT(string v){

	LUT my_LUT(K, v);
	//cout<< "LUT's output : " << v <<endl;
	
	// put nodes in X_bar into LUT
	my_LUT.node.insert(my_LUT.node.begin(), X_bar[v].begin() , X_bar[v].end());
	
	// put the inputs_of_X_bar into LUT's input
	
	my_LUT.input.insert(my_LUT.input.begin(), inputs_of_X_bar[v].begin() , inputs_of_X_bar[v].end());
	
	// sort the sub_circuit in the LUT
	my_LUT.node_topology_sort();
	
	// mark the gate which fanout is v 
	gate[gate_index[v]].is_mapped = 1; 
	
	// generate k_map for my_LUT

	K_LUT.push_back(my_LUT);
	
	if(map_flag){
		cout<< "LUT's input : ";
		for(int i = 0 ; i < my_LUT.input.size() ; i++){
			cout<< my_LUT.input[i] << " ";
		}
		cout<< endl;
		cout<< "Nodes in LUT : ";
		for(int i = 0 ; i < my_LUT.node.size() ; i++){
			cout<< my_LUT.node[i] << " ";
		}

		cout<< endl;
		cout<< "LUT's output : " << my_LUT.output << endl<<endl;
	}

	for(int i = 0 ; i < my_LUT.input.size() ; i++){

		string input_node = my_LUT.input[i];
		//cout<< "input node : " << input_node << endl;
		// check input_node is a PI or not
		int flag = is_PI(input_node);

		// not a PI
		// the gate has not been marked
		if(flag == 0 && gate[gate_index[input_node]].is_mapped == 0){
			generate_LUT(input_node);
		}
		else{
			//cout<< input_node <<" is a PI"<< endl;
		}
	}
}

void calculate_level(void){

	max_level = 0;
	for(int i = 0 ; i < G_input.size() ; i++){
		PI_level = 0;
		propagate(G_input[i],1);
		//cout<< "PI level = " << PI_level << endl;
		if(PI_level > max_level){
			max_level = PI_level;
		}
	}

}


void propagate(string s, int level){

	for(int i = 0 ; i < K_LUT.size() ; i++){
		
		for(int j = 0 ; j < K_LUT[i].input.size() ; j++){

			if( K_LUT[i].input[j] == s){

				// if the LUT's output is PO
				if(is_PO(K_LUT[i].output) ){
					//cout<< "reach PO"<< endl;
					if(level > PI_level){
						PI_level = level;
					}
				}
				else{
					

					if(K_LUT[i].level_count < level){
						K_LUT[i].level_count = level;
						propagate(K_LUT[i].output,level + 1);
					}
					
					
				}


			}

		}
	
	}
}



void Output(char filename[]){

	//clock_t start, finish;
	//double duration;
	//cout<<"Output blif file : "<< filename <<endl;
	//start = clock();
	
	ofstream fout;
	fout.open(filename,ios::out);

	// model name
	fout << ".model " << model_name << endl;

	// inputs
	fout << ".inputs ";
	for(int i = 0 ; i < G_input.size() ; i++){
		fout << G_input[i] << " ";
	}
	fout << endl;

	// outputs
	fout << ".outputs ";
	for(int i = 0 ; i < G_output.size() ; i++){
		fout << G_output[i] << " ";
	}
	fout << endl;

	// LUT
	for(int i = 0 ; i < K_LUT.size() ; i++){
		fout << ".names ";

		for( int j = 0 ; j < K_LUT[i].input.size() ; j++){
			fout << K_LUT[i].input[j] << " ";
		}
		fout << K_LUT[i].output << endl;

		//cout << "current LUT output "<<K_LUT[i].output << endl;
		for(int j = 0 ; j < pow(2,K_LUT[i].input.size()) ; j++){

			bool pattern[K];
			int current = j;

			//cout<< "test : " << j << endl;
			for(int a = 0 ; a < K_LUT[i].input.size() ; a++){
				pattern[a] = current % 2;
				current = current / 2;
			}
			if( K_LUT[i].test_pattern(pattern) ){
				for(int a = 0 ; a < K_LUT[i].input.size() ; a++){
					fout<<pattern[a];
				}
				fout<<" 1"<<endl;
			}
	
		}
	}
	fout << ".end";
	//finish = clock();
	//cout<< "Time duration : " << (double)(finish - start) / CLOCKS_PER_SEC << " sec" << endl;
	//cout<< "====================================================================="  << endl;
}