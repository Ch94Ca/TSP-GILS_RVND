/*
 *   GILS-RVND Algorithm for Traveling Salesman Problem
 *
 *   Author: Carlos Henrique Silva Correia de Araujo
 *   Computer Engineering - UFPB
 *
 */

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <limits>
#include <iterator>
#include <iomanip>
#include <chrono>
#include <string>
#include <sstream>
#include <omp.h>
#include "readData.h"

#define NUM_THREADS 4

using namespace std;

struct insertion // construction algorithm data
{
        double cost;
        int node, k;

}; // end insertion

struct solution // solution data
{
        vector<int> route; // route
        double cost;       // objective funcion

}; // end solution

void printData();
solution start_solution(vector<int> &instance_vertices);                        // Create a sub-route with 3 vertices
void construction(solution &temp);                // Create a initial route
solution perturb(solution start, solution ils);   // Perturb algorithm
solution gils_rvnd(int i_max, int i_ils);         // GILS-RVND algorithm
solution rvnd(solution temp);                     // RVND algorithm
solution rvnd_aux(int chosen, solution temp);     // aux for RVND neighborhood structure selection
solution swap(solution temp);                     // swap algorithm for neighborhood structure
solution reinsertion(solution temp);              // reinsertion algorithm for neighborhood structure
solution or_Opt_2(solution temp);                 // Or-opt2 algorithm for neighborhood structure
solution or_Opt_3(solution temp);                 // Or-opt3 algorithm for neighborhood structure
solution two_Opt(solution temp);                  // 2-opt algorithm for neighborhood structure
bool compare_insertion(insertion a, insertion b); // aux for sort Insertion struct
double verification(solution route);              //

double **matrizAdj; // adjacency matrix
int dimension;      // number of vertices
solution objective; // solution
int perturb_type;

double totalSI;
double totalSwap;
double totalOrOpt;
double totalOrOpt2;
double totalOrOpt3;
double totalTwoOpt;

int improve;

stringstream output_stream;

int main(int argc, char **argv)
{
        auto start = chrono::high_resolution_clock::now();

        srand(time(NULL));

        readData(argc, argv, &dimension, &matrizAdj);

        //printData();

        int i_max;
        int i_ils;

        if (dimension <= 200)
        {
                i_max = 5;
        }
        else
        {
                i_max = 10;

        } // end if/else

        if (dimension >= 150)
        {
                i_ils = dimension / 2;
        }
        else
        {
                i_ils = dimension;

        } // end if/else

        output_stream << endl
                      << "---------------------------------------------------------------" << endl;

        output_stream << endl
                      << "dimension: " << dimension << endl;
        output_stream << endl
                      << "i_max: " << i_max << " | i_ils: " << i_ils << endl;

        output_stream << endl
                      << "---------------------------------------------------------------" << endl;

        if (dimension < 75)
        {

                perturb_type = 0;
        }
        else
        {
                perturb_type = 1;

        } // end if/else

        objective = gils_rvnd(i_max, i_ils);

        double ver = verification(objective);

        auto end = chrono::high_resolution_clock::now();

        double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();

        time_taken *= 1e-9;

        int n_zero;

        output_stream << endl
                      << "---------------------------------------------------------------" << endl;
        output_stream << endl;

        output_stream << "Construction";
        output_stream << setfill(' ') << setw(15) << "Total time: " << fixed << setprecision(9) << totalSI << endl;

        output_stream << "Swap";
        output_stream << setfill(' ') << setw(23) << "Total time: " << fixed << setprecision(9) << totalSwap << endl;

        output_stream << "OrOpt-1";
        output_stream << setfill(' ') << setw(20) << "Total time: " << fixed << setprecision(9) << totalOrOpt << endl;

        output_stream << "OrOpt-2";
        output_stream << setfill(' ') << setw(20) << "Total time: " << fixed << setprecision(9) << totalOrOpt2 << endl;

        output_stream << "OrOpt-3";
        output_stream << setfill(' ') << setw(20) << "Total time: " << fixed << setprecision(9) << totalOrOpt3 << endl;

        output_stream << "TwoOpt";
        output_stream << setfill(' ') << setw(21) << "Total time: " << fixed << setprecision(9) << totalTwoOpt << endl;

        output_stream << endl
                      << "Total time: " << fixed << setprecision(9) << time_taken << " sec" << endl;

        if (dimension > 9)
                n_zero = 2;
        if (dimension > 99)
                n_zero = 3;
        if (dimension > 999)
                n_zero = 4;

        output_stream << endl
                      << "---------------------------------------------------------------" << endl;
        output_stream << endl
                      << "Opt route: " << endl
                      << endl;

        long unsigned int objective_route_size = objective.route.size();

        for (long unsigned int i = 0, j = 0; i < objective_route_size; i++)
        {
                output_stream << setfill('0') << setw(n_zero) << objective.route[i] << " ";

                j++;

                if (j == 15)
                {
                        output_stream << endl;
                        j = 0;

                } // end if

        } // end for

        output_stream << endl
                      << endl
                      << "---------------------------------------------------------------" << endl;

        output_stream << endl
                      << "cost: " << fixed << setprecision(0) << objective.cost;

        output_stream << " | cost (verification): " << fixed << setprecision(0) << ver << endl;

        output_stream << endl
                      << "---------------------------------------------------------------" << endl;
        output_stream << endl;

        string output_string = output_stream.str();

        cout << output_string;

        return 0;

} // end Main

/**
 *  Construction algorithm based on GRASP
 *
 *    Starts with a sub-route and adds the vertices
 *    using the cheapest insertion method
 *    randomly selecting the best until
 *    all vertices are inserted in the route
 *
 *    @return the initial solution
 */
void construction(solution &temp)
{       
        vector<int> instance_vertices;

        temp = start_solution(instance_vertices); // start with 3 vertices

        long unsigned int instance_vertices_size = instance_vertices.size();
        long unsigned int temp_route_size = temp.route.size();

        while (instance_vertices_size)
        {
                vector<insertion> insert((temp_route_size - 1) * instance_vertices_size);

                for (long unsigned int i = 0, j = 0; i < (temp_route_size - 1); i++)
                {
                        for (long unsigned int k = 0; k < instance_vertices_size; k++, j++)
                        {
                                insert[j].node = i + 1;
                                insert[j].k = k;
                                insert[j].cost = matrizAdj[temp.route[i]][instance_vertices[k]] + matrizAdj[temp.route[i + 1]][instance_vertices[k]] - matrizAdj[temp.route[i]][temp.route[i + 1]];
                        } // end for

                } // end for

                sort(insert.begin(), insert.end(), compare_insertion);

                int rand_value = 1 + rand() % 9;
                int chosen = rand() % rand_value;

                temp.cost += insert[chosen].cost;
                temp.route.insert(temp.route.begin() + insert[chosen].node, instance_vertices[insert[chosen].k]);
                instance_vertices.erase(instance_vertices.begin() + insert[chosen].k);
                instance_vertices_size = instance_vertices.size();
                temp_route_size = temp.route.size();

        } // end while

} // end construction

/**
 *  Create a initial sub-route with 3 vertices
 *
 *    @return the sub-route
 */
solution start_solution(vector<int> &instance_vertices)
{
        solution route;

        // fill vector with vertices in the instance
        for (int i = 1; i <= dimension; i++)
        {
                instance_vertices.push_back(i);

        } // end for

        // define a random start vertex
        int i = rand() % instance_vertices.size();
        //int i = 0;
        int start = instance_vertices[i];
        route.route.push_back(start);
        route.route.push_back(start);
        instance_vertices.erase(instance_vertices.begin() + i);
        int current = start;

        int j;

        // start a sub-route with 3 random vertices
        for (int i = 0; i < 2; i++)
        {
                j = rand() % instance_vertices.size();
                route.route.insert(route.route.begin() + i + 1, instance_vertices[j]);
                route.cost += matrizAdj[current][instance_vertices[j]];
                current = instance_vertices[j];
                instance_vertices.erase(instance_vertices.begin() + j);

        } // end for

        // sum to the obj. function distance to start vertex
        route.cost += matrizAdj[current][start];

        return route;

} // end start_solution

/**
 *  Swap algorithm for neighborhood structure
 *
 *    Exchange two vertices in the route and calculate
 *    the cost for the objective function
 *
 *    @param struct with a solution
 *    @return the best neighborhood found
 */
solution swap(solution temp)
{
        auto init = chrono::high_resolution_clock::now();

        solution neighborhood, best_neighborhood;
        best_neighborhood.cost = 0;

        int best_i = 0, best_j = 0;

        long unsigned int temp_route_size = temp.route.size();

        for (long unsigned int i = 1; i < temp_route_size; i++)
        {
                for (long unsigned int j = i + 1; j < temp_route_size - 1; j++)
                {
                        if (j == (i + 1))
                        {
                                neighborhood.cost = matrizAdj[temp.route[i]][temp.route[j + 1]] + matrizAdj[temp.route[j]][temp.route[i - 1]] - matrizAdj[temp.route[i]][temp.route[i - 1]] - matrizAdj[temp.route[j]][temp.route[j + 1]];
                        }
                        else
                        {
                                neighborhood.cost = matrizAdj[temp.route[i]][temp.route[j + 1]] + matrizAdj[temp.route[i]][temp.route[j - 1]] + matrizAdj[temp.route[j]][temp.route[i + 1]] + matrizAdj[temp.route[j]][temp.route[i - 1]] - matrizAdj[temp.route[i]][temp.route[i + 1]] - matrizAdj[temp.route[i]][temp.route[i - 1]] - matrizAdj[temp.route[j]][temp.route[j + 1]] - matrizAdj[temp.route[j]][temp.route[j - 1]];

                        } // end if/else

                        if (neighborhood.cost < best_neighborhood.cost)
                        {
                                best_neighborhood.cost = neighborhood.cost;

                                best_i = i;
                                best_j = j;

                        } // end if

                } // end for

        } // end for

        best_neighborhood.route.assign(temp.route.begin(), temp.route.end());
        iter_swap(best_neighborhood.route.begin() + best_i, best_neighborhood.route.begin() + best_j);

        auto end = chrono::high_resolution_clock::now();

        double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - init).count();

        time_taken *= 1e-9;

        totalSwap += time_taken;

        return best_neighborhood;

} // end swap

/**
 *  Reinsertion algorithm for neighborhood structure
 *
 *    Remove a vertex, reinsert it in the solution and calculate
 *    the cost for the objective function
 *
 *    @param struct with a solution
 *    @return the best neighborhood found
 */
solution reinsertion(solution temp)
{
        auto init = chrono::high_resolution_clock::now();

        solution neighborhood, best_neighborhood;
        best_neighborhood.cost = 0;

        int best_i = 0, best_j = 0;

        long unsigned int temp_route_size = temp.route.size();

        for (long unsigned int i = 1; i < (temp_route_size - 1); i++)
        {
                for (long unsigned int j = 1; j < (temp_route_size); j++)
                {
                        if (j == i || j == (i + 1) || j == (i - 1))
                                continue;

                        neighborhood.cost = matrizAdj[temp.route[i - 1]][temp.route[i + 1]] + matrizAdj[temp.route[j - 1]][temp.route[i]] + matrizAdj[temp.route[i]][temp.route[j]] - matrizAdj[temp.route[i - 1]][temp.route[i]] - matrizAdj[temp.route[i]][temp.route[i + 1]] - matrizAdj[temp.route[j - 1]][temp.route[j]];

                        if (neighborhood.cost < best_neighborhood.cost)
                        {
                                best_neighborhood.cost = neighborhood.cost;

                                best_i = i;

                                if (j < i)
                                {
                                        best_j = j;
                                }
                                else
                                {
                                        best_j = j - 1;

                                } // end if/else

                        } // end if

                } // end for

        } // end for

        best_neighborhood.route.assign(temp.route.begin(), temp.route.end());
        best_neighborhood.route.erase(best_neighborhood.route.begin() + best_i);
        best_neighborhood.route.insert(best_neighborhood.route.begin() + best_j, temp.route[best_i]);

        auto end = chrono::high_resolution_clock::now();

        double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - init).count();

        time_taken *= 1e-9;

        totalOrOpt += time_taken;

        return best_neighborhood;

} // end reinsertion

/**
 *  Or-opt2 algorithm for neighborhood structure
 *
 *    Remove two adjacent vertices, reinsert they in the solution and calculate
 *    the cost for the objective function
 *
 *    @param struct with a solution
 *    @return the best neighborhood found
 */
solution or_Opt_2(solution temp)
{
        auto init = chrono::high_resolution_clock::now();

        solution neighborhood, best_neighborhood;
        best_neighborhood.cost = 0;

        int best_i = 0, best_j = 0;

        long unsigned int temp_route_size = temp.route.size();

        for (long unsigned int i = 1; i < (temp_route_size - 2); i++)
        {
                for (long unsigned int j = 1; j < temp_route_size; j++)
                {
                        if (j == i || j == (i + 1) || j == (i + 2) || j == (i - 2))
                                continue;

                        neighborhood.cost = matrizAdj[temp.route[i - 1]][temp.route[i + 2]] + matrizAdj[temp.route[i]][temp.route[j - 1]] + matrizAdj[temp.route[i + 1]][temp.route[j]] - matrizAdj[temp.route[i - 1]][temp.route[i]] - matrizAdj[temp.route[i + 1]][temp.route[i + 2]] - matrizAdj[temp.route[j - 1]][temp.route[j]];

                        if (neighborhood.cost < best_neighborhood.cost)
                        {
                                best_neighborhood.cost = neighborhood.cost;

                                best_neighborhood.route.assign(temp.route.begin(), temp.route.end());

                                best_i = i;

                                if (j < i)
                                {
                                        best_j = j;
                                }
                                else
                                {
                                        best_j = j - 2;

                                } // end if/else

                        } // end if

                } // end for

        } // end for

        best_neighborhood.route.assign(temp.route.begin(), temp.route.end());
        best_neighborhood.route.erase(best_neighborhood.route.begin() + best_i);
        best_neighborhood.route.erase(best_neighborhood.route.begin() + best_i);
        best_neighborhood.route.insert(best_neighborhood.route.begin() + (best_j), temp.route[best_i + 1]);
        best_neighborhood.route.insert(best_neighborhood.route.begin() + (best_j), temp.route[best_i]);

        auto end = chrono::high_resolution_clock::now();

        double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - init).count();

        time_taken *= 1e-9;

        totalOrOpt2 += time_taken;

        return best_neighborhood;

} // end or_Opt_2

/**
 *  Or-opt3 algorithm for neighborhood structure
 *
 *    Remove three adjacent vertices, reinsert they in the solution and calculate
 *    the cost for the objective function
 *
 *    @param struct with a solution
 *    @return the best neighborhood found
 */
solution or_Opt_3(solution temp)
{
        auto init = chrono::high_resolution_clock::now();

        solution neighborhood, best_neighborhood;
        best_neighborhood.cost = 0;

        int best_i = 0, best_j = 0;

        long unsigned int temp_route_size = temp.route.size();

        for (long unsigned int i = 1; i < (temp_route_size - 3); i++)
        {
                for (long unsigned int j = 1; j < temp_route_size; j++)
                {
                        if (j == i || j == (i + 1) || j == (i + 2) || j == (i + 3) || j == (i - 3))
                                continue;

                        neighborhood.cost = matrizAdj[temp.route[i - 1]][temp.route[i + 3]] + matrizAdj[temp.route[i]][temp.route[j - 1]] + matrizAdj[temp.route[i + 2]][temp.route[j]] - matrizAdj[temp.route[i - 1]][temp.route[i]] - matrizAdj[temp.route[i + 2]][temp.route[i + 3]] - matrizAdj[temp.route[j - 1]][temp.route[j]];

                        if (neighborhood.cost < best_neighborhood.cost)
                        {
                                best_neighborhood.cost = neighborhood.cost;

                                best_i = i;

                                if (j < i)
                                {
                                        best_j = j;
                                }
                                else
                                {
                                        best_j = j - 3;

                                } // end if/else

                        } // end if

                } // end for

        } // end for

        best_neighborhood.route.assign(temp.route.begin(), temp.route.end());

        best_neighborhood.route.erase(best_neighborhood.route.begin() + best_i);
        best_neighborhood.route.erase(best_neighborhood.route.begin() + best_i);
        best_neighborhood.route.erase(best_neighborhood.route.begin() + best_i);
        best_neighborhood.route.insert(best_neighborhood.route.begin() + (best_j), temp.route[best_i + 2]);
        best_neighborhood.route.insert(best_neighborhood.route.begin() + (best_j), temp.route[best_i + 1]);
        best_neighborhood.route.insert(best_neighborhood.route.begin() + (best_j), temp.route[best_i]);

        auto end = chrono::high_resolution_clock::now();

        double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - init).count();

        time_taken *= 1e-9;

        totalOrOpt3 += time_taken;

        return best_neighborhood;

} // end or_Opt_3

/**
 *  2-opt algorithm for neighborhood structure
 *
 *    Remove two non-adjacent arcs, reinsert they in the solution and calculate
 *    the cost for the objective function
 *
 *    @param struct with a solution
 *    @return the best neighborhood found
 */
solution two_Opt(solution temp)
{
        auto init = chrono::high_resolution_clock::now();

        solution neighborhood, best_neighborhood;
        best_neighborhood.cost = 0;

        int best_i = 0, best_j = 0;

        long unsigned int temp_route_size = temp.route.size();

        for (long unsigned int i = 3; i < (temp_route_size - 1); i++)
        {
                for (long unsigned int j = 1; j < (temp_route_size - i); j++)
                {
                        neighborhood.cost = matrizAdj[temp.route[i + (j - 1)]][temp.route[j - 1]] + matrizAdj[temp.route[j]][temp.route[i + j]] - matrizAdj[temp.route[j - 1]][temp.route[j]] - matrizAdj[temp.route[i + (j - 1)]][temp.route[i + j]];

                        if (neighborhood.cost < best_neighborhood.cost)
                        {
                                best_neighborhood.cost = neighborhood.cost;

                                best_i = i;
                                best_j = j;

                        } // end if

                } // end for

        } // end for

        best_neighborhood.route.assign(temp.route.begin(), temp.route.end());
        reverse(best_neighborhood.route.begin() + best_j, best_neighborhood.route.begin() + (best_j + best_i));

        auto end = chrono::high_resolution_clock::now();

        double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - init).count();

        time_taken *= 1e-9;

        totalTwoOpt += time_taken;

        return best_neighborhood;

} // end two_Opt

/**
 *  RVND algorithm
 *
 *    Makes a local search based on RVND
 *
 *    @param struct with a solution
 *    @return the best neighborhood found
 */
solution rvnd(solution temp)
{
        vector<int> n_list = {0, 1, 2, 3, 4};

        while (!n_list.empty())
        {
                int chosen = rand() % n_list.size();
                solution candidate = rvnd_aux(n_list[chosen], temp);

                if (candidate.cost < 0)
                {
                        temp.route = candidate.route;
                        temp.cost += candidate.cost;
                }
                else
                {
                        n_list.erase(n_list.begin() + chosen);

                } // end if/else

        } // end while

        return temp;

} // end rvnd

/**
 *  Aux for RVND neighborhood structure selection
 *
 *    Run only the randomly selected neighborhood structure
 *
 *    @param random value
 *    @return the best neighborhood from the selected structure
 */
solution rvnd_aux(int chosen, solution temp)
{
        solution n;

        switch (chosen)
        {
        case 0:
                n = swap(temp);
                break;
        case 1:
                n = reinsertion(temp);
                break;
        case 2:
                n = or_Opt_2(temp);
                break;
        case 3:
                n = or_Opt_3(temp);
                break;
        case 4:
                n = two_Opt(temp);
                break;

        } // end switch

        return n;

} // end rvnd_aux

/**
 *  GILS-RVND algorithm
 *
 *    Metaheuristic with GRASP, ILS, and RVND to avoid local minimum values
 *
 *    @param Max iterations value
 *    @param Max iterations value for ILS
 *    @return the best solution found
 */
solution gils_rvnd(int i_max, int i_ils)
{       
        vector<solution> start_solution(NUM_THREADS), ils_solution(NUM_THREADS), better_solution(NUM_THREADS);
        vector<double> start_cost(NUM_THREADS, 0.0);
        vector<int> iterILS(NUM_THREADS, 0), improvements(NUM_THREADS, 0);
        vector<stringstream> thread_buffer(NUM_THREADS);

        for(int i = 0; i < NUM_THREADS; i++)
        {
                better_solution[i].cost = numeric_limits<double>::infinity();

        } // for

        omp_set_num_threads(NUM_THREADS);
        
        #pragma omp parallel
        {
                int id;

                id = omp_get_thread_num();

                /*
                        Vari??veis que necessitam de exclus??o m??tua:

                        start_solution
                        ils_solution
                        start_cost
                        iterILS
                        improve
                        init
                        end
                        time_taken
                        output_stream
                        better_solution

                */

                for (int i = id; i < i_max; i += NUM_THREADS)
                {       
                        construction(start_solution[id]);

                        ils_solution[id] = start_solution[id];

                        start_cost[id] = start_solution[id].cost;

                        while (iterILS[id] < i_ils)
                        {
                                start_solution[id] = rvnd(start_solution[id]);

                                if (start_solution[id].cost < ils_solution[id].cost)
                                {
                                        ils_solution[id].cost = start_solution[id].cost;
                                        ils_solution[id].route = start_solution[id].route;

                                        iterILS[id] = 0;

                                        improvements[id]++;
                                }
                                else
                                {
                                        start_solution[id] = ils_solution[id];

                                } // end if/else

                                start_solution[id] = perturb(start_solution[id], start_solution[id]);

                                iterILS[id]++;

                        } // end while

                        if (ils_solution[id].cost < better_solution[id].cost)
                        {
                                better_solution[id].cost = ils_solution[id].cost;
                                better_solution[id].route = ils_solution[id].route;

                                thread_buffer[id] << endl
                                              << "Iter #" << setfill('0') << setw(2) << (i + 1) << " | " << setfill('0') << setw(2) << improvements[id] << " improvements";
                                thread_buffer[id] << " | Start Cost: " << start_cost[id] << " | Best: " << better_solution[id].cost << endl;
                        }

                } // end for

        } // end parallel zone

        solution best = better_solution[0];

        for(int i = 0; i < NUM_THREADS; i++)
        {
                if(better_solution[0].cost < best.cost)
                {
                        best = better_solution[i];

                }
                
                output_stream << thread_buffer[i].str();

        }

        return best;

} // end gils_rvnd

/*
solution gils_rvnd(int i_max, int i_ils)
{
        solution start_solution, ils_solution, better_solution;

        better_solution.cost = numeric_limits<double>::infinity();

        for(int i = 0; i < i_max; i++)
        {
                auto init = chrono::high_resolution_clock::now();

                construction(start_solution);

                ils_solution = start_solution;

                double start_cost = start_solution.cost;

                int iterILS = 0;

                while(iterILS < i_ils)
                {
                        start_solution = rvnd(start_solution);

                        if(start_solution.cost < ils_solution.cost)
                        {
                                ils_solution.cost = start_solution.cost;
                                ils_solution.route = start_solution.route;

                                iterILS = 0;

                                improve++;
                        }
                        else
                        {
                                start_solution = ils_solution;

                        } // end if/else

                        start_solution = perturb(start_solution, start_solution);

                        iterILS++;

                } // end while

                auto end = chrono::high_resolution_clock::now();

                double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - init).count();

                time_taken *= 1e-9;

                if(ils_solution.cost < better_solution.cost)
                {
                        better_solution.cost = ils_solution.cost;
                        better_solution.route = ils_solution.route;

                        output_stream << endl << "Iter #" << setfill('0') << setw(2) << (i + 1) << " | " << setfill('0') << setw(2) << improve << " improvements";
                        output_stream << " | Start Cost: " << start_cost << " | Best: " << better_solution.cost << " | Time: " << time_taken << endl;

                }

        } // end for

        return better_solution;

} // end gils_rvnd
*/

/**
 *  Perturb algorithm
 *
 *    Remove four non-adjacent arcs, reinsert they in the solution and calculate
 *    the cost for the objective function
 *
 *    @param struct with a solution
 *    @return result
 */
solution perturb(solution start, solution ils)
{
        // Perturb algorithm taken from Anand Subramanian code

        if (perturb_type == 0)
        {
                int posicaoAleatoria1, posicaoAleatoria2, posicaoAleatoria3, posicaoAleatoria4;

                posicaoAleatoria1 = 1 + rand() % (dimension - 2);
                posicaoAleatoria2 = posicaoAleatoria1 + 1;

                posicaoAleatoria3 = posicaoAleatoria1;
                posicaoAleatoria4 = posicaoAleatoria1;

                while (posicaoAleatoria3 == posicaoAleatoria1 || posicaoAleatoria3 == posicaoAleatoria2 || posicaoAleatoria4 == posicaoAleatoria1 || posicaoAleatoria4 == posicaoAleatoria2)
                {
                        posicaoAleatoria3 = 1 + rand() % (dimension - 2);
                        posicaoAleatoria4 = posicaoAleatoria3 + 1;

                } // end while

                start.route[posicaoAleatoria1] = ils.route[posicaoAleatoria3];
                start.route[posicaoAleatoria2] = ils.route[posicaoAleatoria4];
                start.route[posicaoAleatoria3] = ils.route[posicaoAleatoria1];
                start.route[posicaoAleatoria4] = ils.route[posicaoAleatoria2];
        }
        else
        {
                int posicaoAleatoria1, posicaoAleatoria2, posicaoAleatoria3, posicaoAleatoria4;

                int length = dimension / 10;

                posicaoAleatoria1 = 1 + rand() % (dimension - 2);

                if (posicaoAleatoria1 < dimension - 2 - length)
                {
                        posicaoAleatoria2 = posicaoAleatoria1 + 1 + rand() % (length - 1);
                }
                else
                {
                        posicaoAleatoria2 = posicaoAleatoria1 + 1 + rand() % (dimension + 1 - posicaoAleatoria1 - 2);

                } // end if/else

                while (true)
                {
                        posicaoAleatoria3 = 1 + rand() % (dimension - 2);

                        if (posicaoAleatoria3 < dimension - 2 - length)
                        {
                                posicaoAleatoria4 = posicaoAleatoria3 + 1 + rand() % (length - 1);
                        }
                        else
                        {
                                posicaoAleatoria4 = posicaoAleatoria3 + 1 + rand() % (dimension + 1 - posicaoAleatoria3 - 2);

                        } // end if/else

                        if ((posicaoAleatoria4 < posicaoAleatoria1) || (posicaoAleatoria3 > posicaoAleatoria2))
                        {
                                break;
                        } // end if

                } // end while

                int ponto1, ponto2, ponto3, ponto4;

                if (posicaoAleatoria1 < posicaoAleatoria3)
                {
                        ponto1 = posicaoAleatoria1;
                        ponto2 = posicaoAleatoria2;
                        ponto3 = posicaoAleatoria3;
                        ponto4 = posicaoAleatoria4;
                }
                else
                {
                        ponto1 = posicaoAleatoria3;
                        ponto2 = posicaoAleatoria4;
                        ponto3 = posicaoAleatoria1;
                        ponto4 = posicaoAleatoria2;

                } // end if/else

                // Executar a troca das duas particoes (emprestado de Bruno)

                int tampart1 = ponto2 - ponto1 + 1; // tamanho particao 1
                int tampart2 = ponto4 - ponto3 + 1; // tamanho particao 2

                for (int i = 0; i < tampart2; i++)
                {
                        start.route[ponto1 + i] = ils.route[ponto3 + i];

                } // end for

                for (int i = 0; i < ponto3 - ponto2 - 1; i++)
                {
                        start.route[ponto1 + tampart2 + i] = ils.route[ponto2 + 1 + i];

                } // end for

                for (int i = 0; i < tampart1; i++)
                {
                        start.route[ponto1 + tampart2 + ponto3 - ponto2 - 1 + i] = ils.route[ponto1 + i];

                } // end for

        } // end if/else

        start.cost = verification(start);

        return start;

} // end perturb

// aux for sort Insertion struct sort
bool compare_insertion(insertion a, insertion b)
{
        return a.cost < b.cost;

} // end compare_insertion

double verification(solution temp)
{
        double cost = 0;

        for (long unsigned int i = 0; i < (temp.route.size() - 1); i++)
        {
                cost += matrizAdj[temp.route[i]][temp.route[i + 1]];

        } // end for

        return cost;

} // end verification

void printData()
{
        output_stream << endl
                      << "dimension: " << dimension << endl
                      << endl;

        for (int i = 1; i <= dimension; i++)
        {
                for (int j = 1; j <= dimension; j++)
                {
                        output_stream << setfill(' ') << setw(4) << matrizAdj[i][j] << " ";

                } // end For

                output_stream << endl;

        } // end For

} // end printData
