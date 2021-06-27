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
#include "readData.h"

using namespace std;

struct insertion // construction algorithm data
{
        double cost;
        int node, k;

}; // end insertion

struct solution // solution data
{
        vector<int> route;      // route
        double cost;            // objective funcion

}; // end solution

void printData();
solution start_solution();                                    // Create a sub-route with 3 vertices
void construction(solution &temp);                                      // Create a initial route
solution perturb(solution start, solution ils);               // Perturb algorithm
solution gils_rvnd(int i_max, int i_ils);                     // GILS-RVND algorithm
solution rvnd(solution temp);                                 // RVND algorithm
solution rvnd_aux(int chosen, solution temp);             // aux for RVND neighborhood structure selection
solution swap(solution temp);                             // swap algorithm for neighborhood structure
solution reinsertion(solution temp);                      // reinsertion algorithm for neighborhood structure
solution or_Opt_2(solution temp);                         // Or-opt2 algorithm for neighborhood structure
solution or_Opt_3(solution temp);                         // Or-opt3 algorithm for neighborhood structure
solution two_Opt(solution temp);                          // 2-opt algorithm for neighborhood structure
bool compare_insertion(insertion a, insertion b);             // aux for sort Insertion struct
double verification (solution route);                         //

double ** matrizAdj;                  // adjacency matrix
int dimension;                        // number of vertices
solution objective;                   // solution
int perturb_type;

vector<int> instance_vertices;

double totalSI;
double totalSwap;
double totalOrOpt;
double totalOrOpt2;
double totalOrOpt3;
double totalTwoOpt;

int improve;

int main(int argc, char** argv)
{
        auto start = chrono::high_resolution_clock::now();

        srand(time(NULL));
        
        readData(argc, argv, &dimension, &matrizAdj);

        //printData();

        int i_max;
        int i_ils;

        if(dimension <= 200)
        {
                i_max = 5;

        }
        else
        {
                i_max = 10;

        } // end if/else

        if(dimension >= 150)
        {
                i_ils = dimension / 2;
        }
        else
        {
                i_ils = dimension;

        } // end if/else

        cout << endl << "---------------------------------------------------------------" << endl;

        cout << endl << "dimension: " << dimension << endl;
        cout << endl << "i_max: " << i_max << " | i_ils: " << i_ils << endl;

        cout << endl << "---------------------------------------------------------------" << endl;

        if(dimension < 75) {

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

        cout << endl << "---------------------------------------------------------------" << endl; cout << endl;

        cout << "Construction";
        cout << setfill(' ') << setw(15) << "Total time: " << fixed << setprecision(9) << totalSI << endl;

        cout << "Swap";
        cout << setfill(' ') << setw(23) << "Total time: " << fixed << setprecision(9) << totalSwap << endl;

        cout << "OrOpt-1";
        cout << setfill(' ') << setw(20) << "Total time: " << fixed << setprecision(9) << totalOrOpt << endl;

        cout << "OrOpt-2";
        cout << setfill(' ') << setw(20) << "Total time: " << fixed << setprecision(9) << totalOrOpt2 << endl;

        cout << "OrOpt-3";
        cout << setfill(' ') << setw(20) << "Total time: " << fixed << setprecision(9) << totalOrOpt3 << endl;

        cout << "TwoOpt";
        cout << setfill(' ') << setw(21) << "Total time: " << fixed << setprecision(9) << totalTwoOpt << endl;

        cout << endl << "Total time: " << fixed << setprecision(9) << time_taken << " sec" << endl;

        if(dimension > 9) n_zero = 2;
        if(dimension > 99) n_zero = 3;
        if(dimension > 999) n_zero = 4;

        cout << endl << "---------------------------------------------------------------" << endl;
        cout << endl << "Opt route: " << endl << endl;
        for(long unsigned int i = 0, j = 0; i < objective.route.size(); i++)
        {
                cout << setfill('0') << setw(n_zero) << objective.route[i] << " ";

                j++;

                if(j == 15)
                {
                        cout << endl;
                        j = 0;

                } // end if

        } // end for

        cout << endl << endl << "---------------------------------------------------------------" << endl;

        cout << endl << "cost: " << fixed << setprecision(0) << objective.cost;

        cout << " | cost (verification): " << fixed << setprecision(0) << ver << endl;

        cout << endl << "---------------------------------------------------------------" << endl; cout << endl;

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
        auto init = chrono::high_resolution_clock::now();

        temp = start_solution(); // start with 3 vertices

        while(!instance_vertices.empty())
        {
                vector<insertion> insert((temp.route.size() - 1) * instance_vertices.size());

                for(long unsigned int i = 0, j = 0; i < (temp.route.size() - 1); i++)
                {
                        for(long unsigned int k = 0; k < instance_vertices.size(); k++, j++)
                        {
                                insert[j].node = i + 1;
                                insert[j].k = k;
                                insert[j].cost = matrizAdj[temp.route[i]][instance_vertices[k]]
                                                 + matrizAdj[temp.route[i + 1]][instance_vertices[k]]
                                                 - matrizAdj[temp.route[i]][temp.route[i + 1]];
                        } // end for

                } // end for

                sort(insert.begin(), insert.end(), compare_insertion);

                int rand_value = 1 + rand() % 9;
                int chosen = rand() % rand_value;

                temp.cost += insert[chosen].cost;
                temp.route.insert(temp.route.begin() + insert[chosen].node, instance_vertices[insert[chosen].k]);
                instance_vertices.erase(instance_vertices.begin() + insert[chosen].k);

        } // end while

        auto end = chrono::high_resolution_clock::now();

        double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - init).count();

        time_taken *= 1e-9;

        totalSI += time_taken;

} // end construction

/**
 *  Create a initial sub-route with 3 vertices
 *
 *    @return the sub-route
 */
solution start_solution()
{
        solution route;

        // fill vector with vertices in the instance
        for(int i = 1; i <= dimension; i++)
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

        // start a sub-route with 3 random vertices
        for(int i = 0; i < 2; i++)
        {
                int j = rand() % instance_vertices.size();
                route.route.insert(route.route.begin() + i+1, instance_vertices[j]);
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

        for(long unsigned int i = 1; i < temp.route.size(); i++)
        {
                for(long unsigned int j = i + 1; j < temp.route.size() - 1; j++)
                {
                        if(j == (i + 1))
                        {
                                neighborhood.cost = matrizAdj[temp.route[i]][temp.route[j + 1]]
                                                    + matrizAdj[temp.route[j]][temp.route[i - 1]]
                                                    - matrizAdj[temp.route[i]][temp.route[i - 1]]
                                                    - matrizAdj[temp.route[j]][temp.route[j + 1]];

                        }
                        else
                        {
                                neighborhood.cost = matrizAdj[temp.route[i]][temp.route[j + 1]]
                                                    + matrizAdj[temp.route[i]][temp.route[j - 1]]
                                                    + matrizAdj[temp.route[j]][temp.route[i + 1]]
                                                    + matrizAdj[temp.route[j]][temp.route[i - 1]]
                                                    - matrizAdj[temp.route[i]][temp.route[i + 1]]
                                                    - matrizAdj[temp.route[i]][temp.route[i - 1]]
                                                    - matrizAdj[temp.route[j]][temp.route[j + 1]]
                                                    - matrizAdj[temp.route[j]][temp.route[j - 1]];

                        } // end if/else

                        if(neighborhood.cost < best_neighborhood.cost)
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

        for(long unsigned int i = 1; i < (temp.route.size() - 1); i++)
        {
                for(long unsigned int j = 1; j < (temp.route.size()); j++)
                {
                        if(j == i || j == (i+1) || j == (i-1)) continue;

                        neighborhood.cost = matrizAdj[temp.route[i-1]][temp.route[i+1]]
                                            + matrizAdj[temp.route[j-1]][temp.route[i]]
                                            + matrizAdj[temp.route[i]][temp.route[j]]
                                            - matrizAdj[temp.route[i-1]][temp.route[i]]
                                            - matrizAdj[temp.route[i]][temp.route[i+1]]
                                            - matrizAdj[temp.route[j-1]][temp.route[j]];

                        if(neighborhood.cost < best_neighborhood.cost)
                        {
                                best_neighborhood.cost = neighborhood.cost;

                                best_i = i;

                                if(j < i)
                                {
                                        best_j = j;
                                }
                                else
                                {
                                        best_j = j - 1;

                                }  // end if/else

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

        for(long unsigned int i = 1; i < (temp.route.size() - 2); i++)
        {
                for(long unsigned int j = 1; j < (temp.route.size()); j++)
                {
                        if(j == i || j == (i+1) || j == (i+2) || j == (i-2)) continue;

                        neighborhood.cost = matrizAdj[temp.route[i-1]][temp.route[i+2]]
                                            + matrizAdj[temp.route[i]][temp.route[j-1]]
                                            + matrizAdj[temp.route[i+1]][temp.route[j]]
                                            - matrizAdj[temp.route[i-1]][temp.route[i]]
                                            - matrizAdj[temp.route[i+1]][temp.route[i+2]]
                                            - matrizAdj[temp.route[j-1]][temp.route[j]];

                        if(neighborhood.cost < best_neighborhood.cost)
                        {
                                best_neighborhood.cost = neighborhood.cost;

                                best_neighborhood.route.assign(temp.route.begin(), temp.route.end());

                                best_i = i;

                                if(j < i)
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
        best_neighborhood.route.insert(best_neighborhood.route.begin() + (best_j), temp.route[best_i+1]);
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

        for(long unsigned int i = 1; i < (temp.route.size() - 3); i++)
        {
                for(long unsigned int j = 1; j < (temp.route.size()); j++)
                {
                        if(j == i || j == (i+1) || j == (i+2) || j == (i+3) || j == (i-3)) continue;

                        neighborhood.cost = matrizAdj[temp.route[i-1]][temp.route[i+3]]
                                            + matrizAdj[temp.route[i]][temp.route[j-1]]
                                            + matrizAdj[temp.route[i+2]][temp.route[j]]
                                            - matrizAdj[temp.route[i-1]][temp.route[i]]
                                            - matrizAdj[temp.route[i+2]][temp.route[i+3]]
                                            - matrizAdj[temp.route[j-1]][temp.route[j]];

                        if(neighborhood.cost < best_neighborhood.cost)
                        {
                                best_neighborhood.cost = neighborhood.cost;

                                best_i = i;

                                if(j < i)
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
        best_neighborhood.route.insert(best_neighborhood.route.begin() + (best_j), temp.route[best_i+2]);
        best_neighborhood.route.insert(best_neighborhood.route.begin() + (best_j), temp.route[best_i+1]);
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

        for(long unsigned int i = 3; i < (temp.route.size() - 1); i++)
        {
                for(long unsigned int j = 1; j < (temp.route.size() - i); j++)
                {
                        neighborhood.cost = matrizAdj[temp.route[i + (j-1)]][temp.route[j-1]]
                                            + matrizAdj[temp.route[j]][temp.route[i + j]]
                                            - matrizAdj[temp.route[j-1]][temp.route[j]]
                                            - matrizAdj[temp.route[i + (j-1)]][temp.route[i + j]];

                        if(neighborhood.cost < best_neighborhood.cost)
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
        vector<int> n_list = {0,1,2,3,4};

        while(!n_list.empty())
        {
                int chosen = rand() % n_list.size();
                solution candidate = rvnd_aux(n_list[chosen], temp);

                if(candidate.cost < 0)
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

        switch(chosen)
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

                        cout << endl << "Iter #" << setfill('0') << setw(2) << (i + 1) << " | " << setfill('0') << setw(2) << improve << " improvements";
                        cout << " | Start Cost: " << start_cost << " | Best: " << better_solution.cost << " | Time: " << time_taken << endl;

                }

                /*
                else
                {
                        cout << endl << "Iter #" << setfill('0') << setw(2) << (i + 1) << " | 0 improvements | Time: " << time_taken << endl;

                } // end if/else
                */

        } // end for

        return better_solution;

} // end gils_rvnd

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

        if(perturb_type == 0)
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

                int length = dimension/10;

                posicaoAleatoria1 = 1 + rand() % (dimension - 2);

                if (posicaoAleatoria1  < dimension - 2 - length)
                {
                        posicaoAleatoria2 = posicaoAleatoria1 + 1 + rand() % (length-1);
                }
                else
                {
                        posicaoAleatoria2 = posicaoAleatoria1 + 1 + rand() % (dimension + 1 - posicaoAleatoria1 - 2);

                } // end if/else

                while (true)
                {
                        posicaoAleatoria3 = 1 + rand() % (dimension - 2);

                        if (posicaoAleatoria3  < dimension - 2 - length)
                        {
                                posicaoAleatoria4 = posicaoAleatoria3 + 1 + rand() % (length-1);
                        }
                        else
                        {
                                posicaoAleatoria4 = posicaoAleatoria3 + 1 + rand() % (dimension + 1 - posicaoAleatoria3 - 2);

                        } // end if/else

                        if ((posicaoAleatoria4 < posicaoAleatoria1) || (posicaoAleatoria3 > posicaoAleatoria2)) {
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
                int tampart2 =  ponto4 - ponto3 + 1;// tamanho particao 2


                for (int i = 0; i < tampart2; i++)
                {
                        start.route[ponto1 + i] = ils.route[ponto3 + i];

                } // end for

                for (int i = 0; i < ponto3-ponto2-1; i++)
                {
                        start.route[ponto1 + tampart2 + i] = ils.route[ponto2 +1 + i];

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

double verification (solution temp)
{
        double cost = 0;

        for(long unsigned int i = 0; i < (temp.route.size() - 1); i++)
        {
                cost += matrizAdj[temp.route[i]][temp.route[i+1]];

        } // end for

        return cost;

} // end verification

void printData()
{
        cout << endl << "dimension: " << dimension << endl << endl;

        for (int i = 1; i <= dimension; i++)
        {
                for (int j = 1; j <= dimension; j++)
                {
                        cout << setfill(' ') << setw(4) << matrizAdj[i][j] << " ";

                } // end For

                cout << endl;

        } // end For

} // end printData
