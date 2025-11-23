#include <iostream>
#include <time.h>

//For superior random float 
//>
#include <random>
#include <chrono>
//<

#include <vector>
#include <algorithm>

//For Matrix Market Reader
//>
#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
//<

//For OpenMP
//>
#ifdef _OPENMP
#include <omp.h>
#endif
//<

//Matrix definition parameters

#define MATRIX_DENSITY 0.1 //10% of elements will be non-null
#define MATRIX_SIZE 100 //Size of the squared matrix [final size: MATRIX_SIZE*MATRIX_SIZE]

#define MULTIPLICATION_ITERATIONS 1 //Amout of times which an SpmV algorithm gets iterated (to calculate ns)

#define ALGORITHM_ITERATIONS 10 //Amout of times which the multiplication timing calculation is repeated 

using namespace std;

enum ProgramMode {
	Debug,
	Production,
};

//Method definition
//>
void ClearScreen(bool afterInput);
void PrintMultiplicationAndTime(const std::vector<float>& multiplication_result, const double& time_elapsed, const bool logMultiplicationResult);

float SuperiorRandomFloat(float a, float b, bool limitDecimals = true);
float RandomFloat(float a, float b);
int RandomInt(int a, int b);

template<int N>
void GenerateRandomMatrix(float (&matrix)[N][N], const float density);
template<int N>
void PrintMatrix(const float (&matrix)[N][N]);
//<

class MatrixFormat {
	public: 
		vector<int> Arow;
		vector<int> Acol;
		vector<float> Aval;
};

class COOFormat : public MatrixFormat{ };

template<int N>
void GetCOOFormat(const float (&matrix)[N][N], COOFormat &coo_format);
void PrintCOOFormat(const COOFormat &coo_format);
void MultiplyCOO(const COOFormat& coo_format, std::vector<float>& randomVector, std::vector<float>& result, double& time);
			
class CSRFormat : public MatrixFormat{ };

template<int N>
void GetCSRFormat(const float (&matrix)[N][N], CSRFormat &csr_format);
void PrintCSRFormat(const CSRFormat &csr_format);
void MultiplyCSR(const CSRFormat& csr_format, std::vector<float>& randomVector, std::vector<float>& result, double& time);

bool SortByFirstAndSecond(const std::tuple<int, int, int>& a, const std::tuple<int, int, int>& b) {
	return std::get<0>(a) < std::get<0>(b) || (std::get<0>(a) == std::get<0>(b) && std::get<1>(a) < std::get<1>(b));
}

int main(int argc, char** argv) {
	srand(time(0));
	
	if (argc < 3){
		fprintf(stderr, "Usage: %s [martix-market-filename][number-of-threads]\n", argv[0]);
		exit(1);
	}else{
		const int numberOfThreads = std::stoi(argv[2]);
		omp_set_num_threads(numberOfThreads);
	}
	
	#pragma omp parallel
	{
	    int thread_id = omp_get_thread_num();
	    int num_threads = omp_get_num_threads();
	    if (thread_id == 0) {
	        cout << "Using " << num_threads << " threads for the algorithm" << endl;
	    }
	}

	
	enum ProgramMode programMode = Production; //"Debug" or "Production"
	bool executeCOO_multiplication = false;
	bool logMultiplicationResult = false;
	bool logCOOandCSR = false;
	
	COOFormat coo_format;
	CSRFormat csr_format;
	
	int M = -1, N = -1;
	
	if (programMode == Debug){
		float sparseMatrix[MATRIX_SIZE][MATRIX_SIZE];
		GenerateRandomMatrix(sparseMatrix, MATRIX_DENSITY);
		GetCOOFormat(sparseMatrix, coo_format);
		GetCSRFormat(sparseMatrix, csr_format);
	
		cout << "Matrix density: " << MATRIX_DENSITY << "\t" << "Matrix size: " << MATRIX_SIZE << endl << endl;
		
		PrintMatrix(sparseMatrix);
		
		cout << endl << "COO Format of the matrix: \t[value_count: " << coo_format.Aval.size() << "]" << endl << endl;
		
		PrintCOOFormat(coo_format);
		
		cout << endl << "CSR Format of the matrix: \t[value_count: " << csr_format.Aval.size() << "]" << endl << endl;
		
		PrintCSRFormat(csr_format);
		
		bool waitForInput = true;
		cout << endl << "Type anything and enter to clear the screen: ";
		ClearScreen(waitForInput);
	}
	else if (programMode == Production){
		int ret_code;
	    MM_typecode matcode;
	    FILE *f;
	    int nz, i, *I, *J;
	    double *val;
	
	    if (argc < 2)
		{
			fprintf(stderr, "Usage: %s [martix-market-filename][number-of-threads]\n", argv[0]);
			exit(1);
		}
	    else if ((f = fopen(argv[1], "r")) == NULL) exit(1);
	    
	    if (mm_read_banner(f, &matcode) != 0)
	    {
	        printf("Could not process Matrix Market banner.\n");
	        exit(1);
	    }
	
	    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
	            mm_is_sparse(matcode) )
	    {
	        printf("Sorry, this application does not support ");
	        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
	        exit(1);
	    }
	
	    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
	        exit(1);
	        
	    if (M != N)
	    {
	        printf("The matrix is not squared. This kind of behaviour is not supported!\n");
	        exit(1);
	    }
	
	    I = (int *) malloc(nz * sizeof(int));
	    J = (int *) malloc(nz * sizeof(int));
	    val = (double *) malloc(nz * sizeof(double));
	
	    for (i=0; i<nz; i++)
	    {
	        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
	        I[i]--;  /* adjust from 1-based to 0-based */
	        J[i]--;
	    }
	
	    if (f != stdin) fclose(f);
	    
		coo_format.Arow.resize(nz);
		coo_format.Acol.resize(nz);
		coo_format.Aval.resize(nz);
		
		//Order COO
	    std::vector<std::tuple<int, int, float>> vector;
	    for (i=0; i<nz; i++) 
    		vector.push_back(std::make_tuple(I[i], J[i], (float)val[i]));
		std::sort(vector.begin(), vector.end(), SortByFirstAndSecond);
		
		//Create COO from ordered values
		int index = 0;
		for (const auto& t : vector) {
		    coo_format.Arow[index] = std::get<0>(t);
		    coo_format.Acol[index] = std::get<1>(t);
		    coo_format.Aval[index] = std::get<2>(t);
		    index++;
		}
	        
		if (logCOOandCSR) {
			cout << endl << "COO Format of the matrix: \t[value_count: " << coo_format.Aval.size() << "]" << endl << endl;
			PrintCOOFormat(coo_format);
		}
		
		//Create CSR
		csr_format.Aval.resize(nz);
	    csr_format.Acol.resize(nz);
	    csr_format.Arow.resize(M + 1, 0);
	    
		int current_row = 0;
		index = 0;
		
		// Convert from COO to CSR
		for (int i = 0; i < nz; i++) {
		    csr_format.Aval[index] = coo_format.Aval[i];
		    csr_format.Acol[index] = coo_format.Acol[i];
		
		    while (coo_format.Arow[i] > current_row) {
		        current_row++;
		        csr_format.Arow[current_row] = index;
		    }
		
		    index++;
		}
	    csr_format.Arow[M] = nz;
	    
		if (logCOOandCSR) {
			cout << endl << "CSR Format of the matrix: \t[value_count: " << csr_format.Aval.size() << "]" << endl << endl;
			PrintCSRFormat(csr_format);
			cout << endl;
		}
	}
	
	//Declare once - use many. I use these only via reference to avoid creating garbage
	//>
	int matrixSize = M == -1 ? MATRIX_SIZE : M;
	
	if (programMode != Debug)
		cout << "Matrix size: " << M << "x" << M << endl << endl;
	
	std::vector<float> randomVector, multiplication_result;
	multiplication_result.resize(matrixSize, 0);
	randomVector.resize(matrixSize, 0);
	
	if (programMode == Debug) for (int i=0; i<matrixSize; i++) randomVector[i] = 1;
	else for (int i=0; i<matrixSize; i++) randomVector[i] = SuperiorRandomFloat(0.5f, 10);
	
	double time_elapsed;
	double cumulated_elapsed;
	std::vector<float> times_vector;
	times_vector.resize(ALGORITHM_ITERATIONS, 0);
	int target_iterations = MULTIPLICATION_ITERATIONS / matrixSize;
	if (target_iterations < 1) target_iterations = 1;
	//<
	
	if (executeCOO_multiplication){
		cout << "--- COO Multiplication ---" << endl;
		cumulated_elapsed=0;
		for (int t=0; t<target_iterations; t++){
			MultiplyCOO(coo_format, randomVector, multiplication_result, time_elapsed);
			cumulated_elapsed += time_elapsed;
		}
		time_elapsed = cumulated_elapsed / target_iterations;
		
		PrintMultiplicationAndTime(multiplication_result, time_elapsed, logMultiplicationResult);
	}
	
	cout << "--- CSR Multiplication ---" << endl;
	cout << "Iterations " << target_iterations << endl;
	
	for (int takes=0; takes<ALGORITHM_ITERATIONS + 1; takes++){
		if (target_iterations > 1) //On my local machine without more iterations the code won't show the time... on the cluster this problem didn't occur
		{
			cumulated_elapsed=0;
			for (int t=0; t<target_iterations; t++){
				MultiplyCSR(csr_format, randomVector, multiplication_result, time_elapsed);
				cumulated_elapsed += time_elapsed;
			}
			time_elapsed = cumulated_elapsed / target_iterations;
		}else
			MultiplyCSR(csr_format, randomVector, multiplication_result, time_elapsed);
		
		if (takes==0) continue; //"hot up" the cache!
		else times_vector[takes-1] = time_elapsed;
		PrintMultiplicationAndTime(multiplication_result, time_elapsed, logMultiplicationResult);
	}
	std::sort(times_vector.begin(), times_vector.end()); 
    int percentilePosition = static_cast<int>(std::ceil(0.9 * times_vector.size())) - 1;
    double percentile = times_vector[percentilePosition];
    cout << endl << "Sorted CPU times [in nanoseconds]: ";
    for (auto t : times_vector) cout << t << "\t";
    cout << endl << endl << "CSR 90/100 percentile: " << percentile / 1000000 << "ms" << endl;

	return 0;
}

//CSR format methods
//>
template<int N>
void GetCSRFormat(const float (&matrix)[N][N], CSRFormat &csr_format) {
	int itemIndex = 0;
	bool first = true;
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            if (matrix[r][c] != 0){
            	if (first) {
            		first = false;
            		csr_format.Arow.push_back(itemIndex);
				}
            	csr_format.Acol.push_back(c);
            	csr_format.Aval.push_back(matrix[r][c]);
            	itemIndex++;
			}
        }
        if (first)
            csr_format.Arow.push_back(itemIndex);
        else
        	first = true;
    }
    csr_format.Arow.push_back(csr_format.Aval.size());
}

void PrintCSRFormat(const CSRFormat &csr_format) {
	cout << "|\t";
	int valueCount = csr_format.Arow.size();
	for (int i=0; i<valueCount; i++){
		cout << csr_format.Arow[i];
        if (i == valueCount - 1) cout << "\t|" << endl;
        else cout << "\t";
	}
	
	cout << "|\t";
	for (int i=0; i < valueCount; i++)
		cout << "- \t";
	cout << "|" << endl;
	
	cout << "|\t";
	valueCount = csr_format.Acol.size();
	for (int i=0; i<valueCount; i++){
		cout << csr_format.Acol[i];
        if (i == valueCount - 1) cout << "\t|" << endl;
        else cout << "\t";
	}
	
	cout << "|\t";
	for (int i=0; i<valueCount; i++){
		cout << csr_format.Aval[i];
        if (i == valueCount - 1) cout << "\t|" << endl;
        else cout << "\t";
	}
}

void MultiplyCSR(const CSRFormat& csr_format, std::vector<float>& randomVector, std::vector<float>& result, double& time){
	fill(result.begin(), result.end(), 0);
	
    int size = csr_format.Arow.size()-1;
	auto start = std::chrono::high_resolution_clock::now();
	
	#pragma omp parallel for
	for (int row=0; row<size; row++){
		int elementIndex = csr_format.Arow[row], nextIndex = csr_format.Arow[row+1];
		float localResult = 0;
		for (int i=elementIndex; i<nextIndex; i++){
    		localResult += csr_format.Aval[i] * randomVector[row];
		}
		if (localResult != 0) result[row] = localResult;
	}
	
	auto end = chrono::high_resolution_clock::now();
	time = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
}
//<

//COO format methods
//>
template<int N>
void GetCOOFormat(const float (&matrix)[N][N], COOFormat &coo_format) {
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            if (matrix[r][c] != 0){
            	coo_format.Arow.push_back(r);
            	coo_format.Acol.push_back(c);
            	coo_format.Aval.push_back(matrix[r][c]);
			}
        }
    }
}

void PrintCOOFormat(const COOFormat &coo_format) {
	cout << "|\t";
	int valueCount = coo_format.Aval.size();
	for (int i=0; i<valueCount; i++){
		cout << coo_format.Arow[i];
        if (i == valueCount - 1) cout << "\t|" << endl;
        else cout << "\t";
	}
	
	cout << "|\t";
	for (int i=0; i<valueCount; i++){
		cout << coo_format.Acol[i];
        if (i == valueCount - 1) cout << "\t|" << endl;
        else cout << "\t";
	}
	
	cout << "|\t";
	for (int i=0; i<valueCount; i++){
		cout << coo_format.Aval[i];
        if (i == valueCount - 1) cout << "\t|" << endl;
        else cout << "\t";
	}
}

void MultiplyCOO(const COOFormat& coo_format, std::vector<float>& randomVector, std::vector<float>& result, double& time){
	fill(result.begin(), result.end(), 0);
	
    int size = coo_format.Aval.size();
	auto start = chrono::high_resolution_clock::now();
    
	for (int v=0; v<size; v++){
		int row = coo_format.Arow[v];
    	result[row] += coo_format.Aval[v] * randomVector[row];
	}
	
	auto end = chrono::high_resolution_clock::now();
	time = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
}
//<

//Matrix utility methods
//>
template<int N>
void GenerateRandomMatrix(float (&matrix)[N][N], const float density) {
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            matrix[r][c] = 0;
        }
    }
    
    int r, c;
    long target_density = density * N * N;
    for (int d=0; d<target_density; d++){
	    r = RandomInt(0, N-1);
		c = RandomInt(0, N-1);
		while (matrix[r][c] != 0)
		{
			r++;
			if (r >= N){
				c = (c + 1) % N;
				r = 0;
			}
		}
    	matrix[r][c] = SuperiorRandomFloat(0.5f, 10);
	}
}

template<int N>
void PrintMatrix(const float (&matrix)[N][N]) {
    for (int r = 0; r < N; ++r) {
		cout << "|\t";
        for (int c = 0; c < N; ++c) {
			if (matrix[r][c] != 0) cout << matrix[r][c];
			else cout << "-";
			if (c == 9) cout << "\t|" << endl;
			else cout << "\t";
        }
    }
}
//<

//Debug utility methods
//>
void ClearScreen(bool afterInput){
	if (afterInput){
		char var;
		cin >> var;
	}
	system("cls");
}

void PrintMultiplicationAndTime(const std::vector<float>& multiplication_result, const double& time_elapsed, const bool logMultiplicationResult){
	printf("Time elapsed: %e milliseconds\n", time_elapsed / 1000000); //%.0f
	if (logMultiplicationResult) {
		cout << "Multiplication result: " << endl << endl << "|\t";
		//for (int i=0; i<multiplication_result.size(); i++) cout << multiplication_result[i] << "\t";
		for (const auto& val : multiplication_result) cout << val << "\t";
		cout << "|" << endl << endl;
	}
}
//<

//Random number generation methods
//>
int RandomInt(int a, int b){
	return (rand() % (b-a+1)) + a; 
}

float RandomFloat(float a, float b) {
    float random = (float)rand() / RAND_MAX;
    float diff = b - a;
    return a + random * diff;
}

//A method for having better randomized numbers, far superior than rand()
float SuperiorRandomFloat(float a, float b, bool limitDecimals) {
    static std::mt19937 rng(
        static_cast<unsigned long long>(
            std::chrono::high_resolution_clock::now().time_since_epoch().count()
        )
    );
    std::uniform_real_distribution<float> dist(a, b);
    return limitDecimals ? round(dist(rng) * 10.0f) / 10.0f : dist(rng);
}
//<
