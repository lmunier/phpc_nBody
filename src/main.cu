/**
 * @file main.cpp
 * @author Munier Louis
 * @date 30.05.19
 * @version 1.0
 *
 * Main file of the nBody project.
 *
 * @mainpage N-Body project
 *
 * Project in the case of the course of Parallel and High Performance Computing at EPFL - Ecole Polytechnique Fédérale
 * de Lausanne.
 *
 * @section Description
 *
 * The n-body problem aims at simulating a dynamical system of particles under the influence of physical forces. We’ll
 * restrain on the gravity field applied on celestial bodies :
 * @f[ F_{ij} = \frac{G m_i m_j (q_j - q_i)}{|| q_j - q_i ||} @f]
 *
 * where @f$ G @f$ is the gravitational constant, @f$ m_i @f$ and @f$ m_j @f$ the masses of the @f$ i @f$-th and
 * @f$j@f$-th bodies and @f$ q_i @f$ and @f$ q_j @f$ their positions.
 *
 * @section Implementation
 *
 * The solution is implemented using Barnes-Hut algorithm and quadtree/octree data-structure. The particularity is that
 * the user can choose differents parameters as define :
 * - NB_DIM : set at DIM_2 or DIM_3 to choose between 2D and 3D problem solution
 * - NB_PARTICLES : number of particles
 * - PRINT : save solution in csv file, at each time-step, to display solution animation in programs (e.g. paraview)
 * - DELTA_T : time-step of each iteration
 */

/**
 * @include constants.hpp which contains all the needed project's constants/includes
 * @include Vector.hpp custom library to have minimal vector implementation
 * @include Tree.hpp library to create a quadtree/octree data structure and interact on different cells/particles
 */
#include "constants.hpp"
#include "Particle.hpp"

#define CUDA_CALL(x) do { if((x)!=cudaSuccess) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)

#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)
  
__global__ void process_pos(float* rnd, int nb_elements){
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < nb_elements) {
        /** Generate random coordinate for a new particle and shift it if it stays in boundaries
         * to stress application */
        rnd[i] *= SIDE_2f;
        rnd[i] -= SIDE_4f;

        if (abs(rnd[i] - SHIFT) > SIDE_4f)
            rnd[i] -= SHIFT;
    }
}
  
__global__ void process_mass(float* rnd, int nb_elements){
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < nb_elements) {
        /** Generate random coordinate for a new particle and shift it if it stays in boundaries
         * to stress application */
        rnd[i] *= MASS_MAX;
    }
}

/**
 * Update the load applied to a particle by implementation a Depth First Search on the quadtree/octree data-structure.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param head pointer to the current cell of the tree
 * @param part_loaded pointer to the current particle for which the load is computed
 */
// void update_load(thrust::device_vector< Particle<Type>* >* particles) {
//     for (auto p_i : *particles) {
//         Type load = Type();

//         for (auto p_j : *particles) {
//             Type tmp = (*p_j).get(POS) - (*p_i).get(POS);
//             float d = max(tmp.norm(), EPSILON);
//             load = load + tmp * (G * (*p_i).get_mass() * (*p_j).get_mass()) / d;
//         }

//         (*p_i).set(LOAD, load);
//     }
// }

/**
 * Update position and velocity for each particle. Generate a csv file with all the position if needed.
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param root pointer on the node of the previous cell
 * @param iter current iteration of the solution
 */
// void update_particles_pos(vector< Particle<Type>* >* particles, Type dim, int iter, const string& dir){
//     for (auto p : *particles) {
//         p->update_vel_pos();

//         if (p->is_out_boundaries(dim))
//             delete p;

// #ifdef PRINT
//         generate_file(p, 1000 * iter * DELTA_T, dir);
// #endif
//     }
// }

/**
 * If PRINT is defined, generate a csv file to display animation of the result in external software (e.g. paraview).
 *
 * @tparam Type of the vector, 2D or 3D (and int, float, etc ...)
 * @param particle pointer on the particle to write in csv file
 * @param millis_time timestep to change filename and save chronology
 */
// #ifdef PRINT
// void generate_file(AbstractType<Type>* particle, int millis_time, const string& dir) {
//     ofstream csv_file;
//     string filename = dir + "/out_" + to_string(millis_time) + ".csv";

//     csv_file.open(filename, ios::app);

//     /** Check if file is empty to write the title of each column */
//     if (csv_file.tellp() == 0) {
// #if NB_DIM == DIM_2
//         csv_file << "x,y\n";
// #elif NB_DIM == DIM_3
//         csv_file << "x,y,z\n";
// #endif
//     }

//     csv_file << particle->get(POS).to_file();
//     csv_file.close();
// }
// #endif

/**
 * Main function, compute time to solve problem and store size of the overall area where particle are studied.
 *
 * @param argc default input in c++ main function
 * @param argv default input in c++ main function
 * @return success if no errors are reached
 */
int main(int argc, char *argv[]) {
    // Utilities to track the time
    cudaEvent_t start;
    cudaEvent_t stop;
    float msecTotal(0.0f);
    
    float part_mass[NB_PARTICLES] = {0.0f};
    float part_pos[NB_DIM * NB_PARTICLES] = {0.0f};
    float part_vel[NB_DIM * NB_PARTICLES] = {0.0f};
    float part_load[NB_DIM * NB_PARTICLES] = {0.0f};

    // Filename to store particles
    std::string dir("");

    if (argv[1])
        dir = argv[1];
    else
        dir = "../output";

    // allocate host memory for pos, vel and load
    unsigned int mem_size_pos = sizeof(part_pos);
    float* h_pos = (float*) malloc(mem_size_pos);
    
    unsigned int mem_size_vel = sizeof(part_vel);
    float* h_vel = (float*) malloc(mem_size_vel);
    
    unsigned int mem_size_mass = sizeof(part_mass);
    float* h_mass = (float*) malloc(mem_size_mass);
    //float flop = 2 * (float)WC * (float)HC * (float)WA;
    
    // allocate device memory
    float* d_pos;
    cudaMalloc((void**) &d_pos, mem_size_pos);
    cudaMemset(d_pos, 0.0f, mem_size_pos);

    float* d_vel;
    cudaMalloc((void**) &d_vel, mem_size_vel);
    cudaMemset(d_vel, 0.0f, mem_size_vel);

    float* d_mass;
    cudaMalloc((void**) &d_mass, mem_size_mass);
    cudaMemset(d_mass, 0.0f, mem_size_mass);

    // allocate device memory for result
    unsigned int mem_size_load = sizeof(part_load);
    float* d_load;
    cudaMalloc((void**) &d_load, mem_size_load);

    // allocate host memory for the result
    float* h_load = (float*) malloc(mem_size_load);

    dim3 threads, grid;

    // create and start timer
    cudaEventCreate(&start);
    cudaEventRecord(start, NULL);

    // copy host memory to device
    cudaMemcpy(d_pos, h_pos, mem_size_pos, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vel, h_vel, mem_size_vel, cudaMemcpyHostToDevice);

    // Generate random numbers on device
    curandGenerator_t gen;

    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);     /* Create pseudo-random number generator */
    curandSetPseudoRandomGeneratorSeed(gen, 1234ULL);           /* Set seed */
    curandGenerateUniform(gen, d_pos, NB_DIM * NB_PARTICLES);   /* Generate floats pos on device */
    curandGenerateUniform(gen, d_mass, NB_PARTICLES);   /* Generate floats mass on device */

    // setup execution parameters
    int threadsPerBlock = THREADS;
    int blocksPerGrid = (NB_DIM * NB_PARTICLES + threadsPerBlock - 1) / threadsPerBlock;
    process_pos<<< blocksPerGrid, threadsPerBlock >>>(d_pos, NB_DIM * NB_PARTICLES);

    blocksPerGrid = (NB_PARTICLES + threadsPerBlock - 1) / threadsPerBlock;
    process_mass<<< blocksPerGrid, threadsPerBlock >>>(d_mass, NB_PARTICLES);

    // copy result from device to host
    cudaMemcpy(h_pos, d_pos, mem_size_pos, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vel, d_vel, mem_size_vel, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_mass, d_mass, mem_size_mass, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_load, d_load, mem_size_load, cudaMemcpyDeviceToHost);

    // stop and destroy timer
    cudaEventCreate(&stop);
    cudaEventRecord(stop, NULL);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&msecTotal, start, stop);
    //printf("Processing time: %f (ms), GFLOPS: %f \n", msecTotal, flop / msecTotal/ 1e+6);
    printf("Processing time: %f (ms)\n", msecTotal);

    /** Print all the parameters */
    std::cout << "Brut force" << std::endl;
    std::cout << "Epsilon " << EPSILON << std::endl;
    std::cout << "Nb particles " << NB_PARTICLES << std::endl;
    std::cout << "Nb dimensions " << NB_DIM << std::endl;
    std::cout << "Side " << SIDE << std::endl;
    std::cout << "Shift " << SHIFT << std::endl;
    std::cout << "Occupation percentage " << OCCUPATION_PERC << std::endl;
    std::cout << "Maximum mass " << MASS_MAX << std::endl;
    std::cout << "Delta t " << DELTA_T << std::endl;
    std::cout << "Nb iterations " << ITERATIONS << std::endl;

    // /* Show result */
    // for(unsigned int i = 0; i < NB_PARTICLES; i++) {
    //     printf("%1.4f %1.4f %1.4f \t %1.4f \n", h_pos[i * NB_DIM], h_pos[i * NB_DIM + 1], h_pos[i * NB_DIM + 2], h_mass[i]);
    // }

    //std::cout << exec_time << std::endl;

    // clean up memory
    free(h_pos);
    free(h_vel);
    free(h_load);
    free(h_mass);
    
    cudaFree(d_pos);
    cudaFree(d_vel);
    cudaFree(d_load);
    cudaFree(d_mass);

    curandDestroyGenerator(gen);
    exit(EXIT_SUCCESS);

    return 0;
}