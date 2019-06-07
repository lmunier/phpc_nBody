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
 * Here the solution is implemented using Brute-force algorithm to have a working parallelization of the algorithm with
 * CUDA and compare with both sequential Barnes-Hut and brute-force without parallel implementation. The user can choose
 * differents parameters as define :
 * - NB_PARTICLES : number of particles
 * - PRINT : save solution in csv file, at each time-step, to display solution animation in programs (e.g. paraview)
 * - DELTA_T : time-step of each iteration
 */

/**
 * @include constants.hpp which contains all the needed project's constants/includes
 */
#include "constants.hpp"

/**
 * __global__ parallelized function to process the random generation of the particle positions. It is done to stress
 * application.
 * 
 * @param rnd values generated for particle positions in x, y and z
 * @param nb_elements to avoid process too much time the same values due to parallelization
 */
__global__ void process_pos(float* rnd, int nb_elements){
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < nb_elements) {
        /** Process random coordinate of particles to have them shifted if they stay in boundaries */
        rnd[i] *= SIDE_2f;
        rnd[i] -= SIDE_4f;

        if (abs(rnd[i] - SHIFT) < SIDE_4f)
            rnd[i] -= SHIFT;
    }
}
 
/**
 * __global__ parallelized function to process the random generation of the mass value. It is done to have value
 * between 0 and MASS_MAX.
 * 
 * @param rnd values generated for particle mass
 * @param nb_elements to avoid process too much time the same values due to parallelization
 */
__global__ void process_mass(float* rnd, int nb_elements){
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < nb_elements) {
        rnd[i] *= MASS_MAX;
    }
}

/**
 * If PRINT is defined, generate a csv file to display animation of the result in external software (e.g. paraview).
 *
 * @param particle pointer on the particle to write in csv file
 * @param iter iteration to change filename
 * @param dir to have filepath where write all the values
 */
#ifdef PRINT
__host__ void generate_file(float* particle, int iter, const std::string& dir) {
    std::ofstream csv_file;
    std::string filename = dir + "/out_" + std::to_string(iter) + ".csv";

    csv_file.open(filename);
    csv_file << "x,y,z\n";

    for (unsigned int i = 0; i < NB_PARTICLES; i++)
        csv_file << particle[DIM_3 * i] << "," << particle[DIM_3 * i + 1] << "," << particle[DIM_3 * i + 2] << std::endl;
    
    csv_file.close();
}
#endif

/**
 * __global__ parallelized function to compute the acceleration for each particle. Each particle have his own thread
 * and compute the load applied on it from each of the others particle.
 *
 * @param pos pointer on the position of each particle in device memory
 * @param acc pointer on the acceleration of each particle in device memory
 * @param mass pointer on the mass of each particle in device memory
 */
__global__ void kernel_compute_acc(float* pos, float* acc, float* mass) {
    float r[DIM_3] = {0.0f};
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int idx = DIM_3 * i;

    if ( i < NB_PARTICLES ) {
        for (unsigned int j = 0; j < NB_PARTICLES; j++) {
            if (i != j) {
                r[0] = pos[DIM_3 * j] - pos[idx];
                r[1] = pos[DIM_3 * j + 1] - pos[idx + 1];
                r[2] = pos[DIM_3 * j + 2] - pos[idx + 2];

                float distSqrd = sqrtf(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
                float s = G * mass[j] / distSqrd;

                acc[idx] += r[0] * s;
                acc[idx + 1] += r[1] * s;
                acc[idx + 2] += r[2] * s;
            }
        }
    }
}

/**
 * __global__ parallelized function to compute the update position and velocity for each particle. Each particle have his
 * own thread and compute the load applied on it from each of the others particle.
 *
 * @param pos pointer on the position of each particle in device memory
 * @param vel pointer on the velocity of each particle in device memory
 * @param acc pointer on the acceleration of each particle in device memory
 * @param mass pointer on the mass of each particle in device memory
 */
__global__ void kernel_update_pos_vel(float* pos, float* vel, float* acc, float* mass) {
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int idx = DIM_3 * i;

    if (i < NB_PARTICLES) {
        /** Update velocity */
        vel[idx] += acc[idx] * DELTA_T;
        vel[idx + 1] += acc[idx + 1] * DELTA_T;
        vel[idx + 2] += acc[idx + 2] * DELTA_T;

        /** Update position */
        pos[idx] += vel[idx] * DELTA_T;
        pos[idx + 1] += vel[idx + 1] * DELTA_T;
        pos[idx + 2] += vel[idx + 2] * DELTA_T;

        /** Reset acceleration */
        acc[idx] = 0.0f;
        acc[idx + 1] = 0.0f;
        acc[idx + 2] = 0.0f;
    }
}

/**
 * __host__ function that call both comput_acc and update_pos_vel functions. It Allocate also the number of block needed
 * to compute all the values. A synchronization function is also called to avoid having unfinished computation on a block
 * before calling the second function.
 *
 * @param pos pointer on the position of each particle in device memory
 * @param vel pointer on the velocity of each particle in device memory
 * @param acc pointer on the acceleration of each particle in device memory
 * @param mass pointer on the mass of each particle in device memory
 */
__host__ void update_particles(float* pos, float* vel, float* acc, float* mass) {
    int threadsPerBlock = min(MAX_THREADS, DIM_3 * NB_PARTICLES);
    int blocksPerGrid = (DIM_3 * NB_PARTICLES + threadsPerBlock - 1) / threadsPerBlock;
    kernel_compute_acc<<< blocksPerGrid, threadsPerBlock >>>(pos, acc, mass);

    cudaDeviceSynchronize();
    kernel_update_pos_vel<<< blocksPerGrid, threadsPerBlock >>>(pos, vel, acc, mass);
}

/**
 * Main function, compute time to solve problem and store size of the overall area where particle are studied.
 * In this function, all the memory is also allocated on the host and device part.
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

    // Filename to store particles
    std::string dir("output");

    // allocate host memory for pos, vel and mass
    unsigned int mem_size_pos = sizeof(float) * DIM_3 * NB_PARTICLES;
    float* h_pos = (float*) malloc(mem_size_pos);
    
    unsigned int mem_size_vel = sizeof(float) * DIM_3 * NB_PARTICLES;
    float* h_vel = (float*) malloc(mem_size_vel);
    
    unsigned int mem_size_mass = sizeof(float) * NB_PARTICLES;
    float* h_mass = (float*) malloc(mem_size_mass);
    //float flop = (float) * (NB_PARTICLES * NB_PARTICLES + NB_PARTICLES) * DIM_3;
    
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
    unsigned int mem_size_acc = sizeof(float) * DIM_3 * NB_PARTICLES;
    float* d_acc;
    cudaMalloc((void**) &d_acc, mem_size_acc);

    // allocate host memory for the result
    float* h_acc = (float*) malloc(mem_size_acc);

    // create and start timer
    cudaEventCreate(&start);
    cudaEventRecord(start, NULL);

    // copy host memory to device
    cudaMemcpy(d_pos, h_pos, mem_size_pos, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vel, h_vel, mem_size_vel, cudaMemcpyHostToDevice);

    // Generate random numbers on device
    curandGenerator_t gen;

    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);                 /* Create pseudo-random number generator */
    curandSetPseudoRandomGeneratorSeed(gen, (unsigned long long) start);    /* Set seed */
    curandGenerateUniform(gen, d_pos, DIM_3 * NB_PARTICLES);                /* Generate floats pos on device */
    curandGenerateUniform(gen, d_mass, NB_PARTICLES);                       /* Generate floats mass on device */

    // setup execution parameters
    int threadsPerBlock = min(MAX_THREADS, DIM_3 * NB_PARTICLES);
    int blocksPerGrid = (DIM_3 * NB_PARTICLES + threadsPerBlock - 1) / threadsPerBlock;
    process_pos<<< blocksPerGrid, threadsPerBlock >>>(d_pos, DIM_3 * NB_PARTICLES);

    blocksPerGrid = (NB_PARTICLES + threadsPerBlock - 1) / threadsPerBlock;
    process_mass<<< blocksPerGrid, threadsPerBlock >>>(d_mass, NB_PARTICLES);

    for (unsigned int k = 0; k < ITERATIONS; k++) {
        update_particles(d_pos, d_vel, d_acc, d_mass);
        cudaDeviceSynchronize();

#ifdef PRINT
        // copy result from device to host
        cudaMemcpy(h_pos, d_pos, mem_size_pos, cudaMemcpyDeviceToHost);
        generate_file(h_pos, k, dir);
#endif
    }

    cudaMemcpy(h_pos, d_pos, mem_size_pos, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vel, d_vel, mem_size_vel, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_mass, d_mass, mem_size_mass, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_acc, d_acc, mem_size_acc, cudaMemcpyDeviceToHost);

    // stop and destroy timer
    cudaEventCreate(&stop);
    cudaEventRecord(stop, NULL);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&msecTotal, start, stop);

    //printf("Processing time: %f (ms), GFLOPS: %f \n", msecTotal, flop / msecTotal/ 1e+6);
    printf("Processing time: %f (ms)\n", msecTotal);

    /** Print all the parameters */
    std::cout << "-- Brut force parallelized --" << std::endl;
    std::cout << "Epsilon " << EPSILON << std::endl;
    std::cout << "Nb particles " << NB_PARTICLES << std::endl;
    std::cout << "Nb dimensions " << DIM_3 << std::endl;
    std::cout << "Side " << SIDE << std::endl;
    std::cout << "Shift " << SHIFT << std::endl;
    std::cout << "Occupation percentage " << OCCUPATION_PERC << std::endl;
    std::cout << "Maximum mass " << MASS_MAX << std::endl;
    std::cout << "Delta t " << DELTA_T << std::endl;
    std::cout << "Nb iterations " << ITERATIONS << std::endl;

    // clean up memory
    free(h_pos);
    free(h_vel);
    free(h_acc);
    free(h_mass);
    
    cudaFree(d_pos);
    cudaFree(d_vel);
    cudaFree(d_acc);
    cudaFree(d_mass);

    curandDestroyGenerator(gen);
    exit(EXIT_SUCCESS);
}