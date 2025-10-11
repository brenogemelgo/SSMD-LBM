#include "deviceFunctions.cuh"
#include "bcsFunctions.cuh"
#include "init.cuh"
#include "lbm.cuh"
#include "bcs.cuh"
#include "../helpers/hostFunctions.cuh"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Error: Usage: " << argv[0] << " <velocity set> <ID>\n";
        return 1;
    }
    const std::string VELOCITY_SET = argv[1];
    const std::string SIM_ID       = argv[2];
    const std::string SIM_DIR = createSimulationDirectory(VELOCITY_SET, SIM_ID);

    if (setDeviceFromEnv() < 0) return 1;
    //#define BENCHMARK
    
    setDeviceFields();
    
    constexpr dim3 block(32u, 2u, 2u); 
    constexpr dim3 blockX(16u, 16u, 1u);
    constexpr dim3 blockY(16u, 16u, 1u);
    constexpr dim3 blockZ(16u, 16u, 1u);

    constexpr dim3 grid(divUp(NX, block.x),
                        divUp(NY, block.y),
                        divUp(NZ, block.z));

    constexpr dim3 gridX(divUp(NY, blockX.x),
                         divUp(NZ, blockX.y),
                         1u);

    constexpr dim3 gridY(divUp(NX, blockY.x),
                         divUp(NZ, blockY.y),
                         1u);

    constexpr dim3 gridZ(divUp(NX, blockZ.x),
                         divUp(NY, blockZ.y),
                         1u);

    constexpr size_t dynamic = 0;

    cudaStream_t queue{};
    checkCudaErrors(cudaStreamCreate(&queue));

    // =========================== INITIALIZATION =========================== //

        setFields<<<grid, block, dynamic, queue>>>(fields);
        setOilJet<<<grid, block, dynamic, queue>>>(fields);
        setWaterJet<<<grid, block, dynamic, queue>>>(fields);
        setDistros<<<grid, block, dynamic, queue>>>(fields);
    
    // ===================================================================== //

    const auto START_TIME = std::chrono::high_resolution_clock::now();
    for (int STEP = 0; STEP <= NSTEPS; ++STEP) {

        // ======================== LATTICE BOLTZMANN RELATED ======================== //

            computePhase<<<grid, block, dynamic, queue>>>(fields);
            forceStreamCollide<<<grid, block, dynamic, queue>>>(fields);

        // ========================================================================== //


        // ============================== BOUNDARY CONDITIONS ============================== //
        
            applyOilInflow<<<gridZ, blockZ, dynamic, queue>>>(fields);
            applyWaterInflow<<<gridY, blockY, dynamic, queue>>>(fields);
            applyOutflowZ<<<gridZ, blockZ, dynamic, queue>>>(fields);
            applyOutflowY<<<gridY, blockY, dynamic, queue>>>(fields);
            periodicX   <<<gridX, blockX, dynamic, queue>>>(fields);
            //periodicY   <<<gridY, blockY, dynamic, queue>>>(fields);

        // ================================================================================= //

        #if !defined(BENCHMARK)

        checkCudaErrors(cudaDeviceSynchronize());

        if (STEP % MACRO_SAVE == 0) {

            //copyAndSaveToBinary(fields.rho, PLANE, SIM_DIR, SIM_ID, STEP, "rho");
            copyAndSaveToBinary(fields.phi, PLANE, SIM_DIR, SIM_ID, STEP, "phi");
            copyAndSaveToBinary(fields.uz,  PLANE, SIM_DIR, SIM_ID, STEP, "uz");
            std::cout << "Step " << STEP << ": bins in " << SIM_DIR << "\n";

        }

        #endif
    }

    const auto END_TIME = std::chrono::high_resolution_clock::now();
    checkCudaErrorsOutline(cudaStreamDestroy(queue));

    cudaFree(fields.f);
    cudaFree(fields.g);
    cudaFree(fields.phi);
    cudaFree(fields.rho);
    cudaFree(fields.normx);
    cudaFree(fields.normy);
    cudaFree(fields.normz);
    cudaFree(fields.ux);
    cudaFree(fields.uy);
    cudaFree(fields.uz);
    cudaFree(fields.pxx);
    cudaFree(fields.pyy);
    cudaFree(fields.pzz);
    cudaFree(fields.pxy);
    cudaFree(fields.pxz);
    cudaFree(fields.pyz);

    #if defined(D_FIELDS)
    cudaFree(dfields.vorticity_mag);
    cudaFree(dfields.velocity_mag);
    #endif 

    const std::chrono::duration<double> ELAPSED_TIME = END_TIME - START_TIME;
    const uint64_t TOTAL_CELLS = static_cast<uint64_t>(NX) * NY * NZ * static_cast<uint64_t>(NSTEPS ? NSTEPS : 1);
    const double   MLUPS       = static_cast<double>(TOTAL_CELLS) / (ELAPSED_TIME.count() * 1e6);

    std::cout << "\n// =============================================== //\n";
    std::cout << "     Total execution time    : " << ELAPSED_TIME.count() << " s\n";
    std::cout << "     Performance             : " << MLUPS << " MLUPS\n";
    std::cout << "// =============================================== //\n\n";

    generateSimulationInfoFile(SIM_DIR, SIM_ID, VELOCITY_SET, MLUPS);
    getLastCudaErrorOutline("Final sync");

    return 0;
}
