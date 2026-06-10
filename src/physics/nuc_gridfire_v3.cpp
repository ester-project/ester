#ifndef WITH_CMAKE
#include "ester-config.h"
#endif

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "physics.h"

#include "gridfire/gridfire.h"
#include "gridfire/solver/strategies/PointSolver.h"
#include "gridfire/utils/gf_omp.h"
#include "fourdst/atomic/species.h"
#include "fourdst/composition/composition.h"
#include "fourdst/composition/utils.h"

namespace {

fourdst::composition::Composition build_test_composition(const composition_map& comp,
                                                         int i,
                                                         int j)
{
    const std::vector<std::string> symbols = {
        "H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"
    };

    std::vector<double> X = {
        comp["H1"](i,j), comp["He3"](i,j), comp["He4"](i,j), comp["C12"](i,j),
        comp["N14"](i,j), comp["O16"](i,j), comp["Ne20"](i,j), comp["Mg24"](i,j)
    };

    double sum_selected = 0.0;
    for (const double x : X) sum_selected += x;

    const double remainder = 1.0 - sum_selected;
    if (remainder > 0.0) X.back() += remainder;  // temporary Mg-24 metal bucket

    return fourdst::composition::buildCompositionFromMassFractions(symbols, X);
}

gridfire::NetIn make_gridfire_netin(const composition_map& comp,
                                    const matrix& T,
                                    const matrix& rho,
                                    int i,
                                    int j,
                                    double tMax)
{
    gridfire::NetIn netIn;
    netIn.composition = build_test_composition(comp, i, j);
    netIn.temperature = T(i,j);
    netIn.density = rho(i,j);
    netIn.energy = 0.0;
    netIn.tMax = tMax;
    netIn.dt0 = 1.0e-12;
    return netIn;
}

struct GridfireDiag {
    double eps_inst = std::numeric_limits<double>::quiet_NaN();
    double dlneps_lnrho = std::numeric_limits<double>::quiet_NaN();
    double dlneps_lnT = std::numeric_limits<double>::quiet_NaN();
    double deps_drho = std::numeric_limits<double>::quiet_NaN();
    double deps_dT = std::numeric_limits<double>::quiet_NaN();
    int num_steps = -1;
};

/*
struct GridfireRuntime {
    std::unique_ptr<gridfire::policy::MainSequencePolicy> policy;
    std::unique_ptr<gridfire::policy::ConstructionResults> constructed;

    const gridfire::engine::DynamicEngine* engine = nullptr;

    std::unique_ptr<gridfire::solver::PointSolver> localSolver;
    std::unique_ptr<gridfire::solver::GridSolverContext> solverCtx;
    std::unique_ptr<gridfire::solver::GridSolver> gridSolver;

    explicit GridfireRuntime(const fourdst::composition::Composition& seedComposition)
    {

        
        policy = std::make_unique<gridfire::policy::MainSequencePolicy>(seedComposition);

        constructed = std::make_unique<gridfire::policy::ConstructionResults>(
            policy->construct()
        );

        auto& [constructed_engine, ctx_template] = *constructed;

        engine = &constructed_engine;
        


       
        std::unique_ptr<gridfire::engine::GraphEngine> base_engine = std::make_unique<gridfire::engine::GraphEngine>(seedComposition, 3); // 3 is depth parameter
        auto blob = base_engine->constructStateBlob(nullptr);

        auto* state = gridfire::engine::scratch::get_state<gridfire::engine::scratch::MultiscalePartitioningEngineViewScratchPad, true>(*blob);

        localSolver = std::make_unique<gridfire::solver::PointSolver>(*base_engine);

        gridfire::solver::PointSolverContext ctx_template(*blob);
        solverCtx = std::make_unique<gridfire::solver::GridSolverContext>(*ctx_template);

        solverCtx->zone_completion_logging = false;
        solverCtx->set_stdout_logging(false);
        solverCtx->set_detailed_logging(false);

        gridSolver = std::make_unique<gridfire::solver::GridSolver>(*base_engine, *localSolver);
    }
};*/

/*
struct GridfireRuntime {
    std::unique_ptr<gridfire::policy::MainSequencePolicy> policy;
    const gridfire::engine::DynamicEngine* engine;
    std::unique_ptr<gridfire::solver::PointSolver> localSolver;
    std::unique_ptr<gridfire::solver::GridSolverContext> solverCtx;
    std::unique_ptr<gridfire::solver::GridSolver> gridSolver;


    explicit GridfireRuntime(const fourdst::composition::Composition& seedComposition) : engine() // constructor, have initiliser lists 
    {
        // EMB speedup: construct the expensive policy/engine/solver stack once.
        policy = std::make_unique<gridfire::policy::MainSequencePolicy>(seedComposition);

        auto [engine_local,ctx_template] = policy->construct();
        engine = &engine_local;
        //engine = &std::get<0>(constructed);
        //auto& ctx_template = *std::get<1>(constructed);

        localSolver = std::make_unique<gridfire::solver::PointSolver>(*engine);
        solverCtx = std::make_unique<gridfire::solver::GridSolverContext>(ctx_template);
        solverCtx->zone_completion_logging = false;
        solverCtx->set_stdout_logging(false);
        solverCtx->set_detailed_logging(false);

        gridSolver = std::make_unique<gridfire::solver::GridSolver>(*engine, *localSolver);
    }

    GridfireRuntime(const GridfireRuntime&) = delete;
    GridfireRuntime& operator=(const GridfireRuntime&) = delete;
};*/



struct GridfireRuntime {
    std::unique_ptr<gridfire::policy::MainSequencePolicy> policy;
    std::unique_ptr<gridfire::policy::ConstructionResults> constructed;

    const gridfire::engine::DynamicEngine* engine = nullptr;

    std::unique_ptr<gridfire::solver::PointSolver> localSolver;
    std::unique_ptr<gridfire::solver::GridSolverContext> solverCtx;
    std::unique_ptr<gridfire::solver::GridSolver> gridSolver;

    explicit GridfireRuntime(const fourdst::composition::Composition& seedComposition)
        : policy(nullptr),
        constructed(nullptr),
        engine(nullptr),
        localSolver(nullptr),
        solverCtx(nullptr),
        gridSolver(nullptr)
    {
        std::cerr << "[GridFireRuntime] start\n" << std::flush;

        policy = std::make_unique<gridfire::policy::MainSequencePolicy>(seedComposition);
        std::cerr << "[GridFireRuntime] policy made\n" << std::flush;

        constructed = std::make_unique<gridfire::policy::ConstructionResults>(
            policy->construct()
        );
        std::cerr << "[GridFireRuntime] construct done\n" << std::flush;

        auto& [constructed_engine, ctx_template] = *constructed;
        engine = &constructed_engine;
        std::cerr << "[GridFireRuntime] engine pointer set\n" << std::flush;

        localSolver = std::make_unique<gridfire::solver::PointSolver>(*engine);
        std::cerr << "[GridFireRuntime] point solver made\n" << std::flush;

        solverCtx = std::make_unique<gridfire::solver::GridSolverContext>(*ctx_template);
        std::cerr << "[GridFireRuntime] grid solver ctx made\n" << std::flush;

        solverCtx->zone_completion_logging = false;
        solverCtx->set_stdout_logging(false);
        solverCtx->set_detailed_logging(false);

        gridSolver = std::make_unique<gridfire::solver::GridSolver>(*engine, *localSolver);
        std::cerr << "[GridFireRuntime] grid solver made\n" << std::flush;
    }

    void reset_context()
    {
        auto& [constructed_engine, ctx_template] = *constructed;

        solverCtx = std::make_unique<gridfire::solver::GridSolverContext>(
            *ctx_template
        );

        solverCtx->zone_completion_logging = false;
        solverCtx->set_stdout_logging(false);
        solverCtx->set_detailed_logging(false);
    }

    GridfireRuntime(const GridfireRuntime&) = delete;
    GridfireRuntime& operator=(const GridfireRuntime&) = delete;
};

/*struct GridfireRuntime {
    std::unique_ptr<gridfire::engine::GraphEngine> engine;
    std::unique_ptr<gridfire::engine::scratch::StateBlob> ctx_template;

    std::unique_ptr<gridfire::solver::PointSolver> localSolver;
    std::unique_ptr<gridfire::solver::GridSolverContext> solverCtx;
    std::unique_ptr<gridfire::solver::GridSolver> gridSolver;

    explicit GridfireRuntime(const fourdst::composition::Composition& seedComposition)
        : engine(nullptr),
          ctx_template(nullptr),
          localSolver(nullptr),
          solverCtx(nullptr),
          gridSolver(nullptr)
    {
        std::cerr << "[GridFireRuntime GraphEngine] start\n" << std::flush;

        engine = std::make_unique<gridfire::engine::GraphEngine>(
            seedComposition,
            3
        );

        ctx_template = engine->constructStateBlob(nullptr);
        localSolver = std::make_unique<gridfire::solver::PointSolver>(*engine);

        solverCtx = std::make_unique<gridfire::solver::GridSolverContext>(
            *ctx_template
        );

        solverCtx->zone_completion_logging = false;
        solverCtx->set_stdout_logging(false);
        solverCtx->set_detailed_logging(false);

        gridSolver = std::make_unique<gridfire::solver::GridSolver>(
            *engine,
            *localSolver
        );

        std::cerr << "[GridFireRuntime GraphEngine] done\n" << std::flush;
    }

    void reset_context()
    {
        solverCtx = std::make_unique<gridfire::solver::GridSolverContext>(
            *ctx_template
        );

        solverCtx->zone_completion_logging = false;
        solverCtx->set_stdout_logging(false);
        solverCtx->set_detailed_logging(false);
    }

    GridfireRuntime(const GridfireRuntime&) = delete;
    GridfireRuntime& operator=(const GridfireRuntime&) = delete;
}; */ // trying something out with out MS policy 

GridfireRuntime& gridfire_runtime(const composition_map& comp,
                                  const matrix& T,
                                  const matrix& rho,
                                  double tMax)
{
    static std::unique_ptr<GridfireRuntime> runtime;
    if (!runtime) {
        const gridfire::NetIn seed = make_gridfire_netin(comp, T, rho, 0, 0, tMax);
        runtime = std::make_unique<GridfireRuntime>(seed.composition);
    }
    return *runtime;
}

double test_pointsolver_eps_inst(GridfireRuntime& runtime,
                                 const gridfire::NetIn& netIn)
{
    double eps_inst = std::numeric_limits<double>::quiet_NaN();
    int callback_count = 0;

    auto& [constructed_engine, ctx_template] = *runtime.constructed;
    gridfire::solver::PointSolverContext point_ctx(*ctx_template);

    //gridfire::solver::PointSolverContext point_ctx(*runtime.ctx_template); // for without MS policy

    point_ctx.set_stdout_logging(true);
    point_ctx.set_detailed_logging(true);

    /*point_ctx.callback =
        [&eps_inst, &callback_count]
        (const gridfire::solver::PointSolverTimestepContext& ctx)
        {
            callback_count++;

            auto rhs_calc =
                ctx.engine.getMostRecentRHSCalculation(ctx.state_ctx);

            if (rhs_calc.has_value()) {
                eps_inst = rhs_calc->nuclearEnergyGenerationRate;
            }
        };*/

    /*point_ctx.callback =
            [&callback_count]
            (const gridfire::solver::PointSolverTimestepContext& ctx)
            {
                callback_count++;

                if (callback_count <= 5) {
                    std::cerr << "[PointSolver callback] count="
                            << callback_count
                            << " t=" << ctx.t
                            << " num_steps=" << ctx.num_steps
                            << " T9=" << ctx.T9
                            << " rho=" << ctx.rho
                            << "\n" << std::flush;
                }
            };*/


    
    point_ctx.callback =
    [&eps_inst, &callback_count]
    (const gridfire::solver::PointSolverTimestepContext& ctx)
    {
        callback_count++;

        std::cerr << "[PointSolver callback] before RHS count="
                  << callback_count
                  << "\n" << std::flush;

        auto rhs_calc =
            ctx.engine.getMostRecentRHSCalculation(ctx.state_ctx);

        std::cerr << "[PointSolver callback] after RHS count="
                  << callback_count
                  << " has_value="
                  << rhs_calc.has_value()
                  << "\n" << std::flush;

        if (rhs_calc.has_value()) {
            eps_inst = rhs_calc->nuclearEnergyGenerationRate;
        }
    };

    const gridfire::NetOut out =
        runtime.localSolver->evaluate(point_ctx, netIn);

    std::cerr << "[PointSolver test] callbacks="
              << callback_count
              << " eps_inst="
              << eps_inst
              << " NetOut.energy="
              << out.energy
              << " num_steps="
              << out.num_steps
              << "\n" << std::flush;

    return eps_inst;
}

std::vector<GridfireDiag> run_gridfire_profile(GridfireRuntime& runtime,
                                               const composition_map& comp,
                                               const matrix& T,
                                               const matrix& rho,
                                               double tMax)
{
    const int nr = T.nrows();
    const int nth = T.ncols();
    const int ncell = nr * nth;

    std::cerr << "[run_gridfire_profile] entered\n" << std::flush;
    std::cerr << "[run_gridfire_profile] nr=" << nr << " nth=" << nth << " ncell=" << ncell << "\n" << std::flush;

    std::vector<GridfireDiag> diags(ncell);
    std::vector<gridfire::NetIn> netIns;
    netIns.reserve(ncell);

    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nth; ++j) {
            netIns.push_back(make_gridfire_netin(comp, T, rho, i, j, tMax));
        }
    }

    //std::cerr << "[run_gridfire_profile] netIns built: "
    //      << netIns.size() << "\n" << std::flush;

    if (!netIns.empty()) {
        test_pointsolver_eps_inst(runtime, netIns[0]);
    }

    try {

        //runtime.reset_context();
        
        // EMB speedup: one GridSolver call for the whole profile, not one solver construction per cell.
        //std::cerr << "[run_gridfire_profile] before reset_context\n" << std::flush;
        runtime.reset_context();
        //std::cerr << "[run_gridfire_profile] after reset_context\n" << std::flush;

        //std::cerr << "[run_gridfire_profile] before evaluate\n" << std::flush;
        //const std::vector<gridfire::NetOut> netOuts =
        //    runtime.gridSolver->evaluate(*runtime.solverCtx, netIns); //replacing this with below

        std::vector<double> eps_inst_cache(
            netIns.size(),
            std::numeric_limits<double>::quiet_NaN()
        );

        std::vector<int> callback_count(
            netIns.size(),
            0
        );

        runtime.solverCtx->timestep_callbacks.clear();
        runtime.solverCtx->timestep_callbacks.resize(netIns.size());

        for (size_t kk = 0; kk < netIns.size(); ++kk) {
            runtime.solverCtx->set_callback(
                [&eps_inst_cache, &callback_count, kk]
                (const gridfire::solver::TimestepContextBase& base_ctx)
                {
                    callback_count[kk]++;

                    const auto* ctx =
                        dynamic_cast<
                            const gridfire::solver::PointSolverTimestepContext*
                        >(&base_ctx);

                    if (!ctx) return;

                    auto rhs_calc =
                        ctx->engine.getMostRecentRHSCalculation(
                            ctx->state_ctx
                        );

                    if (rhs_calc.has_value()) {
                        eps_inst_cache[kk] =
                            rhs_calc->nuclearEnergyGenerationRate;
                    }

                    if (kk == 0 && callback_count[kk] <= 5) {
                        std::cerr
                            << "[callback] kk=0 count="
                            << callback_count[kk]
                            << " rhs="
                            << rhs_calc.has_value()
                            << "\n"
                            << std::flush;
                    }
                },
                kk
            );
        }

        //std::cerr << "[run_gridfire_profile] before evaluate\n" << std::flush;
        const std::vector<gridfire::NetOut> netOuts =
            runtime.gridSolver->evaluate(*runtime.solverCtx, netIns);

        /*{
            const int ktest = 0;

            std::cerr << "[central debug] netOut energy="
                    << netOuts[ktest].energy
                    << " num_steps="
                    << netOuts[ktest].num_steps
                    << "\n" << std::flush;

            if (ktest < static_cast<int>(runtime.solverCtx->solver_workspaces.size())) {
                auto* zone_ctx =
                    dynamic_cast<gridfire::solver::PointSolverContext*>(
                        runtime.solverCtx->solver_workspaces[ktest].get()
                    );

                std::cerr << "[central debug] zone_ctx="
                        << zone_ctx
                        << "\n" << std::flush;

                if (zone_ctx) {
                    std::cerr << "[central debug] engine_ctx="
                            << zone_ctx->engine_ctx.get()
                            << "\n" << std::flush;
                }

                if (zone_ctx && zone_ctx->engine_ctx) {
                    std::cerr << "[central debug] before getNetworkReactions\n"
                            << std::flush;

                    const auto& reactions =
                        runtime.engine->getNetworkReactions(*zone_ctx->engine_ctx);

                    std::cerr << "[central debug] after getNetworkReactions n="
                            << reactions.size()
                            << "\n" << std::flush;

                    std::cerr << "[central debug] before RHS\n"
                            << std::flush;

                    auto rhs_calc =
                        runtime.engine->getMostRecentRHSCalculation(*zone_ctx->engine_ctx);

                    std::cerr << "[central debug] after RHS has_value="
                            << rhs_calc.has_value()
                            << "\n" << std::flush;

                    if (rhs_calc.has_value()) {
                        std::cerr << "[central debug] eps_inst="
                                << rhs_calc->nuclearEnergyGenerationRate
                                << "\n" << std::flush;
                    }
                }
            }
        }*/

        std::cerr
            << "[callback summary] count[0]="
            << callback_count[0]
            << " count[last]="
            << callback_count.back()
            << "\n"
            << std::flush;

        //std::cerr << "[run_gridfire_profile] after evaluate, netOuts.size()="
        //        << netOuts.size() << "\n" << std::flush;

        //std::cerr << "[post] solver_workspaces.size()="
        //        << runtime.solverCtx->solver_workspaces.size()
        //        << "\n" << std::flush;

        const int nout = std::min<int>(static_cast<int>(netOuts.size()), ncell);
        for (int k = 0; k < nout; ++k) {
            const int i = k / nth;
            const int j = k % nth;

            //std::cerr << "[post] k=" << k << "\n" << std::flush;

            GridfireDiag& diag = diags[k];
            const auto& out = netOuts[k];
            diag.num_steps = out.num_steps;
            diag.deps_dT = out.dEps_dT;
            diag.deps_drho = out.dEps_dRho;


            /*
            std::cerr << "[NetOut] k=" << k
                    << " num_steps=" << out.num_steps
                    << " energy=" << out.energy
                    << " dEps_dT=" << out.dEps_dT
                    << " dEps_dRho=" << out.dEps_dRho
                    << " eps_cache=" << eps_inst_cache[k]
                    << "\n" << std::flush;
            */

            /*
            std::cerr << "k static cast rhs_calc before\n" << std::endl;

            if (k < static_cast<int>(runtime.solverCtx->solver_workspaces.size())) {
                auto* zone_ctx = dynamic_cast<gridfire::solver::PointSolverContext*>(
                    runtime.solverCtx->solver_workspaces[k].get()
                );
                if (zone_ctx && zone_ctx->engine_ctx) {

                    std::cout << "inside rhs_calc if statement\n" << std::endl;
                    auto rhs_calc = runtime.engine->getMostRecentRHSCalculation(*zone_ctx->engine_ctx);
                    if (rhs_calc.has_value()) diag.eps_inst = rhs_calc->nuclearEnergyGenerationRate;
                }
            }

            */
            diag.eps_inst = eps_inst_cache[k];




            //std::cerr << "k static cast rhs_calc after:" << eps_inst_cache[k] << std::endl;


            if (std::isfinite(diag.eps_inst) && diag.eps_inst > 0.0) {
                diag.dlneps_lnrho = (rho(i,j) / diag.eps_inst) * diag.deps_drho;
                diag.dlneps_lnT = (T(i,j) / diag.eps_inst) * diag.deps_dT;
            }

            if (!std::isfinite(diag.eps_inst) || diag.eps_inst < 0.0 || diag.eps_inst > 1.0e30) {
                diag.eps_inst = std::numeric_limits<double>::quiet_NaN();
                diag.dlneps_lnrho = std::numeric_limits<double>::quiet_NaN();
                diag.dlneps_lnT = std::numeric_limits<double>::quiet_NaN();
                diag.deps_drho = std::numeric_limits<double>::quiet_NaN();
                diag.deps_dT = std::numeric_limits<double>::quiet_NaN();
                diag.num_steps = -998;
            }
        }



    }
    catch (const gridfire::exceptions::GridFireError& e) {
        std::cerr << "GridFire profile evaluation failed: " << e.what() << "\n";
        for (auto& diag : diags) diag.num_steps = -999;
    }
    catch (const std::exception& e) {
        std::cerr << "Non-GridFire exception during GridFire profile evaluation: " << e.what() << "\n";
        for (auto& diag : diags) diag.num_steps = -1000;
    }

    //std::cerr << "[run_gridfire_profile] returning diags\n" << std::flush;

    return diags;
}

std::string gridfire_log_filename(double tMax)
{
    std::ostringstream fname;
    fname << "gridfire_vs_simple_profile_tMax_" << std::scientific << std::setprecision(0) << tMax << ".csv";
    std::string filename = fname.str();
    std::replace(filename.begin(), filename.end(), '+', 'p');
    std::replace(filename.begin(), filename.end(), '-', 'm');
    return filename;
}

void write_gridfire_header_once(const std::string& filename)
{
    static bool header_written = false;
    if (header_written) return;

    std::ofstream log(filename, std::ios::out);
    log << "call_id,i,j,T,rho,num_steps,"
           "eps_gridfire,dlneps_lnrho_gridfire,dlneps_lnT_gridfire,"
           "deps_drho_gridfire,deps_dT_gridfire,"
           "eps_simple_CN,dlneps_lnrho_simple_CN,dlneps_lnT_simple_CN,"
           "eps_simple_ON,dlneps_lnrho_simple_ON,dlneps_lnT_simple_ON\n";
    header_written = true;
}

} // namespace

int nuc_gridfire(const composition_map& comp,
                 const matrix& T,
                 const matrix& rho,
                 nuc_struct& nuc)
{
    // EMB note: ideally this is called once at ESTER startup. Until there is a
    // dedicated ESTER GridFire init hook, guard it here.
    static bool gf_parallel_initialized = false;
    if (!gf_parallel_initialized) {
        GF_PAR_INIT();
        gf_parallel_initialized = true;
    }

    nuc_struct simple_cn;
    std::strncpy(simple_cn.name, "simple_CN", sizeof(simple_cn.name) - 1);
    simple_cn.name[sizeof(simple_cn.name) - 1] = '\0';

    const int simple_err = nuc_simple(comp, T, rho, simple_cn);
    if (simple_err != 0) {
        std::cout << "nuc_simple simple_CN returned error code " << simple_err << "\n";
        return simple_err;
    }

    nuc_struct simple_on;
    std::strncpy(simple_on.name, "simple_ON", sizeof(simple_on.name) - 1);
    simple_on.name[sizeof(simple_on.name) - 1] = '\0';

    const int simple_on_err = nuc_simple(comp, T, rho, simple_on);
    if (simple_on_err != 0) {
        std::cout << "nuc_simple simple_ON returned error code " << simple_on_err << "\n";
        return simple_on_err;
    }

    static int call_id = 0;
    ++call_id;

    const double tMax = 1; // seconds; same as v2 for direct comparison.
    const std::string filename = gridfire_log_filename(tMax);
    write_gridfire_header_once(filename);

    //GridfireRuntime& runtime = gridfire_runtime(comp, T, rho, tMax); // really would need to be moved outside of this nuc function yk, make separate "intilisation script" for external modules
    
    std::cerr << "[nuc_gridfire] before gridfire_runtime\n" << std::flush;
    GridfireRuntime& runtime = gridfire_runtime(comp, T, rho, tMax);
    std::cerr << "[nuc_gridfire] after gridfire_runtime\n" << std::flush;    
        const std::vector<GridfireDiag> diags = run_gridfire_profile(runtime, comp, T, rho, tMax);

    {
        std::ofstream log(filename, std::ios::app);
        const int nr = T.nrows();
        const int nth = T.ncols();

        for (int i = 0; i < nr; ++i) {
            for (int j = 0; j < nth; ++j) {
                const int k = i * nth + j;
                const GridfireDiag& diag = diags[k];

                log << call_id << "," << i << "," << j << ","
                    << T(i,j) << "," << rho(i,j) << ","
                    << diag.num_steps << "," << diag.eps_inst << ","
                    << diag.dlneps_lnrho << "," << diag.dlneps_lnT << ","
                    << diag.deps_drho << "," << diag.deps_dT << ","
                    << simple_cn.eps(i,j) << "," << simple_cn.dlneps_lnrho(i,j) << ","
                    << simple_cn.dlneps_lnT(i,j) << ","
                    << simple_on.eps(i,j) << "," << simple_on.dlneps_lnrho(i,j) << ","
                    << simple_on.dlneps_lnT(i,j) << "\n";
            }
            log.flush(); // does for every row, incase it crashes and logs nothing
        }
    }

    // Keep ESTER stable for now: log GridFire but return simple_CN values.
    nuc.pp = simple_cn.pp;
    nuc.cno = simple_cn.cno;
    nuc.eps = simple_cn.eps;
    nuc.dlneps_lnrho = simple_cn.dlneps_lnrho;
    nuc.dlneps_lnT = simple_cn.dlneps_lnT;

    return 0;
}
