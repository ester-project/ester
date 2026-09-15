#ifndef WITH_CMAKE
#include "ester-config.h"
#endif

#include <algorithm>
#include <chrono>
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
#include "gridfire/engine/scratchpads/engine_graph_scratchpad.h"
#include "gridfire/engine/scratchpads/utils.h"
#include <cstdlib>

namespace {

fourdst::composition::Composition build_test_composition(const composition_map& comp,
                                                         int i,
                                                         int j)
{
    //const std::vector<std::string> symbols = {
    //    "H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"
    //};

    const std::vector<std::string> symbols = {
        "H-1", "H-2", "H-3", "n-1",
        "He-3", "He-4", "C-12", "N-14",
        "O-16", "Ne-20", "Mg-24"
    };

    //
    //std::vector<double> X = {
    //    comp["H1"](i,j), comp["He3"](i,j), comp["He4"](i,j), comp["C12"](i,j),
    //    comp["N14"](i,j), comp["O16"](i,j), comp["Ne20"](i,j), comp["Mg24"](i,j)
    //};

    std::vector<double> X = {
        comp["H1"](i,j),
        0.0,              // H-2 diagnostic zero
        0.0,              // H-3 diagnostic zero
        0.0,              // n-1 diagnostic zero
        comp["He3"](i,j),
        comp["He4"](i,j),
        comp["C12"](i,j),
        comp["N14"](i,j),
        comp["O16"](i,j),
        comp["Ne20"](i,j),
        comp["Mg24"](i,j)
    };

    double sum_selected = 0.0;
    for (const double x : X) sum_selected += x;

    const double remainder = 1.0 - sum_selected;
    if (remainder > 0.0) X.back() += remainder;  // temporary Mg-24 metal "bucket", sorting out tMax memory issues then will try full composition or limited input composition to match id's exactly.

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

    // raw GridFire values retained before failed-zone cleaning
    double raw_netout_energy = std::numeric_limits<double>::quiet_NaN();
    double raw_netout_deps_drho = std::numeric_limits<double>::quiet_NaN();
    double raw_netout_deps_dT = std::numeric_limits<double>::quiet_NaN();
    int raw_num_steps = -1;
    bool rhs_available = false;

    double eps_avg = std::numeric_limits<double>::quiet_NaN(); // time-averaged energy generation: NetOut.energy / tMax

    // Exact composition passed to GridFire for this zone (after the temporary
    // Mg-24 remainder bucket has been applied).
    double Xin_H1   = std::numeric_limits<double>::quiet_NaN();
    double Xin_He3  = std::numeric_limits<double>::quiet_NaN();
    double Xin_He4  = std::numeric_limits<double>::quiet_NaN();
    double Xin_C12  = std::numeric_limits<double>::quiet_NaN();
    double Xin_N14  = std::numeric_limits<double>::quiet_NaN();
    double Xin_O16  = std::numeric_limits<double>::quiet_NaN();
    double Xin_Ne20 = std::numeric_limits<double>::quiet_NaN();
    double Xin_Mg24 = std::numeric_limits<double>::quiet_NaN();

    // ESTER composition before the temporary Mg-24 remainder bucket.
    double Xcomp_H1   = std::numeric_limits<double>::quiet_NaN();
    double Xcomp_He3  = std::numeric_limits<double>::quiet_NaN();
    double Xcomp_He4  = std::numeric_limits<double>::quiet_NaN();
    double Xcomp_C12  = std::numeric_limits<double>::quiet_NaN();
    double Xcomp_N14  = std::numeric_limits<double>::quiet_NaN();
    double Xcomp_O16  = std::numeric_limits<double>::quiet_NaN();
    double Xcomp_Ne20 = std::numeric_limits<double>::quiet_NaN();
    double Xcomp_Mg24 = std::numeric_limits<double>::quiet_NaN();

    // Composition returned by GridFire at tMax.
    double Xout_H1   = std::numeric_limits<double>::quiet_NaN();
    double Xout_He3  = std::numeric_limits<double>::quiet_NaN();
    double Xout_He4  = std::numeric_limits<double>::quiet_NaN();
    double Xout_C12  = std::numeric_limits<double>::quiet_NaN();
    double Xout_N14  = std::numeric_limits<double>::quiet_NaN();
    double Xout_O16  = std::numeric_limits<double>::quiet_NaN();
    double Xout_Ne20 = std::numeric_limits<double>::quiet_NaN();
    double Xout_Mg24 = std::numeric_limits<double>::quiet_NaN();

    // Per-profile timings. These are identical for every row belonging to
    // one nuc_gridfire call, which makes grouping by call_id straightforward.
    double evaluate_s = std::numeric_limits<double>::quiet_NaN();
    double rhs_lookup_s = std::numeric_limits<double>::quiet_NaN();
    double postprocess_s = std::numeric_limits<double>::quiet_NaN();

};

template <typename CompositionT>
double safe_gridfire_mass_fraction(
    const CompositionT& composition,
    const char* symbol)
{
    try {
        return composition.getMassFraction(symbol);
    }
    catch (const std::exception& e) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

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
    std::unique_ptr<gridfire::policy::MainSequencePolicy> policy; // need to replace with own policy
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


struct EmbeddedBucketPoint { // test to investigate tMax = 1s "spikes" MG edit: spikes have disappeared with pip install version...
    int i;
    double T;
    double rho;
};

fourdst::composition::Composition build_embedded_bucket_composition()
{
    const std::vector<std::string> symbols = {
        "H-1", "H-2", "H-3", "n-1",
        "He-3", "He-4", "C-12", "N-14",
        "O-16", "Ne-20", "Mg-24"
    };

    std::vector<double> X = {
        0.699973,        // H-1
        0.0,             // H-2
        0.0,             // H-3
        0.0,             // n-1
        3.53616e-05,     // He-3
        0.282708,        // He-4
        0.0029579,       // C-12
        0.000914679,     // N-14
        0.00830978,      // O-16
        0.00156936,      // Ne-20
        0.000506014      // Mg-24 before bucket remainder
    };

    double sum_selected = 0.0;
    for (const double x : X) sum_selected += x;

    const double remainder = 1.0 - sum_selected;
    if (remainder > 0.0) {
        X.back() += remainder;
    }

    std::cerr << std::setprecision(16)
              << "[embedded bucket] sum_selected_raw=" << sum_selected
              << " remainder_added_to_Mg24=" << std::max(0.0, remainder)
              << " X_Mg24_gridfire_input=" << X.back()
              << "\n" << std::flush;

    return fourdst::composition::buildCompositionFromMassFractions(symbols, X);
}

gridfire::NetIn make_embedded_bucket_netin(const EmbeddedBucketPoint& p,
                                           const fourdst::composition::Composition& comp,
                                           double tMax)
{
    gridfire::NetIn netIn;
    netIn.composition = comp;
    netIn.temperature = p.T;
    netIn.density = p.rho;
    netIn.energy = 0.0;
    netIn.tMax = tMax;
    netIn.dt0 = 1.0e-12;
    return netIn;
}

void run_embedded_bucket_probe_once(double tMax)
{
    static bool already_ran = false;
    if (already_ran) return;
    already_ran = true;

    std::cerr << "[embedded bucket] start tMax=" << tMax << "\n" << std::flush;

    const auto comp = build_embedded_bucket_composition();

    const std::vector<EmbeddedBucketPoint> points = {
        {150, 4373540.0, 0.100452},
        {151, 4365590.0, 0.099734},
        {152, 4341900.0, 0.0976187},
        {153, 4302920.0, 0.0942168},
        {154, 4249420.0, 0.0897005},
        {155, 4182380.0, 0.0842874},
        {156, 4103030.0, 0.0782204},
        {157, 4012790.0, 0.071748},
        {158, 3913200.0, 0.0651069},
        {159, 3805750.0, 0.0585093},
        {160, 3691790.0, 0.0521344},
        {999, 2378163.0, 0.00744586}
    };

    GridfireRuntime probe_runtime(comp);

    std::vector<gridfire::NetIn> netIns;
    netIns.reserve(points.size());
    for (const auto& p : points) {
        netIns.push_back(make_embedded_bucket_netin(p, comp, tMax));
    }

    std::cerr << "[embedded bucket] before reset_context\n" << std::flush;
    probe_runtime.reset_context();

    std::cerr << "[embedded bucket] before GridSolver evaluate n="
              << netIns.size() << "\n" << std::flush;

    const std::vector<gridfire::NetOut> netOuts =
        probe_runtime.gridSolver->evaluate(*probe_runtime.solverCtx, netIns);

    // keeping track of number of steps for evaulation
    
    long long total_steps = 0;
    int max_steps = -1;
    int max_steps_i = -1;

    for (std::size_t k = 0; k < netOuts.size(); ++k) {
        const int steps = netOuts[k].num_steps;

        if (steps > 0)
            total_steps += steps;

        if (steps > max_steps) {
            max_steps = steps;
            max_steps_i = static_cast<int>(k);
        }
    }

    std::cerr
        << " total_steps=" << total_steps
        << " max_steps=" << max_steps
        << " max_steps_i=" << max_steps_i
        << "\n";


    std::cerr << "[embedded bucket] after GridSolver evaluate netOuts.size="
              << netOuts.size()
              << " solver_workspaces.size="
              << probe_runtime.solverCtx->solver_workspaces.size()
              << "\n" << std::flush;

    const int nout = std::min<int>(
        static_cast<int>(netOuts.size()),
        static_cast<int>(points.size())
    );

    std::cout
        << "embedded_bucket_header,"
        << "i,T,rho,num_steps,eps_rhs,eps_avg_from_netout,netout_energy,"
        << "deps_drho,deps_dT\n";

    for (int k = 0; k < nout; ++k) {
        const auto& p = points[k];
        const auto& out = netOuts[k];

        double eps_rhs = std::numeric_limits<double>::quiet_NaN();

        std::cerr << "[embedded bucket] post k=" << k
                  << " i=" << p.i << "\n" << std::flush;

        if (k < static_cast<int>(probe_runtime.solverCtx->solver_workspaces.size())) {
            auto* zone_ctx =
                dynamic_cast<gridfire::solver::PointSolverContext*>(
                    probe_runtime.solverCtx->solver_workspaces[k].get()
                );

            std::cerr << "[embedded bucket] zone_ctx=" << zone_ctx;
            if (zone_ctx) {
                std::cerr << " engine_ctx=" << zone_ctx->engine_ctx.get();
            }
            std::cerr << "\n" << std::flush;

            if (zone_ctx && zone_ctx->engine_ctx) {
                std::cerr << "[embedded bucket] before RHS k=" << k << "\n" << std::flush;

                auto rhs_calc =
                    probe_runtime.engine->getMostRecentRHSCalculation(
                        *zone_ctx->engine_ctx
                    );

                std::cerr << "[embedded bucket] after RHS k=" << k
                          << " has_value=" << rhs_calc.has_value()
                          << "\n" << std::flush;

                if (rhs_calc.has_value()) {
                    eps_rhs = rhs_calc->nuclearEnergyGenerationRate;
                }
            }
        }

        std::cout
            << "embedded_bucket,"
            << p.i << ","
            << std::setprecision(12) << p.T << ","
            << std::setprecision(12) << p.rho << ","
            << out.num_steps << ","
            << std::setprecision(16) << eps_rhs << ","
            << std::setprecision(16) << (out.energy / tMax) << ","
            << std::setprecision(16) << out.energy << ","
            << std::setprecision(16) << out.dEps_dRho << ","
            << std::setprecision(16) << out.dEps_dT
            << "\n";
    }

    std::cerr << "[embedded bucket] done\n" << std::flush;
}

void run_embedded_full_profile_probe_once(double tMax, int zones)
{
    static bool already_ran = false;
    if (already_ran) return;
    already_ran = true;

    // initial profile
    
    const std::vector<EmbeddedBucketPoint> points = {
        {0, 30000000, 6.9617899999999997},
        {1, 30000000, 6.9617899999999997},
        {2, 30000000, 6.9617899999999997},
        {3, 29999900, 6.9618000000000002},
        {4, 29999800, 6.9618099999999998},
        {5, 29999500, 6.9618200000000003},
        {6, 29998900, 6.9618500000000001},
        {7, 29998000, 6.9618900000000004},
        {8, 29996800, 6.9619600000000004},
        {9, 29995000, 6.9620499999999996},
        {10, 29992600, 6.9621700000000004},
        {11, 29989700, 6.9623200000000001},
        {12, 29986000, 6.96251},
        {13, 29981700, 6.9627299999999996},
        {14, 29976700, 6.9629799999999999},
        {15, 29971100, 6.9632699999999996},
        {16, 29964900, 6.9635899999999999},
        {17, 29958200, 6.9639300000000004},
        {18, 29951100, 6.9642799999999996},
        {19, 29943800, 6.9646499999999998},
        {20, 29936500, 6.9650299999999996},
        {21, 29929300, 6.9653900000000002},
        {22, 29922400, 6.9657400000000003},
        {23, 29916000, 6.9660700000000002},
        {24, 29910200, 6.9663599999999999},
        {25, 29905200, 6.9666100000000002},
        {26, 29901200, 6.9668200000000002},
        {27, 29898300, 6.9669699999999999},
        {28, 29896400, 6.96706},
        {29, 29895800, 6.9670899999999998},
        {30, 29895800, 6.9670899999999998},
        {31, 29895200, 6.9671200000000004},
        {32, 29893400, 6.9672099999999997},
        {33, 29890300, 6.9673699999999998},
        {34, 29886000, 6.9675900000000004},
        {35, 29880400, 6.9678699999999996},
        {36, 29873500, 6.9682199999999996},
        {37, 29865300, 6.9686300000000001},
        {38, 29855900, 6.9691099999999997},
        {39, 29845100, 6.9696600000000002},
        {40, 29833100, 6.9702599999999997},
        {41, 29819900, 6.9709300000000001},
        {42, 29805600, 6.9716500000000003},
        {43, 29790200, 6.9724199999999996},
        {44, 29774000, 6.9732399999999997},
        {45, 29757100, 6.9740900000000003},
        {46, 29739700, 6.9749600000000003},
        {47, 29722000, 6.9758500000000003},
        {48, 29704200, 6.9767400000000004},
        {49, 29686700, 6.9776100000000003},
        {50, 29669700, 6.9784600000000001},
        {51, 29653600, 6.9792699999999996},
        {52, 29638500, 6.9800199999999997},
        {53, 29624700, 6.9807100000000002},
        {54, 29612600, 6.9813099999999997},
        {55, 29602400, 6.9818199999999999},
        {56, 29594200, 6.9822300000000004},
        {57, 29588200, 6.9825200000000001},
        {58, 29584600, 6.9827000000000004},
        {59, 29583300, 6.9827599999999999},
        {60, 29583300, 6.9827599999999999},
        {61, 29582100, 6.9828200000000002},
        {62, 29578400, 6.9830100000000002},
        {63, 29572400, 6.9833100000000004},
        {64, 29563900, 6.9837300000000004},
        {65, 29553000, 6.9842700000000004},
        {66, 29539800, 6.9849199999999998},
        {67, 29524300, 6.98569},
        {68, 29506600, 6.9865599999999999},
        {69, 29486900, 6.9875299999999996},
        {70, 29465200, 6.9885999999999999},
        {71, 29441800, 6.9897499999999999},
        {72, 29416800, 6.9909800000000004},
        {73, 29390400, 6.9922800000000001},
        {74, 29363000, 6.9936199999999999},
        {75, 29334800, 6.9950000000000001},
        {76, 29306100, 6.9963899999999999},
        {77, 29277400, 6.9977900000000002},
        {78, 29249000, 6.9991700000000003},
        {79, 29221300, 7.0005100000000002},
        {80, 29194600, 7.0018000000000002},
        {81, 29169500, 7.0030200000000002},
        {82, 29146200, 7.0041399999999996},
        {83, 29125100, 7.0051600000000001},
        {84, 29106700, 7.0060399999999996},
        {85, 29091200, 7.0067899999999996},
        {86, 29078800, 7.0073800000000004},
        {87, 29069800, 7.0078100000000001},
        {88, 29064300, 7.00807},
        {89, 29062500, 7.0081600000000002},
        {90, 29062500, 7.0081600000000002},
        {91, 29060700, 7.0082500000000003},
        {92, 29055200, 7.0085100000000002},
        {93, 29046100, 7.0089499999999996},
        {94, 29033400, 7.0095599999999996},
        {95, 29017200, 7.0103299999999997},
        {96, 28997700, 7.01126},
        {97, 28974900, 7.0123499999999996},
        {98, 28949100, 7.0135800000000001},
        {99, 28920400, 7.0149400000000002},
        {100, 28889000, 7.0164200000000001},
        {101, 28855300, 7.0180199999999999},
        {102, 28819600, 7.0197000000000003},
        {103, 28782300, 7.0214600000000003},
        {104, 28743600, 7.0232700000000001},
        {105, 28704100, 7.0251200000000003},
        {106, 28664300, 7.02698},
        {107, 28624600, 7.0288300000000001},
        {108, 28585400, 7.0306499999999996},
        {109, 28547500, 7.0324099999999996},
        {110, 28511200, 7.0340800000000003},
        {111, 28477000, 7.03566},
        {112, 28445500, 7.0370999999999997},
        {113, 28417200, 7.0384000000000002},
        {114, 28392400, 7.0395399999999997},
        {115, 28371600, 7.0404900000000001},
        {116, 28355100, 7.0412400000000002},
        {117, 28343100, 7.0417899999999998},
        {118, 28335800, 7.0421199999999997},
        {119, 28333300, 7.0422399999999996},
        {120, 28333300, 7.0422399999999996},
        {121, 28330900, 7.0423499999999999},
        {122, 28323600, 7.0426799999999998},
        {123, 28311500, 7.0432300000000003},
        {124, 28294600, 7.0439999999999996},
        {125, 28273200, 7.0449700000000002},
        {126, 28247300, 7.0461400000000003},
        {127, 28217200, 7.0475000000000003},
        {128, 28183200, 7.0490399999999998},
        {129, 28145500, 7.0507400000000002},
        {130, 28104500, 7.0525799999999998},
        {131, 28060600, 7.0545400000000003},
        {132, 28014200, 7.05661},
        {133, 27965800, 7.05877},
        {134, 27915900, 7.0609799999999998},
        {135, 27865200, 7.0632200000000003},
        {136, 27814100, 7.0654599999999999},
        {137, 27763400, 7.0676899999999998},
        {138, 27713600, 7.0698600000000003},
        {139, 27665300, 7.0719599999999998},
        {140, 27619400, 7.0739599999999996},
        {141, 27576200, 7.0758200000000002},
        {142, 27536600, 7.0775300000000003},
        {143, 27501000, 7.0790699999999998},
        {144, 27469900, 7.0804},
        {145, 27443700, 7.0815200000000003},
        {146, 27423000, 7.0823999999999998},
        {147, 27408000, 7.0830500000000001},
        {148, 27398900, 7.0834299999999999},
        {149, 27395800, 7.0835600000000003},
        {150, 27395800, 7.0835600000000003},
        {151, 27392800, 7.0836899999999998},
        {152, 27383600, 7.0840800000000002},
        {153, 27368500, 7.0847300000000004},
        {154, 27347500, 7.0856199999999996},
        {155, 27320700, 7.0867599999999999},
        {156, 27288500, 7.08812},
        {157, 27251200, 7.0896999999999997},
        {158, 27208900, 7.0914799999999998},
        {159, 27162300, 7.0934400000000002},
        {160, 27111600, 7.0955500000000002},
        {161, 27057400, 7.09781},
        {162, 27000400, 7.1001799999999999},
        {163, 26941000, 7.1026400000000001},
        {164, 26879900, 7.1051500000000001},
        {165, 26817900, 7.1076899999999998},
        {166, 26755600, 7.1102299999999996},
        {167, 26693800, 7.11273},
        {168, 26633300, 7.1151799999999996},
        {169, 26574900, 7.1175300000000004},
        {170, 26519200, 7.1197600000000003},
        {171, 26467100, 7.1218399999999997},
        {172, 26419300, 7.1237399999999997},
        {173, 26376400, 7.1254400000000002},
        {174, 26338900, 7.1269200000000001},
        {175, 26307500, 7.1281499999999998},
        {176, 26282600, 7.12913},
        {177, 26264600, 7.1298399999999997},
        {178, 26253700, 7.1302700000000003},
        {179, 26250000, 7.1304100000000004},
        {180, 26250000, 7.1304100000000004},
        {181, 26246300, 7.1305500000000004},
        {182, 26235400, 7.1309800000000001},
        {183, 26217200, 7.1316899999999999},
        {184, 26192000, 7.1326799999999997},
        {185, 26160000, 7.1339199999999998},
        {186, 26121500, 7.1354199999999999},
        {187, 26076800, 7.1371500000000001},
        {188, 26026400, 7.1391},
        {189, 25970700, 7.1412399999999998},
        {190, 25910400, 7.1435500000000003},
        {191, 25846000, 7.1459999999999999},
        {192, 25778200, 7.1485700000000003},
        {193, 25707800, 7.1512200000000004},
        {194, 25635500, 7.1539299999999999},
        {195, 25562200, 7.1566599999999996},
        {196, 25488700, 7.1593900000000001},
        {197, 25416000, 7.1620699999999999},
        {198, 25344800, 7.1646799999999997},
        {199, 25276100, 7.1671800000000001},
        {200, 25210800, 7.1695500000000001},
        {201, 25149700, 7.1717599999999999},
        {202, 25093700, 7.1737700000000002},
        {203, 25043400, 7.1755699999999996},
        {204, 24999700, 7.17713},
        {205, 24963000, 7.1784299999999996},
        {206, 24933900, 7.1794599999999997},
        {207, 24912900, 7.1802099999999998},
        {208, 24900100, 7.1806599999999996},
        {209, 24895800, 7.1808100000000001},
        {210, 24895800, 7.1808100000000001},
        {211, 24891600, 7.1809599999999998},
        {212, 24878800, 7.1814099999999996},
        {213, 24857600, 7.1821599999999997},
        {214, 24828200, 7.1831899999999997},
        {215, 24790900, 7.1844999999999999},
        {216, 24746100, 7.1860600000000003},
        {217, 24694100, 7.1878700000000002},
        {218, 24635500, 7.1898999999999997},
        {219, 24570900, 7.1921299999999997},
        {220, 24500900, 7.1945300000000003},
        {221, 24426200, 7.1970799999999997},
        {222, 24347800, 7.1997400000000003},
        {223, 24266300, 7.2024800000000004},
        {224, 24182800, 7.2052699999999996},
        {225, 24098300, 7.2080799999999998},
        {226, 24013600, 7.2108800000000004},
        {227, 23929800, 7.2136199999999997},
        {228, 23847900, 7.2162899999999999},
        {229, 23769000, 7.2188400000000001},
        {230, 23694000, 7.2212500000000004},
        {231, 23623900, 7.22349},
        {232, 23559700, 7.22553},
        {233, 23502200, 7.2273500000000004},
        {234, 23452100, 7.2289199999999996},
        {235, 23410100, 7.2302400000000002},
        {236, 23376900, 7.2312700000000003},
        {237, 23352800, 7.23203},
        {238, 23338200, 7.2324799999999998},
        {239, 23333300, 7.2326300000000003},
        {240, 23333300, 7.2326300000000003},
        {241, 23328400, 7.23278},
        {242, 23313800, 7.2332400000000003},
        {243, 23289600, 7.2339900000000004},
        {244, 23256100, 7.2350199999999996},
        {245, 23213500, 7.2363400000000002},
        {246, 23162300, 7.2379100000000003},
        {247, 23103000, 7.2397200000000002},
        {248, 23036300, 7.2417499999999997},
        {249, 22962600, 7.2439799999999996},
        {250, 22883000, 7.2463699999999998},
        {251, 22798100, 7.2488999999999999},
        {252, 22709000, 7.2515400000000003},
        {253, 22616500, 7.2542499999999999},
        {254, 22521800, 7.2570100000000002},
        {255, 22425900, 7.2597800000000001},
        {256, 22330000, 7.2625299999999999},
        {257, 22235200, 7.2652200000000002},
        {258, 22142700, 7.26783},
        {259, 22053500, 7.2703300000000004},
        {260, 21968900, 7.2726699999999997},
        {261, 21889800, 7.2748499999999998},
        {262, 21817400, 7.2768300000000004},
        {263, 21752600, 7.2785900000000003},
        {264, 21696200, 7.2801200000000001},
        {265, 21648900, 7.28139},
        {266, 21611500, 7.2823900000000004},
        {267, 21584400, 7.2831200000000003},
        {268, 21568000, 7.2835599999999996},
        {269, 21562500, 7.2836999999999996},
        {270, 21562500, 7.2836999999999996},
        {271, 21557000, 7.2838500000000002},
        {272, 21540600, 7.2842900000000004},
        {273, 21513400, 7.2850099999999998},
        {274, 21475600, 7.2860100000000001},
        {275, 21427800, 7.2872700000000004},
        {276, 21370300, 7.28878},
        {277, 21303700, 7.2905199999999999},
        {278, 21228700, 7.2924699999999998},
        {279, 21146100, 7.2946},
        {280, 21056800, 7.2968799999999998},
        {281, 20961700, 7.2992900000000001},
        {282, 20861800, 7.3018000000000001},
        {283, 20758400, 7.3043699999999996},
        {284, 20652500, 7.3069800000000003},
        {285, 20545300, 7.3095999999999997},
        {286, 20438200, 7.3121900000000002},
        {287, 20332400, 7.3147200000000003},
        {288, 20229100, 7.31717},
        {289, 20129700, 7.3194999999999997},
        {290, 20035400, 7.3216900000000003},
        {291, 19947400, 7.3237199999999998},
        {292, 19866800, 7.3255699999999999},
        {293, 19794700, 7.32721},
        {294, 19731900, 7.3286199999999999},
        {295, 19679400, 7.3297999999999996},
        {296, 19637800, 7.33073},
        {297, 19607700, 7.3314000000000004},
        {298, 19589400, 7.3318000000000003},
        {299, 19583300, 7.3319400000000003},
        {300, 19583300, 7.3319400000000003},
        {301, 19577200, 7.3320699999999999},
        {302, 19559000, 7.3324800000000003},
        {303, 19528700, 7.3331499999999998},
        {304, 19486800, 7.3340699999999996},
        {305, 19433700, 7.3352300000000001},
        {306, 19369800, 7.3366199999999999},
        {307, 19296000, 7.3382199999999997},
        {308, 19212800, 7.3399999999999999},
        {309, 19121200, 7.3419499999999998},
        {310, 19022200, 7.3440399999999997},
        {311, 18916900, 7.3462399999999999},
        {312, 18806400, 7.3485199999999997},
        {313, 18691900, 7.3508500000000003},
        {314, 18574800, 7.3532200000000003},
        {315, 18456300, 7.3555799999999998},
        {316, 18338000, 7.3579100000000004},
        {317, 18221200, 7.3601799999999997},
        {318, 18107200, 7.3623700000000003},
        {319, 17997600, 7.3644600000000002},
        {320, 17893600, 7.3664100000000001},
        {321, 17796600, 7.36822},
        {322, 17707800, 7.3698600000000001},
        {323, 17628400, 7.3713100000000003},
        {324, 17559300, 7.3725699999999996},
        {325, 17501500, 7.3736100000000002},
        {326, 17455700, 7.3744300000000003},
        {327, 17422600, 7.3750200000000001},
        {328, 17402500, 7.3753799999999998},
        {329, 17395800, 7.3754999999999997},
        {330, 17395800, 7.3754999999999997},
        {331, 17389100, 7.3756199999999996},
        {332, 17369000, 7.3759699999999997},
        {333, 17335800, 7.3765599999999996},
        {334, 17289700, 7.37737},
        {335, 17231300, 7.3784000000000001},
        {336, 17161100, 7.3796200000000001},
        {337, 17079900, 7.3810200000000004},
        {338, 16988600, 7.3825799999999999},
        {339, 16888000, 7.3842800000000004},
        {340, 16779400, 7.3860999999999999},
        {341, 16663800, 7.3879999999999999},
        {342, 16542600, 7.3899800000000004},
        {343, 16417100, 7.3920000000000003},
        {344, 16288700, 7.3940400000000004},
        {345, 16159000, 7.3960699999999999},
        {346, 16029500, 7.3980600000000001},
        {347, 15901600, 7.40001},
        {348, 15777000, 7.4018800000000002},
        {349, 15657100, 7.4036499999999998},
        {350, 15543500, 7.4053100000000001},
        {351, 15437500, 7.4068399999999999},
        {352, 15340600, 7.4082299999999996},
        {353, 15253800, 7.4094600000000002},
        {354, 15178400, 7.4105100000000004},
        {355, 15115300, 7.4113899999999999},
        {356, 15065400, 7.4120799999999996},
        {357, 15029200, 7.4125800000000002},
        {358, 15007300, 7.4128800000000004},
        {359, 15000000, 7.4129800000000001},
    };

    // Intermediate tMax=1 s ESTER profile: call_id=18, j=0
    /*const std::vector<EmbeddedBucketPoint> points = {
        {0, 27015800, 8.4334100000000003},
        {1, 27015700, 8.4334199999999999},
        {2, 27015100, 8.4335400000000007},
        {3, 27012300, 8.4340700000000002},
        {4, 27004900, 8.4354499999999994},
        {5, 26989900, 8.4382900000000003},
        {6, 26963300, 8.4432799999999997},
        {7, 26921300, 8.4511500000000002},
        {8, 26860000, 8.4625699999999995},
        {9, 26775800, 8.4780999999999995},
        {10, 26666100, 8.4981000000000009},
        {11, 26529200, 8.5226500000000005},
        {12, 26364800, 8.5515100000000004},
        {13, 26173900, 8.5841499999999993},
        {14, 25959000, 8.6196800000000007},
        {15, 25723800, 8.6570099999999996},
        {16, 25473000, 8.6949000000000005},
        {17, 25212000, 8.7320600000000006},
        {18, 24946600, 8.7673000000000005},
        {19, 24682800, 8.7995999999999999},
        {20, 24426300, 8.8282000000000007},
        {21, 24182400, 8.8526399999999992},
        {22, 23956100, 8.8727999999999998},
        {23, 23751600, 8.8887800000000006},
        {24, 23572600, 8.9009599999999995},
        {25, 23422100, 8.9098299999999995},
        {26, 23302700, 8.9159600000000001},
        {27, 23216000, 8.9198699999999995},
        {28, 23163500, 8.9220299999999995},
        {29, 23145900, 8.9227100000000004},
        {30, 23145900, 8.9227100000000004},
        {31, 23134900, 8.9231400000000001},
        {32, 23102000, 8.9243699999999997},
        {33, 23047500, 8.9262599999999992},
        {34, 22972100, 8.9285800000000002},
        {35, 22876400, 8.9309999999999992},
        {36, 22761600, 8.9331300000000002},
        {37, 22628900, 8.9344999999999999},
        {38, 22479900, 8.9346200000000007},
        {39, 22316300, 8.9329599999999996},
        {40, 22140100, 8.9290199999999995},
        {41, 21953500, 8.9223099999999995},
        {42, 21758900, 8.9124599999999994},
        {43, 21558500, 8.8991600000000002},
        {44, 21355100, 8.8822700000000001},
        {45, 21150900, 8.8617799999999995},
        {46, 20948700, 8.8378700000000006},
        {47, 20750800, 8.8109000000000002},
        {48, 20559500, 8.7813800000000004},
        {49, 20377100, 8.7499900000000004},
        {50, 20205700, 8.7175399999999996},
        {51, 20047100, 8.6849299999999996},
        {52, 19903000, 8.6531300000000009},
        {53, 19775100, 8.6231100000000005},
        {54, 19664600, 8.5958100000000002},
        {55, 19572600, 8.5721299999999996},
        {56, 19500200, 8.5528399999999998},
        {57, 19447900, 8.5385899999999992},
        {58, 19416300, 8.5298400000000001},
        {59, 19405800, 8.5268899999999999},
        {60, 19405800, 8.5268899999999999},
        {61, 19393800, 8.5235299999999992},
        {62, 19358100, 8.5134299999999996},
        {63, 19299100, 8.4964600000000008},
        {64, 19217700, 8.47241},
        {65, 19114900, 8.4410600000000002},
        {66, 18992400, 8.4021600000000003},
        {67, 18851600, 8.3554600000000008},
        {68, 18694500, 8.3008199999999999},
        {69, 18523400, 8.2381700000000002},
        {70, 18340300, 8.1676099999999998},
        {71, 18147900, 8.0894100000000009},
        {72, 17948400, 8.0040600000000008},
        {73, 17744500, 7.9122599999999998},
        {74, 17538700, 7.8149600000000001},
        {75, 17333400, 7.7133200000000004},
        {76, 17131100, 7.6086999999999998},
        {77, 16934100, 7.5026200000000003},
        {78, 16744500, 7.3967200000000002},
        {79, 16564500, 7.2927200000000001},
        {80, 16395900, 7.1923599999999999},
        {81, 16240400, 7.0973300000000004},
        {82, 16099600, 7.0092699999999999},
        {83, 15974800, 6.9296899999999999},
        {84, 15867200, 6.8599399999999999},
        {85, 15777900, 6.8012100000000002},
        {86, 15707600, 6.7544899999999997},
        {87, 15656900, 6.7205599999999999},
        {88, 15626300, 6.69998},
        {89, 15616100, 6.6930800000000001},
        {90, 15616100, 6.6930800000000001},
        {91, 15602700, 6.6840299999999999},
        {92, 15562700, 6.6569500000000001},
        {93, 15496800, 6.6120000000000001},
        {94, 15406000, 6.5494599999999998},
        {95, 15291700, 6.4697899999999997},
        {96, 15155700, 6.3735900000000001},
        {97, 14999900, 6.2617399999999996},
        {98, 14826900, 6.1353299999999997},
        {99, 14639100, 5.9957200000000004},
        {100, 14439400, 5.84457},
        {101, 14230400, 5.68377},
        {102, 14015200, 5.5154699999999997},
        {103, 13796600, 5.3419999999999996},
        {104, 13577300, 5.1658099999999996},
        {105, 13360300, 4.9893900000000002},
        {106, 13147800, 4.8152100000000004},
        {107, 12942400, 4.6456400000000002},
        {108, 12746100, 4.4828799999999998},
        {109, 12561000, 4.3289400000000002},
        {110, 12388600, 4.1855599999999997},
        {111, 12230600, 4.0542199999999999},
        {112, 12088300, 3.93615},
        {113, 11962800, 3.8323},
        {114, 11855100, 3.74342},
        {115, 11765900, 3.6700699999999999},
        {116, 11696000, 3.61267},
        {117, 11645700, 3.5714899999999998},
        {118, 11615400, 3.5467200000000001},
        {119, 11605300, 3.5384500000000001},
        {120, 11605300, 3.5384500000000001},
        {121, 11597200, 3.5318700000000001},
        {122, 11573200, 3.5122399999999998},
        {123, 11533600, 3.4799199999999999},
        {124, 11479100, 3.4354800000000001},
        {125, 11410500, 3.3797100000000002},
        {126, 11328900, 3.3136100000000002},
        {127, 11235500, 3.2383199999999999},
        {128, 11131700, 3.1551499999999999},
        {129, 11019200, 3.0655100000000002},
        {130, 10899400, 2.9708600000000001},
        {131, 10774100, 2.8727299999999998},
        {132, 10645000, 2.7726099999999998},
        {133, 10513700, 2.67197},
        {134, 10381900, 2.57219},
        {135, 10251200, 2.4745400000000002},
        {136, 10123100, 2.3801800000000002},
        {137, 9999010, 2.29013},
        {138, 9880220, 2.2052299999999998},
        {139, 9767930, 2.12622},
        {140, 9663200, 2.0536599999999998},
        {141, 9567010, 1.9880199999999999},
        {142, 9480200, 1.92964},
        {143, 9403510, 1.8787499999999999},
        {144, 9337580, 1.83552},
        {145, 9282930, 1.8000700000000001},
        {146, 9239990, 1.7724500000000001},
        {147, 9209090, 1.75271},
        {148, 9190460, 1.7408600000000001},
        {149, 9184230, 1.73691},
        {150, 9184230, 1.73691},
        {151, 9178560, 1.7333099999999999},
        {152, 9161650, 1.72258},
        {153, 9133760, 1.70495},
        {154, 9095330, 1.68082},
        {155, 9046930, 1.6506799999999999},
        {156, 8989310, 1.61517},
        {157, 8923330, 1.575},
        {158, 8849950, 1.53095},
        {159, 8770210, 1.48386},
        {160, 8685230, 1.4345600000000001},
        {161, 8596170, 1.38388},
        {162, 8504190, 1.3326199999999999},
        {163, 8410480, 1.28152},
        {164, 8316190, 1.2312799999999999},
        {165, 8222450, 1.18249},
        {166, 8130350, 1.1356900000000001},
        {167, 8040910, 1.0913299999999999},
        {168, 7955090, 1.04976},
        {169, 7873790, 1.01128},
        {170, 7797820, 0.97612299999999996},
        {171, 7727920, 0.94444899999999998},
        {172, 7664740, 0.91637599999999997},
        {173, 7608840, 0.89198299999999997},
        {174, 7560740, 0.87131700000000001},
        {175, 7520830, 0.85440099999999997},
        {176, 7489440, 0.84124600000000005},
        {177, 7466840, 0.83185399999999998},
        {178, 7453210, 0.82621999999999995},
        {179, 7448650, 0.82434300000000005},
        {180, 7448650, 0.82434300000000005},
        {181, 7444320, 0.82255199999999995},
        {182, 7431390, 0.81722300000000003},
        {183, 7410050, 0.80847800000000003},
        {184, 7380630, 0.79652000000000001},
        {185, 7343580, 0.78161700000000001},
        {186, 7299440, 0.76409800000000005},
        {187, 7248860, 0.74433199999999999},
        {188, 7192560, 0.722723},
        {189, 7131340, 0.69968900000000001},
        {190, 7066050, 0.67564999999999997},
        {191, 6997560, 0.65101699999999996},
        {192, 6926760, 0.62618099999999999},
        {193, 6854570, 0.60150099999999995},
        {194, 6781860, 0.57730400000000004},
        {195, 6709500, 0.55387500000000001},
        {196, 6638330, 0.53145900000000001},
        {197, 6569150, 0.51026099999999996},
        {198, 6502710, 0.49044500000000002},
        {199, 6439690, 0.47214},
        {200, 6380750, 0.45544400000000002},
        {201, 6326460, 0.44042700000000001},
        {202, 6277340, 0.42713699999999999},
        {203, 6233850, 0.41560200000000003},
        {204, 6196380, 0.40583900000000001},
        {205, 6165270, 0.39785500000000001},
        {206, 6140790, 0.39165100000000003},
        {207, 6123150, 0.38722200000000001},
        {208, 6112510, 0.38456699999999999},
        {209, 6108950, 0.383683},
        {210, 6108950, 0.383683},
        {211, 6105540, 0.38283099999999998},
        {212, 6095340, 0.38029800000000002},
        {213, 6078520, 0.37614399999999998},
        {214, 6055310, 0.37046600000000002},
        {215, 6026060, 0.36339399999999999},
        {216, 5991190, 0.35508800000000001},
        {217, 5951200, 0.34572599999999998},
        {218, 5906640, 0.33550099999999999},
        {219, 5858130, 0.32461299999999998},
        {220, 5806320, 0.31326300000000001},
        {221, 5751880, 0.301647},
        {222, 5695520, 0.28994700000000001},
        {223, 5637940, 0.278335},
        {224, 5579840, 0.266961},
        {225, 5521900, 0.25596000000000002},
        {226, 5464800, 0.24544299999999999},
        {227, 5409190, 0.23550499999999999},
        {228, 5355660, 0.22622100000000001},
        {229, 5304800, 0.21765100000000001},
        {230, 5257130, 0.209838},
        {231, 5213140, 0.20281399999999999},
        {232, 5173270, 0.196599},
        {233, 5137910, 0.19120699999999999},
        {234, 5107410, 0.186644},
        {235, 5082060, 0.18291399999999999},
        {236, 5062090, 0.18001400000000001},
        {237, 5047690, 0.17794599999999999},
        {238, 5039000, 0.176705},
        {239, 5036100, 0.176292},
        {240, 5036100, 0.176292},
        {241, 5033340, 0.1759},
        {242, 5025100, 0.174734},
        {243, 5011500, 0.172821},
        {244, 4992730, 0.170207},
        {245, 4969050, 0.16694999999999999},
        {246, 4940790, 0.16312499999999999},
        {247, 4908330, 0.15881200000000001},
        {248, 4872120, 0.15410099999999999},
        {249, 4832620, 0.14908299999999999},
        {250, 4790350, 0.14385000000000001},
        {251, 4745870, 0.138492},
        {252, 4699710, 0.13309199999999999},
        {253, 4652460, 0.12773000000000001},
        {254, 4604680, 0.122474},
        {255, 4556940, 0.11738700000000001},
        {256, 4509760, 0.112521},
        {257, 4463670, 0.10792},
        {258, 4419180, 0.10362},
        {259, 4376790, 0.099648100000000003},
        {260, 4336950, 0.096024700000000004},
        {261, 4300090, 0.092765},
        {262, 4266620, 0.089879399999999998},
        {263, 4236860, 0.087374300000000002},
        {264, 4211150, 0.085253499999999996},
        {265, 4189740, 0.083518700000000001},
        {266, 4172860, 0.082170099999999996},
        {267, 4160670, 0.081207399999999999},
        {268, 4153310, 0.080630099999999996},
        {269, 4150850, 0.080437800000000004},
        {270, 4150850, 0.080437800000000004},
        {271, 4148580, 0.080262},
        {272, 4141800, 0.079738600000000007},
        {273, 4130590, 0.078879900000000003},
        {274, 4115110, 0.077705399999999994},
        {275, 4095550, 0.076241500000000004},
        {276, 4072170, 0.0745199},
        {277, 4045250, 0.072577000000000003},
        {278, 4015140, 0.070451799999999995},
        {279, 3982220, 0.068184800000000004},
        {280, 3946880, 0.065816799999999995},
        {281, 3909570, 0.063387700000000005},
        {282, 3870740, 0.060935200000000002},
        {283, 3830850, 0.058494600000000001},
        {284, 3790380, 0.056097399999999999},
        {285, 3749820, 0.053771800000000002},
        {286, 3709630, 0.0515421},
        {287, 3670280, 0.049428699999999999},
        {288, 3632240, 0.047448299999999999},
        {289, 3595930, 0.045614399999999999},
        {290, 3561760, 0.0439376},
        {291, 3530090, 0.042425999999999998},
        {292, 3501240, 0.041085799999999999},
        {293, 3475540, 0.039920600000000001},
        {294, 3453280, 0.038932899999999999},
        {295, 3434710, 0.038123999999999998},
        {296, 3420050, 0.037494600000000003},
        {297, 3409450, 0.037045099999999997},
        {298, 3403040, 0.036775299999999997},
        {299, 3400900, 0.0366854},
        {300, 3400900, 0.0366854},
        {301, 3399010, 0.036607899999999999},
        {302, 3393350, 0.036377199999999998},
        {303, 3384000, 0.0359984},
        {304, 3371060, 0.0354796},
        {305, 3354690, 0.034831899999999999},
        {306, 3335090, 0.034068899999999999},
        {307, 3312480, 0.033205699999999998},
        {308, 3287160, 0.032259000000000003},
        {309, 3259410, 0.031246199999999998},
        {310, 3229570, 0.030184699999999998},
        {311, 3198010, 0.029092099999999999},
        {312, 3165110, 0.027984700000000001},
        {313, 3131270, 0.026878200000000001},
        {314, 3096920, 0.025786799999999999},
        {315, 3062460, 0.024723499999999999},
        {316, 3028290, 0.023699700000000001},
        {317, 2994820, 0.022725200000000001},
        {318, 2962470, 0.021808299999999999},
        {319, 2931620, 0.020955600000000001},
        {320, 2902640, 0.020172800000000001},
        {321, 2875860, 0.019464200000000001},
        {322, 2851570, 0.018833300000000001},
        {323, 2830050, 0.018282699999999999},
        {324, 2811490, 0.017814400000000001},
        {325, 2796090, 0.017429699999999999},
        {326, 2783980, 0.017129599999999998},
        {327, 2775260, 0.0169149},
        {328, 2770000, 0.016785899999999999},
        {329, 2768240, 0.016742900000000002},
        {330, 2768240, 0.016742900000000002},
        {331, 2766740, 0.016707699999999999},
        {332, 2762270, 0.0166029},
        {333, 2754900, 0.016430500000000001},
        {334, 2744730, 0.016193900000000001},
        {335, 2731930, 0.0158978},
        {336, 2716690, 0.015547699999999999},
        {337, 2699230, 0.0151503},
        {338, 2679820, 0.014712599999999999},
        {339, 2658750, 0.014242400000000001},
        {340, 2636330, 0.013747199999999999},
        {341, 2612910, 0.013235},
        {342, 2588840, 0.0127133},
        {343, 2564460, 0.0121894},
        {344, 2540130, 0.011670399999999999},
        {345, 2516170, 0.011162399999999999},
        {346, 2492910, 0.0106712},
        {347, 2470630, 0.0102019},
        {348, 2449590, 0.0097589100000000008},
        {349, 2429990, 0.0093458699999999992},
        {350, 2412020, 0.0089659200000000005},
        {351, 2395790, 0.0086215600000000003},
        {352, 2381400, 0.0083147199999999994},
        {353, 2368920, 0.0080468699999999994},
        {354, 2358350, 0.0078190599999999992},
        {355, 2349730, 0.0076320199999999998},
        {356, 2343040, 0.00748622},
        {357, 2338280, 0.0073819200000000001},
        {358, 2335430, 0.0073192999999999999},
        {359, 2334480, 0.00729841},
    };*/

    // final converged profile
    /*const std::vector<EmbeddedBucketPoint> points = {
        {0, 25241900, 16.4331},
        {1, 25241900, 16.433},
        {2, 25240800, 16.431899999999999},
        {3, 25236000, 16.4269},
        {4, 25223600, 16.413900000000002},
        {5, 25198000, 16.387},
        {6, 25152800, 16.339700000000001},
        {7, 25081200, 16.264800000000001},
        {8, 24976300, 16.1554},
        {9, 24831600, 16.004899999999999},
        {10, 24641700, 15.8085},
        {11, 24402800, 15.562799999999999},
        {12, 24113000, 15.2669},
        {13, 23772500, 14.922499999999999},
        {14, 23384000, 14.5336},
        {15, 22952600, 14.1068},
        {16, 22485200, 13.650600000000001},
        {17, 21991100, 13.1751},
        {18, 21480400, 12.6913},
        {19, 20964700, 12.2103},
        {20, 20455700, 11.7431},
        {21, 19964900, 11.299899999999999},
        {22, 19503600, 10.8896},
        {23, 19082100, 10.520200000000001},
        {24, 18709600, 10.1981},
        {25, 18394000, 9.9283300000000008},
        {26, 18141700, 9.7148699999999995},
        {27, 17958000, 9.5605200000000004},
        {28, 17846300, 9.4672000000000001},
        {29, 17808800, 9.4359699999999993},
        {30, 17808800, 9.4359699999999993},
        {31, 17801200, 9.4296100000000003},
        {32, 17778400, 9.4106100000000001},
        {33, 17740700, 9.3791399999999996},
        {34, 17688800, 9.3355200000000007},
        {35, 17623100, 9.2801799999999997},
        {36, 17544700, 9.2136999999999993},
        {37, 17454500, 9.1367700000000003},
        {38, 17353800, 9.0502099999999999},
        {39, 17243900, 8.9549800000000008},
        {40, 17126100, 8.8521300000000007},
        {41, 17002000, 8.7428699999999999},
        {42, 16873100, 8.6284500000000008},
        {43, 16741100, 8.5102499999999992},
        {44, 16607500, 8.3896999999999995},
        {45, 16474000, 8.2682599999999997},
        {46, 16342100, 8.1474299999999999},
        {47, 16213300, 8.0286799999999996},
        {48, 16089200, 7.9134599999999997},
        {49, 15971000, 7.8031899999999998},
        {50, 15860100, 7.6991800000000001},
        {51, 15757600, 7.6026699999999998},
        {52, 15664600, 7.5147899999999996},
        {53, 15582000, 7.4365399999999999},
        {54, 15510800, 7.3688200000000004},
        {55, 15451500, 7.3123800000000001},
        {56, 15404800, 7.26783},
        {57, 15371200, 7.2356699999999998},
        {58, 15350800, 7.2162300000000004},
        {59, 15344000, 7.2097300000000004},
        {60, 15344000, 7.2097300000000004},
        {61, 15333000, 7.1991500000000004},
        {62, 15300000, 7.1675599999999999},
        {63, 15245600, 7.1154000000000002},
        {64, 15170600, 7.0433599999999998},
        {65, 15076200, 6.9524499999999998},
        {66, 14963800, 6.8439399999999999},
        {67, 14835100, 6.7193399999999999},
        {68, 14692000, 6.5804299999999998},
        {69, 14536700, 6.4291900000000002},
        {70, 14371300, 6.2677899999999998},
        {71, 14198200, 6.0985199999999997},
        {72, 14019700, 5.9237799999999998},
        {73, 13838200, 5.7460100000000001},
        {74, 13655900, 5.5676100000000002},
        {75, 13475200, 5.39093},
        {76, 13298100, 5.2181899999999999},
        {77, 13126500, 5.0514299999999999},
        {78, 12962300, 4.8925299999999998},
        {79, 12807200, 4.7431200000000002},
        {80, 12662500, 4.6046199999999997},
        {81, 12529700, 4.4782299999999999},
        {82, 12409900, 4.3649100000000001},
        {83, 12304100, 4.26546},
        {84, 12213100, 4.1804699999999997},
        {85, 12137800, 4.1104000000000003},
        {86, 12078700, 4.0556000000000001},
        {87, 12036100, 4.0163000000000002},
        {88, 12010500, 3.9926699999999999},
        {89, 12001900, 3.9847800000000002},
        {90, 12001900, 3.9847800000000002},
        {91, 11987500, 3.9714999999999998},
        {92, 11944400, 3.9319999999999999},
        {93, 11873700, 3.8672800000000001},
        {94, 11776500, 3.7789700000000002},
        {95, 11654900, 3.6692800000000001},
        {96, 11511200, 3.5409099999999998},
        {97, 11347900, 3.39696},
        {98, 11168000, 3.2407699999999999},
        {99, 10974700, 3.07585},
        {100, 10771000, 2.9057200000000001},
        {101, 10560300, 2.7337400000000001},
        {102, 10345700, 2.5630500000000001},
        {103, 10129900, 2.3964400000000001},
        {104, 9915860, 2.2362899999999999},
        {105, 9705940, 2.0845199999999999},
        {106, 9502430, 1.94259},
        {107, 9307360, 1.81152},
        {108, 9122500, 1.6919299999999999},
        {109, 8949430, 1.58412},
        {110, 8789490, 1.4881},
        {111, 8643810, 1.4036999999999999},
        {112, 8513360, 1.3306},
        {113, 8398940, 1.26841},
        {114, 8301170, 1.2166999999999999},
        {115, 8220560, 1.1750499999999999},
        {116, 8157490, 1.1430899999999999},
        {117, 8112250, 1.1205000000000001},
        {118, 8085040, 1.10704},
        {119, 8075950, 1.1025700000000001},
        {120, 8075950, 1.1025700000000001},
        {121, 8059000, 1.09426},
        {122, 8008600, 1.06978},
        {123, 7926040, 1.03043},
        {124, 7813430, 0.97823099999999996},
        {125, 7673580, 0.91576599999999997},
        {126, 7509850, 0.84592299999999998},
        {127, 7325980, 0.77167699999999995},
        {128, 7125990, 0.69586099999999995},
        {129, 6913920, 0.62099099999999996},
        {130, 6693790, 0.54913299999999998},
        {131, 6469380, 0.48185499999999998},
        {132, 6244180, 0.420209},
        {133, 6021260, 0.364782},
        {134, 5803270, 0.31575999999999999},
        {135, 5592520, 0.27301599999999998},
        {136, 5390970, 0.23619899999999999},
        {137, 5200310, 0.20482},
        {138, 5021910, 0.17832200000000001},
        {139, 4856830, 0.15612899999999999},
        {140, 4705750, 0.13769400000000001},
        {141, 4569180, 0.122507},
        {142, 4447560, 0.11011},
        {143, 4341330, 0.100109},
        {144, 4250850, 0.092170500000000002},
        {145, 4176420, 0.086023100000000005},
        {146, 4118300, 0.081452999999999998},
        {147, 4076660, 0.078299499999999994},
        {148, 4051630, 0.076451099999999994},
        {149, 4043280, 0.075842199999999999},
        {150, 4043280, 0.075842199999999999},
        {151, 4035890, 0.075306499999999998},
        {152, 4013870, 0.073727899999999999},
        {153, 3977650, 0.071188299999999996},
        {154, 3927920, 0.067815500000000001},
        {155, 3865600, 0.063771099999999997},
        {156, 3791740, 0.059236900000000002},
        {157, 3707490, 0.054399299999999998},
        {158, 3614060, 0.049436500000000001},
        {159, 3512660, 0.044506499999999997},
        {160, 3404530, 0.039740200000000003},
        {161, 3290960, 0.0352371},
        {162, 3173390, 0.031064899999999999},
        {163, 3053480, 0.027262499999999999},
        {164, 2932770, 0.0238471},
        {165, 2812700, 0.020818900000000001},
        {166, 2694650, 0.018164799999999998},
        {167, 2579980, 0.015861799999999999},
        {168, 2469990, 0.013882},
        {169, 2366200, 0.0121932},
        {170, 2270170, 0.010762799999999999},
        {171, 2183040, 0.0095620600000000007},
        {172, 2105510, 0.0085652000000000002},
        {173, 2037940, 0.00774955},
        {174, 1980470, 0.0070950900000000001},
        {175, 1933240, 0.0065842299999999999},
        {176, 1896420, 0.0062021400000000001},
        {177, 1870080, 0.00593736},
        {178, 1854260, 0.0057817199999999997},
        {179, 1848990, 0.0057303900000000001},
        {180, 1848990, 0.0057303900000000001},
        {181, 1845460, 0.0056961700000000004},
        {182, 1834940, 0.0055948999999999999},
        {183, 1817640, 0.0054306299999999997},
        {184, 1793870, 0.0052097799999999998},
        {185, 1764110, 0.0049407100000000001},
        {186, 1728900, 0.0046331699999999998},
        {187, 1688890, 0.0042976400000000001},
        {188, 1644800, 0.0039448599999999997},
        {189, 1597350, 0.0035851699999999999},
        {190, 1547260, 0.0032281100000000002},
        {191, 1495230, 0.00288193},
        {192, 1441910, 0.00255331},
        {193, 1387970, 0.00224727},
        {194, 1334040, 0.0019671200000000002},
        {195, 1280750, 0.00171463},
        {196, 1228750, 0.0014902699999999999},
        {197, 1178610, 0.00129349},
        {198, 1130860, 0.0011230400000000001},
        {199, 1085960, 0.00097712199999999997},
        {200, 1044330, 0.00085362999999999997},
        {201, 1006340, 0.000750316},
        {202, 972329, 0.00066494500000000001},
        {203, 942595, 0.00059540099999999996},
        {204, 917319, 0.00053981400000000001},
        {205, 896579, 0.00049658799999999996},
        {206, 880408, 0.000464389},
        {207, 868836, 0.00044215600000000003},
        {208, 861885, 0.000429123},
        {209, 859567, 0.00042482899999999999},
        {210, 859567, 0.00042482899999999999},
        {211, 858146, 0.00042221000000000001},
        {212, 853908, 0.00041446000000000001},
        {213, 846928, 0.000401889},
        {214, 837331, 0.00038499600000000003},
        {215, 825284, 0.000364426},
        {216, 810993, 0.00034093799999999999},
        {217, 794694, 0.00031535100000000001},
        {218, 776646, 0.00028850000000000002},
        {219, 757115, 0.00026119099999999997},
        {220, 736383, 0.00023415600000000001},
        {221, 714738, 0.00020802400000000001},
        {222, 692475, 0.00018329699999999999},
        {223, 669873, 0.00016035600000000001},
        {224, 647192, 0.00013945099999999999},
        {225, 624702, 0.00012071},
        {226, 602653, 0.00010416099999999999},
        {227, 581252, 8.9755100000000004e-05},
        {228, 560715, 7.7374100000000002e-05},
        {229, 541240, 6.6861599999999997e-05},
        {230, 522988, 5.8041500000000003e-05},
        {231, 506116, 5.0727199999999998e-05},
        {232, 490798, 4.4733800000000002e-05},
        {233, 477184, 3.9891700000000003e-05},
        {234, 465402, 3.6051099999999998e-05},
        {235, 455576, 3.3083000000000002e-05},
        {236, 447816, 3.0882e-05},
        {237, 442212, 2.9366999999999999e-05},
        {238, 438824, 2.84806e-05},
        {239, 437691, 2.81888e-05},
        {240, 437691, 2.81888e-05},
        {241, 437023, 2.8018099999999999e-05},
        {242, 435028, 2.7512999999999999e-05},
        {243, 431731, 2.6693899999999999e-05},
        {244, 427170, 2.5593500000000002e-05},
        {245, 421402, 2.4253999999999999e-05},
        {246, 414496, 2.27248e-05},
        {247, 406537, 2.1058599999999999e-05},
        {248, 397623, 1.9309000000000001e-05},
        {249, 387870, 1.7526599999999999e-05},
        {250, 377413, 1.5757299999999999e-05},
        {251, 366401, 1.4040600000000001e-05},
        {252, 354976, 1.24093e-05},
        {253, 343274, 1.0889e-05},
        {254, 331435, 9.4966699999999997e-06},
        {255, 319593, 8.2422599999999996e-06},
        {256, 307865, 7.1295199999999998e-06},
        {257, 296342, 6.1569700000000001e-06},
        {258, 285113, 5.3185300000000004e-06},
        {259, 274272, 4.6048299999999999e-06},
        {260, 263896, 4.0049499999999997e-06},
        {261, 254064, 3.50706e-06},
        {262, 244899, 3.09879e-06},
        {263, 236524, 2.7687899999999998e-06},
        {264, 229070, 2.5069800000000002e-06},
        {265, 222697, 2.3044700000000001e-06},
        {266, 217558, 2.1540999999999998e-06},
        {267, 213783, 2.05045e-06},
        {268, 211475, 1.98974e-06},
        {269, 210699, 1.9697500000000001e-06},
        {270, 210699, 1.9697500000000001e-06},
        {271, 210336, 1.9604899999999998e-06},
        {272, 209249, 1.933e-06},
        {273, 207442, 1.88823e-06},
        {274, 204923, 1.8276300000000001e-06},
        {275, 201704, 1.7531100000000001e-06},
        {276, 197802, 1.66692e-06},
        {277, 193243, 1.5714399999999999e-06},
        {278, 188066, 1.4691099999999999e-06},
        {279, 182321, 1.36226e-06},
        {280, 176076, 1.2530999999999999e-06},
        {281, 169427, 1.1434900000000001e-06},
        {282, 162508, 1.0349799999999999e-06},
        {283, 155438, 9.2920499999999996e-07},
        {284, 148362, 8.2753800000000001e-07},
        {285, 141420, 7.3123600000000002e-07},
        {286, 134688, 6.4176100000000002e-07},
        {287, 128230, 5.60175e-07},
        {288, 122101, 4.8713400000000005e-07},
        {289, 116332, 4.2296399999999997e-07},
        {290, 110965, 3.6753899999999999e-07},
        {291, 106045, 3.20437e-07},
        {292, 101613, 2.81077e-07},
        {293, 97712.300000000003, 2.4876400000000002e-07},
        {294, 94382.5, 2.2278599999999999e-07},
        {295, 91647.800000000003, 2.0249799999999999e-07},
        {296, 89516, 1.8735099999999999e-07},
        {297, 87990.100000000006, 1.7688500000000001e-07},
        {298, 87073, 1.7075099999999999e-07},
        {299, 86767.100000000006, 1.68731e-07},
        {300, 86767.100000000006, 1.68731e-07},
        {301, 86628.5, 1.6782e-07},
        {302, 86215, 1.6511900000000001e-07},
        {303, 85533.199999999997, 1.60717e-07},
        {304, 84593.600000000006, 1.5475999999999999e-07},
        {305, 83410.699999999997, 1.4744e-07},
        {306, 82002.600000000006, 1.38989e-07},
        {307, 80389.800000000003, 1.2966499999999999e-07},
        {308, 78595.5, 1.1974000000000001e-07},
        {309, 76644.199999999997, 1.0948700000000001e-07},
        {310, 74561.100000000006, 9.9168200000000003e-08},
        {311, 72371.699999999997, 8.9024200000000001e-08},
        {312, 70100.899999999994, 7.9262100000000004e-08},
        {313, 67771.399999999994, 7.0053600000000005e-08},
        {314, 65401.300000000003, 6.1530299999999999e-08},
        {315, 63009.699999999997, 5.3776499999999998e-08},
        {316, 60612.199999999997, 4.6838300000000002e-08},
        {317, 58215.300000000003, 4.07323e-08},
        {318, 55827.599999999999, 3.5439399999999999e-08},
        {319, 53453.800000000003, 3.0921600000000002e-08},
        {320, 51084.900000000001, 2.7135300000000001e-08},
        {321, 48739.099999999999, 2.4009099999999999e-08},
        {322, 46451, 2.1467099999999998e-08},
        {323, 44268.099999999999, 1.94308e-08},
        {324, 42275.599999999999, 1.7801799999999999e-08},
        {325, 40564.800000000003, 1.6497900000000001e-08},
        {326, 39196.099999999999, 1.5480699999999998e-08},
        {327, 38208.099999999999, 1.47453e-08},
        {328, 37616.199999999997, 1.4299e-08},
        {329, 37419.300000000003, 1.41495e-08},
        {330, 37419.300000000003, 1.41495e-08},
        {331, 37342.800000000003, 1.4091100000000001e-08},
        {332, 37114.699999999997, 1.39167e-08},
        {333, 36739.900000000001, 1.3628200000000001e-08},
        {334, 36226.099999999999, 1.3229400000000001e-08},
        {335, 35583.599999999999, 1.2726399999999999e-08},
        {336, 34824, 1.21285e-08},
        {337, 33961, 1.1447900000000001e-08},
        {338, 33011, 1.07e-08},
        {339, 31990.599999999999, 9.9018699999999998e-09},
        {340, 30916.5, 9.07269e-09},
        {341, 29804.700000000001, 8.23191e-09},
        {342, 28671.200000000001, 7.3985099999999999e-09},
        {343, 27529.200000000001, 6.5902799999999997e-09},
        {344, 26388.400000000001, 5.8232200000000001e-09},
        {345, 25257.900000000001, 5.1096099999999997e-09},
        {346, 24153.5, 4.4565899999999998e-09},
        {347, 23091.200000000001, 3.8683799999999998e-09},
        {348, 22082.900000000001, 3.3473700000000002e-09},
        {349, 21136.700000000001, 2.89381e-09},
        {350, 20260.099999999999, 2.5053600000000001e-09},
        {351, 19462.5, 2.17798e-09},
        {352, 18752.799999999999, 1.9061699999999998e-09},
        {353, 18136.799999999999, 1.68446e-09},
        {354, 17616.599999999999, 1.50743e-09},
        {355, 17192.700000000001, 1.36997e-09},
        {356, 16865.200000000001, 1.2676300000000001e-09},
        {357, 16632.799999999999, 1.1969799999999999e-09},
        {358, 16494.299999999999, 1.1555100000000001e-09},
        {359, 16448.299999999999, 1.1418400000000001e-09},
};*/

    if (zones < 1 || zones > static_cast<int>(points.size())) {
        std::cerr
            << "[embedded profile] zones must be in [1,"
            << points.size() << "]\n";
        std::exit(EXIT_FAILURE);
    }

    std::cerr
        << "[embedded profile] start tMax=" << tMax
        << " zones=" << zones
        << " using initial call_id=1 T-rho profile\n"
        << std::flush;

    // Deliberately keep one fixed composition for every zone so that this
    // test isolates the GridSolver response to the heterogeneous T-rho profile.
    const auto comp = build_embedded_bucket_composition();

    GridfireRuntime probe_runtime(comp);

    std::vector<gridfire::NetIn> netIns;
    netIns.reserve(zones);

    for (int i = 0; i < zones; ++i) {
        netIns.push_back(
            make_embedded_bucket_netin(points[i], comp, tMax)
        );
    }

    std::cerr
        << "[embedded profile] before reset_context\n"
        << std::flush;

    probe_runtime.reset_context();

    std::cerr
        << "[embedded profile] before GridSolver evaluate n="
        << netIns.size() << "\n"
        << std::flush;

    const std::vector<gridfire::NetOut> netOuts =
        probe_runtime.gridSolver->evaluate(
            *probe_runtime.solverCtx,
            netIns
        );

    long long total_steps = 0;

    int zero_step_zones = 0;
    int first_zero_step_i = -1;

    int max_steps = -1;
    int max_steps_i = -1;

    for (std::size_t k = 0; k < netOuts.size(); ++k) {

        const int steps = netOuts[k].num_steps;

        if (steps > 0) {
            total_steps += steps;
        }
        else {
            ++zero_step_zones;

            if (first_zero_step_i < 0) {
                first_zero_step_i = points[k].i;
            }
        }

        if (steps > max_steps) {
            max_steps = steps;
            max_steps_i = points[k].i;
        }
    }

    std::cerr
        << "[embedded profile] SUCCESS"
        << " netOuts.size()=" << netOuts.size()
        << " solver_workspaces.size()="
        << probe_runtime.solverCtx->solver_workspaces.size()
        << " total_steps=" << total_steps
        << " max_steps=" << max_steps
        << " max_steps_i=" << max_steps_i
        << " zero_step_zones=" << zero_step_zones
        << " first_zero_step_i=" << first_zero_step_i
        << "\n"
        << std::flush;
}


std::vector<GridfireDiag> run_gridfire_profile(GridfireRuntime& runtime,
                                               const composition_map& comp,
                                               const matrix& T,
                                               const matrix& rho,
                                               double tMax)
{
    const int nr = T.nrows();
    const int nth = T.ncols();

    const char* ester_species[] = {
        "H1",
        "He3",
        "He4",
        "C12",
        "N14",
        "O16",
        "Ne20",
        "Mg24"
    };

    const char* gridfire_species[] = {
        "H-1",
        "He-3",
        "He-4",
        "C-12",
        "N-14",
        "O-16",
        "Ne-20",
        "Mg-24"
    };

    constexpr int n_diag_species = 8;


    //const int ncell = nr * nth;
    //const int ncell = std::min(nr, 180) * nth;

    //std::cerr << "[run_gridfire_profile] entered\n" << std::flush;
    //std::cerr << "[run_gridfire_profile] nr=" << nr << " nth=" << nth << " ncell=" << ncell << "\n" << std::flush;

    //std::vector<GridfireDiag> diags(ncell);

    const int i_start = 0;
//    const int i_end = 180;

//    const int i_start = 180; //0;
    const int i_end = nr; //std::min(nr, 180);
    const int ncell = (i_end - i_start) * nth;

    std::vector<GridfireDiag> diags(nr * nth);
    std::vector<gridfire::NetIn> netIns;
    netIns.reserve(ncell);

    std::cerr << "[run_gridfire_profile] entered\n" << std::flush;
    std::cerr << "[run_gridfire_profile] nr=" << nr << " nth=" << nth << " ncell=" << ncell << "\n" << std::flush;

    //for (int i = 0; i < nr; ++i) {
    for (int i = i_start; i < i_end; ++i) { 
        for (int j = 0; j < nth; ++j) {
            netIns.push_back(make_gridfire_netin(comp, T, rho, i, j, tMax));
        }
    }

    try {


        for (int k = 0; k < ncell; ++k) {


            //const int i = k / nth;
            //const int j = k % nth;

            const int i = i_start + k / nth;
            const int j = k % nth;

            if (!std::isfinite(T(i,j)) || !std::isfinite(rho(i,j)) ||
                T(i,j) <= 0.0 || rho(i,j) <= 0.0) {
                std::cerr << "[bad input before evaluate] k=" << k
                        << " i=" << i << " j=" << j
                        << " T=" << T(i,j)
                        << " rho=" << rho(i,j)
                        << "\n" << std::flush;
            }
        }        


        //runtime.reset_context();
        
        // EMB speedup: one GridSolver call for the whole profile, not one solver construction per cell.
        //std::cerr << "[run_gridfire_profile] before reset_context\n" << std::flush;
        runtime.reset_context();
        //std::cerr << "[run_gridfire_profile] after reset_context\n" << std::flush;

        //std::cerr << "[run_gridfire_profile] before evaluate\n" << std::flush;
        //const std::vector<gridfire::NetOut> netOuts =
        //    runtime.gridSolver->evaluate(*runtime.solverCtx, netIns); //replacing this with below

        //std::cerr << "[run_gridfire_profile] before evaluate\n" << std::flush;
        //const std::vector<gridfire::NetOut> netOuts =
        //    runtime.gridSolver->evaluate(*runtime.solverCtx, netIns);

        double Tmin = std::numeric_limits<double>::infinity();
        double Tmax = -std::numeric_limits<double>::infinity();
        double rhomin = std::numeric_limits<double>::infinity();
        double rhomax = -std::numeric_limits<double>::infinity();

        for (int k = 0; k < ncell; ++k) {
            //const int i = k / nth;
            //const int j = k % nth;

            const int i = i_start + k / nth;
            const int j = k % nth;

            Tmin = std::min(Tmin, T(i,j));
            Tmax = std::max(Tmax, T(i,j));
            rhomin = std::min(rhomin, rho(i,j));
            rhomax = std::max(rhomax, rho(i,j));
        }

            std::cerr << "[profile range before evaluate] "
          << "Tmin=" << Tmin << " Tmax=" << Tmax
          << " rhomin=" << rhomin << " rhomax=" << rhomax
          << "\n" << std::flush;

        std::cerr << "[run_gridfire_profile] before evaluate\n" << std::flush;
        const auto evaluate_t0 = std::chrono::steady_clock::now();
        const std::vector<gridfire::NetOut> netOuts =
            runtime.gridSolver->evaluate(*runtime.solverCtx, netIns);
        const auto evaluate_t1 = std::chrono::steady_clock::now();
        const double evaluate_s =
            std::chrono::duration<double>(evaluate_t1 - evaluate_t0).count();
        std::cerr << "[run_gridfire_profile] after evaluate netOuts.size()="
                << netOuts.size()
                << " evaluate_s=" << evaluate_s
                << "\n" << std::flush;

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

        //std::cerr << "[run_gridfire_profile] after evaluate, netOuts.size()="
        //        << netOuts.size() << "\n" << std::flush;

        //std::cerr << "[post] solver_workspaces.size()="
        //        << runtime.solverCtx->solver_workspaces.size()
        //        << "\n" << std::flush;

        const int nout = std::min<int>(static_cast<int>(netOuts.size()), ncell);
        const auto postprocess_t0 = std::chrono::steady_clock::now();
        double rhs_lookup_s = 0.0;

        for (int k = 0; k < nout; ++k) {
            const int i = i_start + k / nth;
            const int j = k % nth;
            const int global_k = i * nth + j;

            GridfireDiag& diag = diags[global_k];
            bool rhs_available = false;

            //std::cerr << "[post] k=" << k << "\n" << std::flush;

            const auto rhs_t0 = std::chrono::steady_clock::now();

            auto* zone_ctx = dynamic_cast<gridfire::solver::PointSolverContext*>(
                runtime.solverCtx->solver_workspaces[k].get()
            );

            if (zone_ctx && zone_ctx->engine_ctx) {

                // NOTE: Emily has flagged that getMostRecentRHSCalculation() may
                // refer to a CVODE trial state rather than the final accepted
                // state. Keep this path unchanged for the present baseline; once
                // GridFire exposes an accepted-state epsilon, log both side by side.
                auto rhs_calc =
                    runtime.engine->getMostRecentRHSCalculation(*zone_ctx->engine_ctx);

                rhs_available = rhs_calc.has_value();

                if (rhs_available) {
                    diag.eps_inst = rhs_calc->nuclearEnergyGenerationRate;
                }
            }

            const auto rhs_t1 = std::chrono::steady_clock::now();
            rhs_lookup_s +=
                std::chrono::duration<double>(rhs_t1 - rhs_t0).count();

            const auto& in = netIns[k];
            const auto& out = netOuts[k];

            diag.eps_avg = out.energy / tMax; // average energy generation over tMax

            // Log the exact composition GridFire received, including the
            // temporary Mg-24 remainder bucket used by build_test_composition().
            diag.Xin_H1   = in.composition.getMassFraction("H-1");
            diag.Xin_He3  = in.composition.getMassFraction("He-3");
            diag.Xin_He4  = in.composition.getMassFraction("He-4");
            diag.Xin_C12  = in.composition.getMassFraction("C-12");
            diag.Xin_N14  = in.composition.getMassFraction("N-14");
            diag.Xin_O16  = in.composition.getMassFraction("O-16");
            diag.Xin_Ne20 = in.composition.getMassFraction("Ne-20");
            diag.Xin_Mg24 = in.composition.getMassFraction("Mg-24");

            // Tabulate the output compositions to track their changes.

            double* Xcomp_fields[] = {
                &diag.Xcomp_H1,
                &diag.Xcomp_He3,
                &diag.Xcomp_He4,
                &diag.Xcomp_C12,
                &diag.Xcomp_N14,
                &diag.Xcomp_O16,
                &diag.Xcomp_Ne20,
                &diag.Xcomp_Mg24
            };

            double* Xin_fields[] = {
                &diag.Xin_H1,
                &diag.Xin_He3,
                &diag.Xin_He4,
                &diag.Xin_C12,
                &diag.Xin_N14,
                &diag.Xin_O16,
                &diag.Xin_Ne20,
                &diag.Xin_Mg24
            };

            double* Xout_fields[] = {
                &diag.Xout_H1,
                &diag.Xout_He3,
                &diag.Xout_He4,
                &diag.Xout_C12,
                &diag.Xout_N14,
                &diag.Xout_O16,
                &diag.Xout_Ne20,
                &diag.Xout_Mg24
            };

            for (int s = 0; s < n_diag_species; ++s) {

                // Raw ESTER composition before GridFire input mapping/bucketing.
                *Xcomp_fields[s] = comp[ester_species[s]](i,j);

                // Exact composition actually supplied to GridFire.
                *Xin_fields[s] =
                    safe_gridfire_mass_fraction(
                        netIns[k].composition,
                        gridfire_species[s]
                    );

                // Final composition returned by GridFire.
                // Failed/incomplete NetOut compositions remain NaN.
                *Xout_fields[s] =
                    safe_gridfire_mass_fraction(
                        out.composition,
                        gridfire_species[s]
                    );
            }

            // Preserve GridFire's raw result before the existing validity
            // guard below replaces failed-zone diagnostics with NaNs.
            diag.raw_num_steps = out.num_steps;
            diag.raw_netout_energy = out.energy;
            diag.raw_netout_deps_dT = out.dEps_dT;
            diag.raw_netout_deps_drho = out.dEps_dRho;
            diag.rhs_available = rhs_available;

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
            //diag.eps_inst = eps_inst_cache[k];




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

        const auto postprocess_t1 = std::chrono::steady_clock::now();
        const double postprocess_s =
            std::chrono::duration<double>(postprocess_t1 - postprocess_t0).count();

        // Store the per-call timings in each successful/attempted zone row so
        // the CSV can be grouped by call_id without needing a second file.
        for (int k = 0; k < nout; ++k) {
            const int i = i_start + k / nth;
            const int j = k % nth;
            GridfireDiag& diag = diags[i * nth + j];
            diag.evaluate_s = evaluate_s;
            diag.rhs_lookup_s = rhs_lookup_s;
            diag.postprocess_s = postprocess_s;
        }

        std::cerr
            << "[GridFire timing] evaluate_s=" << evaluate_s
            << " rhs_lookup_s=" << rhs_lookup_s
            << " postprocess_s=" << postprocess_s
            << " total_profile_s=" << (evaluate_s + postprocess_s)
            << "\n";


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
    //fname << "gridfire_vs_simple_profile_tMax_" << std::scientific << std::setprecision(0) << tMax << ".csv";
    //fname << "heaptrack_tests_tMax_" << std::scientific << std::setprecision(0) << tMax << ".csv";
    fname << "upgraded_pip_tests_tMax_" << std::scientific << std::setprecision(0) << tMax << "w_comp_diags.csv";

    //fname << "simple_cn_v4_independent_run" << std::scientific << std::setprecision(0) << tMax << ".csv";
    //fname << "simple_on_v4_independent_run" << std::scientific << std::setprecision(0) << tMax << ".csv";


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
        "raw_num_steps,rhs_available,raw_netout_energy,"
        "raw_netout_deps_drho,raw_netout_deps_dT,"
        "eps_gridfire,dlneps_lnrho_gridfire,dlneps_lnT_gridfire,"
        "deps_drho_gridfire,deps_dT_gridfire,"
        "eps_to_ester,dlneps_lnrho_to_ester,dlneps_lnT_to_ester,"
        "eps_simple_CN,dlneps_lnrho_simple_CN,dlneps_lnT_simple_CN,"
        "eps_simple_ON,dlneps_lnrho_simple_ON,dlneps_lnT_simple_ON,"
        "eps_avg_gridfire,"
        "Xcomp_H1,Xcomp_He3,Xcomp_He4,Xcomp_C12,Xcomp_N14,Xcomp_O16,Xcomp_Ne20,Xcomp_Mg24,"
        "Xin_H1,Xin_He3,Xin_He4,Xin_C12,Xin_N14,Xin_O16,Xin_Ne20,Xin_Mg24,"
        "Xout_H1,Xout_He3,Xout_He4,Xout_C12,Xout_N14,Xout_O16,Xout_Ne20,Xout_Mg24,"
        "evaluate_s,rhs_lookup_s,postprocess_s\n";
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

    const double tMax = 1e7; // seconds; same as v2 for direct comparison.

    // One-off embedded reproduction of the standalone GridFire bucket probe.
    // This tests whether the same post-GridSolver RHS retrieval works inside
    // the ESTER executable/linking environment before the live ESTER profile runs.
    //run_embedded_bucket_probe_once(tMax);

    // ----- TEMPORARY TEST, IT WILL STOP THE CODE ON FIRST CALL ----
    const double probe_tMax = 1.0e2;  // embedded isolation test only

    //run_embedded_bucket_probe_once(probe_tMax);
    //run_embedded_full_profile_probe_once(probe_tMax,360); 
    
    // second number says where to slice the profile
    // this is temporary, because if 180 is done, then only half of the netin would be filled
    // this is just a memory consumption test
    
    //std::cerr
    //    << "[embedded isolation] completed tMax="
    //    << probe_tMax
    //    << "; terminating before live ESTER GridFire profile\n";

    //std::exit(EXIT_SUCCESS);
    // ----- END OF TEMPROARY TEST -----




    const std::string filename = gridfire_log_filename(tMax);
    write_gridfire_header_once(filename);

    //GridfireRuntime& runtime = gridfire_runtime(comp, T, rho, tMax); // really would need to be moved outside of this nuc function yk, make separate "intilisation script" for external modules
    
    std::cerr << "[nuc_gridfire] before gridfire_runtime\n" << std::flush;
    GridfireRuntime& runtime = gridfire_runtime(comp, T, rho, tMax);
    std::cerr << "[nuc_gridfire] after gridfire_runtime\n" << std::flush;    
        const std::vector<GridfireDiag> diags = run_gridfire_profile(runtime, comp, T, rho, tMax);

    std::cout << "debug print statement 1\n" << std::endl;

    // nuc logger was in this line before test, now its going after - MG 13/08/26

    
    /*
    
    // ------------------------------------------------------------------
    // TEST: fixed nuclear-source cutoff at an existing domain boundary.
    //
    // Keep the logged simple_CN profile unchanged by applying this only
    // after the CSV output has been written.
    // ------------------------------------------------------------------
    constexpr int zero_tail_start = 150;

    const int nr = T.nrows();
    const int nth = T.ncols();

    for (int i = zero_tail_start; i < nr; ++i) {
        for (int j = 0; j < nth; ++j) {
            simple_cn.pp(i,j) = 0.0;
            simple_cn.cno(i,j) = 0.0;
            simple_cn.eps(i,j) = 0.0;
            simple_cn.dlneps_lnrho(i,j) = 0.0;
            simple_cn.dlneps_lnT(i,j) = 0.0;
        }
    }

    std::cerr
        << "[simple zero-tail test] zeroed indices "
        << zero_tail_start << ".." << nr - 1
        << " for all " << nth << " angular columns\n";

    // Keep ESTER stable for now: log GridFire but return simple_CN values.
    nuc.pp = simple_cn.pp;
    nuc.cno = simple_cn.cno;
    nuc.eps = simple_cn.eps;
    nuc.dlneps_lnrho = simple_cn.dlneps_lnrho;
    nuc.dlneps_lnT = simple_cn.dlneps_lnT;

    */


    const std::string nuc_output = "gridfire"; 
    //const std::string nuc_output = "simple_cn"; 
    //const std::string nuc_output = "simple_on"; 


    if (nuc_output=="gridfire") {

    // ------------------------------------------------------------------
    // TEST: feed the GridFire nuclear profile directly into ESTER.
    //
    // Preserve all successful GridFire points inward of the first failed
    // radial point.  From the first failed radial point outward, set the
    // nuclear source and its logarithmic derivatives to zero.
    //
    // This is the GridFire equivalent of the simple_CN zero-tail test.
    // ------------------------------------------------------------------

    const int nr = T.nrows();
    const int nth = T.ncols();

    int first_failed_i = nr;

    // ------------------------------------------------------------
    // Find the first radial index at which GridFire did not produce
    // a usable result in at least one angular column.
    // ------------------------------------------------------------
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nth; ++j) {

            const int k = i * nth + j;
            const GridfireDiag& diag = diags[k];

            const bool valid =
                diag.rhs_available &&
                std::isfinite(diag.eps_inst) &&
                std::isfinite(diag.dlneps_lnrho) &&
                std::isfinite(diag.dlneps_lnT);

            if (!valid) {
                first_failed_i = i;
                break;
            }
        }

        if (first_failed_i != nr) {
            break;
        }
    }


    // ------------------------------------------------------------
    // Allocate ESTER nuclear-output matrices.
    //
    // GridFire values are copied element-by-element, so unlike
    // nuc_simple() the destination matrices must be explicitly
    // dimensioned before writing into them.
    //
    // They are initialized to zero, which also automatically gives
    // the failed outer tail the desired zero values.
    // ------------------------------------------------------------
    nuc.pp = zeros(nr, nth);
    nuc.cno = zeros(nr, nth);
    nuc.eps = zeros(nr, nth);
    nuc.dlneps_lnrho = zeros(nr, nth);
    nuc.dlneps_lnT = zeros(nr, nth);

    std::cout
        << "[nuc matrix dimensions] "
        << "eps=" << nuc.eps.nrows() << "x" << nuc.eps.ncols()
        << " dlneps_lnrho=" << nuc.dlneps_lnrho.nrows()
        << "x" << nuc.dlneps_lnrho.ncols()
        << " dlneps_lnT=" << nuc.dlneps_lnT.nrows()
        << "x" << nuc.dlneps_lnT.ncols()
        << std::endl;


    // ------------------------------------------------------------
    // Successful GridFire region.
    //
    // Everything from first_failed_i outward remains zero because
    // the matrices were initialized with zeros above.
    // ------------------------------------------------------------
    for (int i = 0; i < first_failed_i; ++i) {
        for (int j = 0; j < nth; ++j) {

            const int k = i * nth + j;
            const GridfireDiag& diag = diags[k];

            nuc.eps(i,j) = diag.eps_inst;
            nuc.dlneps_lnrho(i,j) = diag.dlneps_lnrho;
            nuc.dlneps_lnT(i,j) = diag.dlneps_lnT;
        }
    }


    // ------------------------------------------------------------
    // Report the dynamically selected cutoff.
    // ------------------------------------------------------------
    if (first_failed_i < nr) {

        std::cerr
            << "[GridFire zero-tail] first failed radial index = "
            << first_failed_i
            << "; using GridFire for 0.."
            << first_failed_i - 1
            << " and zeroing "
            << first_failed_i << ".." << nr - 1
            << "\n";

    } else {

        std::cerr
            << "[GridFire zero-tail] all "
            << nr
            << " radial points succeeded; no cutoff applied\n";
    }

    std::cout << "before log flush end of gridfire nuc run \n" << std::endl;



    } else if (nuc_output=="simple_cn") {

        nuc.pp = simple_cn.pp;
        nuc.cno = simple_cn.cno;
        nuc.eps = simple_cn.eps;
        nuc.dlneps_lnrho = simple_cn.dlneps_lnrho;
        nuc.dlneps_lnT = simple_cn.dlneps_lnT;

    } else if (nuc_output=="simple_on") {

        nuc.pp = simple_on.pp;
        nuc.cno = simple_on.cno;
        nuc.eps = simple_on.eps;
        nuc.dlneps_lnrho = simple_on.dlneps_lnrho;
        nuc.dlneps_lnT = simple_on.dlneps_lnT;    

    }

    else {

        std::cerr << "Unknown nuc_output: " << nuc_output << "\n";
        return 1;
    }

    // inserting logger

    {
        std::ofstream log(filename, std::ios::app);

        // Enough precision to resolve small composition changes without
        // unnecessarily printing the full maximum decimal representation.
        log << std::setprecision(12);

        const int nr = T.nrows();
        const int nth = T.ncols();

        for (int i = 0; i < nr; ++i) {
        //for (int i = 0; i < std::min(nr, 180); ++i) {
            for (int j = 0; j < nth; ++j) {
                const int k = i * nth + j;
                const GridfireDiag& diag = diags[k];

                log << call_id << "," << i << "," << j << ","
                    << T(i,j) << "," << rho(i,j) << ","
                    << diag.num_steps << ","
                    << diag.raw_num_steps << ","
                    << static_cast<int>(diag.rhs_available) << ","
                    << diag.raw_netout_energy << ","
                    << diag.raw_netout_deps_drho << ","
                    << diag.raw_netout_deps_dT << ","
                    << diag.eps_inst << ","
                    << diag.dlneps_lnrho << "," << diag.dlneps_lnT << ","
                    << diag.deps_drho << "," << diag.deps_dT << ","
                    << nuc.eps(i,j) << ","
                    << nuc.dlneps_lnrho(i,j) << ","
                    << nuc.dlneps_lnT(i,j) << ","
                    << simple_cn.eps(i,j) << "," << simple_cn.dlneps_lnrho(i,j) << ","
                    << simple_cn.dlneps_lnT(i,j) << ","
                    << simple_on.eps(i,j) << "," << simple_on.dlneps_lnrho(i,j) << ","
                    << simple_on.dlneps_lnT(i,j) << ","
                    << diag.eps_avg << ","
                    << diag.Xcomp_H1 << "," << diag.Xcomp_He3 << ","
                    << diag.Xcomp_He4 << "," << diag.Xcomp_C12 << ","
                    << diag.Xcomp_N14 << "," << diag.Xcomp_O16 << ","
                    << diag.Xcomp_Ne20 << "," << diag.Xcomp_Mg24 << ","
                    << diag.Xin_H1 << "," << diag.Xin_He3 << ","
                    << diag.Xin_He4 << "," << diag.Xin_C12 << ","
                    << diag.Xin_N14 << "," << diag.Xin_O16 << ","
                    << diag.Xin_Ne20 << "," << diag.Xin_Mg24 << ","
                    << diag.Xout_H1 << "," << diag.Xout_He3 << ","
                    << diag.Xout_He4 << "," << diag.Xout_C12 << ","
                    << diag.Xout_N14 << "," << diag.Xout_O16 << ","
                    << diag.Xout_Ne20 << "," << diag.Xout_Mg24 << ","
                    << diag.evaluate_s << "," << diag.rhs_lookup_s << ","
                    << diag.postprocess_s << "\n";
            }
        }
        // No explicit flush here: std::ofstream flushes/closes on scope exit.
    }


    std::cout << "after log flush end of gridfire nuc run \n" << std::endl;

    return 0;

}
